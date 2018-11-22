import gzip
import hashlib
from io import StringIO
from collections import Counter
from multiprocessing import Pool
import pandas as pd
import biom
import skbio
import yaml


from q2_types.per_sample_sequences import (
            SingleLanePerSampleSingleEndFastqDirFmt,
            FastqManifestFormat, YamlFormat, FastqGzFormat, QIIME1DemuxDirFmt
            )


from ._abundance_filter import abundance_filter_threshold_Wang_et_al


STATS_COLUMNS = ['sample-id', 'threshold', 'n_input_seqs', 'n_seqs_kept', 'n_seqs_unique_kept']


def return_sample_ids(demux: SingleLanePerSampleSingleEndFastqDirFmt):
    demux_manifest = demux.manifest.view(demux.manifest.format)
    
    with demux_manifest.open() as f:
        manifest_df = pd.read_csv(f, dtype=str, comment='#')
    
    return list(manifest_df['sample-id'])


def return_fastqgz_path_for_sample(demux: SingleLanePerSampleSingleEndFastqDirFmt,
                                   sample_id: str):
    #barcode ID, lane number and read number are not relevant here
    return demux.sequences.path_maker(sample_id=sample_id,
                                        barcode_id=1,
                                        lane_number=1,
                                        read_number=1)


def return_fastq_seqs_for_sample(demux: SingleLanePerSampleSingleEndFastqDirFmt,
                            sample_id: str):
    
    demux_manifest = demux.manifest.view(demux.manifest.format)
    
    with demux_manifest.open() as f:
        manifest_df = pd.read_csv(f, dtype=str, comment='#')
        manifest_df.set_index('sample-id', inplace=True)
    
    filename = manifest_df.loc[str(sample_id)]['filename']
    
    fastqz_pth = demux.path / filename
    
    seqs_iterator = skbio.io.read(str(fastqz_pth), format='fastq')
    
    return seqs_iterator


def sample_counts_series(demux: SingleLanePerSampleSingleEndFastqDirFmt,
                     sample_id: str) -> pd.Series:
    
    c = Counter([str(seq) for seq in return_fastq_seqs_for_sample(demux, sample_id)])
    
    return pd.Series(c, name=sample_id).sort_values(ascending=False)


def abundance_filter_sample(demux_path_str: SingleLanePerSampleSingleEndFastqDirFmt,
                               sample_id: str,
                               new_demux_path_str: SingleLanePerSampleSingleEndFastqDirFmt):
    
    print("Commencing abundance-filtering for sample {}".format(sample_id))

    demux = SingleLanePerSampleSingleEndFastqDirFmt(demux_path_str, mode='r')
    new_demux = SingleLanePerSampleSingleEndFastqDirFmt(new_demux_path_str, mode='r')
            
    
    stats_dict = {}
    stats_dict['sample-id'] = sample_id
    
    counts_series = sample_counts_series(demux, sample_id)
    
    stats_dict['n_input_seqs'] = counts_series.sum()    
    
    threshold = abundance_filter_threshold_Wang_et_al(counts_series)
    stats_dict['threshold'] = threshold
    
    abundance_filtered_series = counts_series[counts_series > threshold].copy()
    
    n_seqs_kept = abundance_filtered_series.sum()
    stats_dict['n_seqs_kept'] = n_seqs_kept
    
    n_seqs_unique_kept = len(abundance_filtered_series.index)
    stats_dict['n_seqs_unique_kept'] = n_seqs_unique_kept
    
    demux_metadata_view = demux.metadata.view(YamlFormat)
    with open(str(demux_metadata_view)) as demux_metadata_fh:
        demux_metadata_dict = yaml.load(demux_metadata_fh)
        phred_offset = demux_metadata_dict['phred-offset']
    
    new_fastqgz_path = new_demux.sequences.path_maker(sample_id=sample_id,
                                                      barcode_id=1,
                                                      lane_number=1,
                                                      read_number=1)
    
    
    with gzip.open(str(new_fastqgz_path), mode='w') as writer:
        for seq in return_fastq_seqs_for_sample(demux, sample_id):
            if str(seq) in abundance_filtered_series.index:
                seq_strIO = StringIO()
                skbio.io.write(seq, format='fastq', phred_offset=phred_offset, into=seq_strIO)
                
                writer.write(bytes(seq_strIO.read(), encoding="UTF-8"))
        
        
    return str(new_fastqgz_path), stats_dict



def finalize_result(original_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                    result_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                    sample_fastqgz_mapping: dict,
                    stats_df: pd.DataFrame) -> SingleLanePerSampleSingleEndFastqDirFmt:
    
    print("in finalize result")
    
    demux = original_sequences
    
    result = result_sequences
    
    #exclude those sample with resulting zero sequences after filtering
    #since a fastq file with zero sequences is invalid
    sample_ids_to_include = list(stats_df[stats_df['n_seqs_kept'] > 0]['sample-id'])
    
    #exit with error here
    if len(sample_ids_to_include) == 0:
        raise ValueError("All sequences from all samples were filtered out through abundance-filtering.")
    
    
    #manifest
    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')
    manifest_fh.write('# direction is not meaningful in this file as these\n')
    manifest_fh.write('# data may be derived from forward, reverse, or \n')
    manifest_fh.write('# joined reads\n')
    
    for sample_id in sample_ids_to_include:
        path = return_fastqgz_path_for_sample(demux, sample_id=sample_id)
        manifest_fh.write('{sample_id},{filename},{direction}\n'.format(sample_id=sample_id,
                                                               filename=path.name,
                                                               direction='forward'))
    
    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)
    
    ###metadata
    demux_metadata_view = demux.metadata.view(YamlFormat)
    with open(str(demux_metadata_view)) as demux_metadata_fh:
        demux_metadata_dict = yaml.load(demux_metadata_fh)
    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump(demux_metadata_dict))
    result.metadata.write_data(metadata, YamlFormat)
    
    
    return result

    
def abundance_filter_single_thread(sequences:SingleLanePerSampleSingleEndFastqDirFmt) -> (SingleLanePerSampleSingleEndFastqDirFmt, 
                                                                            pd.DataFrame):
    list_of_stats_dicts = []
    
    sample_fastqgz_mapping = {}
    
    result = SingleLanePerSampleSingleEndFastqDirFmt()
    
    sample_ids = return_sample_ids(sequences)
    
    for sample_id in sample_ids:
        sample_result = abundance_filter_sample(demux_path_str=str(sequences),
                                                   sample_id=sample_id,
                                                   new_demux_path_str=str(result))
        
        fastqgz_pth_str, stats_dict = sample_result[0], sample_result[1]
        fastqgz = FastqGzFormat(fastqgz_pth_str, mode='r')
            
        sample_fastqgz_mapping[stats_dict['sample-id']] = fastqgz
        list_of_stats_dicts.append(stats_dict)
        
    
    stats_df = pd.DataFrame(list_of_stats_dicts)
    stats_df = stats_df[STATS_COLUMNS]
    
    abundance_filtered_result = finalize_result(sequences,
                                                result,
                                                sample_fastqgz_mapping,
                                                stats_df)
    
    return abundance_filtered_result, stats_df


def abundance_filter_pool(sequences:SingleLanePerSampleSingleEndFastqDirFmt,
                          threads:int) -> (SingleLanePerSampleSingleEndFastqDirFmt, 
                                             pd.DataFrame):
    
    list_of_stats_dicts = []
    
    sample_fastqgz_mapping = {}
    
    result = SingleLanePerSampleSingleEndFastqDirFmt()
    
    sample_ids = return_sample_ids(sequences)
    

    with Pool(processes=threads) as pool:
        iterable = [(str(sequences), sample_id, str(result))
                    for sample_id in sample_ids]
                    
        sample_results = pool.starmap(abundance_filter_sample,
                        iterable)
        
        for fastqgz_pth_str, stats_dict in sample_results:
            fastqgz = FastqGzFormat(fastqgz_pth_str, mode='r')
            
            sample_fastqgz_mapping[stats_dict['sample-id']] = fastqgz
            list_of_stats_dicts.append(stats_dict)
            
    
    stats_df = pd.DataFrame(list_of_stats_dicts)
    stats_df = stats_df[STATS_COLUMNS]
    
    abundance_filtered_result = finalize_result(sequences,
                                                result,
                                                sample_fastqgz_mapping,
                                                stats_df)
    
    return abundance_filtered_result, stats_df


def abundance_filter(sequences:SingleLanePerSampleSingleEndFastqDirFmt,
                          threads:int = 1) -> (SingleLanePerSampleSingleEndFastqDirFmt, 
                                             pd.DataFrame):
    if threads == 1:
        return abundance_filter_single_thread(sequences)
    else:
        return abundance_filter_pool(sequences, threads)

