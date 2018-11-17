import gzip
import hashlib
from io import StringIO

import pandas as pd
import biom
import skbio
import yaml

from q2_types.per_sample_sequences import (
            SingleLanePerSampleSingleEndFastqDirFmt,
            FastqManifestFormat, YamlFormat, FastqGzFormat, QIIME1DemuxDirFmt)



def biomtable_to_dataframe(biom_table_object):
  _bt = biom_table_object
  data = _bt.matrix_data.todense()
  df = pd.DataFrame(data, index=_bt.ids('observation'),
                           columns=_bt.ids('sample'))
  return df


def _rarefy_series(series, n, num_reps=1):

    series_name = series.name

    if n == 0:
        sampled_series = pd.Series(0, index=series.index, name=series.name)
        return sampled_series

    sampled_series = series.sample(n, replace=True, weights=series)
    sampled_value_counts = sampled_series.index.value_counts()
    sampled_value_counts.name = series_name

    if num_reps == 1:
        return sampled_value_counts
    elif num_reps > 1:
        for r in range(num_reps - 1):
            next_sampled_series = series.sample(n, replace=True, weights=series)
            next_sampled_value_counts = next_sampled_series.index.value_counts()

            sampled_value_counts = pd.concat([sampled_value_counts, next_sampled_value_counts], axis=1).fillna(0).sum(axis=1)
            sampled_value_counts.name = series_name
            print(sampled_value_counts.sort_values(ascending=False))

            del next_sampled_series
            del next_sampled_value_counts

        sampled_value_counts = sampled_value_counts.sample(n, replace=True, weights=sampled_value_counts).index.value_counts()
        sampled_value_counts.name = series_name

        return sampled_value_counts


def _bootstrap_series_concatenate(series, num_reps=1):
    series_sum = series.sum()
    
    series_list = [_rarefy_series(series, series_sum, num_reps=1) 
                   for i in range(num_reps)]
    
    return pd.concat(series_list, axis=1).fillna(0)


def abundance_filter_threshold_Wang_et_al(counts_df):
    counts_series = counts_df[counts_df.columns[0]].sort_values(ascending=False)

    bootstrap_df = _bootstrap_series_concatenate(counts_series, num_reps=1000)
    bootstrap_df_transposed = bootstrap_df.transpose()


    abund_real = counts_series.copy()

    abund_boot = bootstrap_df_transposed.mean().sort_values(ascending=False)
    abund_995 = bootstrap_df_transposed.quantile(0.995)
    abund_005 = bootstrap_df_transposed.quantile(0.005)

    abund_adj = (2 * abund_real) - abund_boot
    abund_adj = abund_adj.sort_values(ascending=False).fillna(0)
    
    ci99_higher = abund_adj + (abund_995 - abund_boot)
    ci99_higher = ci99_higher.sort_values(ascending=False).fillna(0)

    ci99_lower = abund_adj - (abund_boot - abund_005)
    ci99_lower = ci99_lower.sort_values(ascending=False)
    
    unreliable = ci99_lower[ci99_lower <= 0].index

    threshold = int(counts_series[unreliable].sort_values(ascending=False)[0])
    
    return threshold


#adapted with modifications from q2-vsearch and q2-quality-filter
def abundance_filter_seqs(demux: SingleLanePerSampleSingleEndFastqDirFmt,
                          filtered_counts_df: pd.DataFrame) -> SingleLanePerSampleSingleEndFastqDirFmt:
    
    result = SingleLanePerSampleSingleEndFastqDirFmt()

    manifest = FastqManifestFormat()
    manifest_fh = manifest.open()
    manifest_fh.write('sample-id,filename,direction\n')
    manifest_fh.write('# direction is not meaningful in this file as these\n')
    manifest_fh.write('# data may be derived from forward, reverse, or \n')
    manifest_fh.write('# joined reads\n')

    log_records_totalread_counts = {}

    demux_metadata_view = demux.metadata.view(YamlFormat)
    
    with open(str(demux_metadata_view)) as demux_metadata_fh:
        demux_metadata_dict = yaml.load(demux_metadata_fh)
        phred_offset = demux_metadata_dict['phred-offset']
    
    demux_manifest = demux.manifest.view(demux.manifest.format)
    demux_manifest = pd.read_csv(demux_manifest.open(), dtype=str, comment='#')
    demux_manifest.set_index('filename', inplace=True)
    
    iterator = demux.sequences.iter_views(FastqGzFormat)
    
    for i, (fname, fastqz) in enumerate(iterator):
        sample_id = demux_manifest.loc[str(fname)]['sample-id']

        # per q2-demux, barcode ID, lane number and read number are not
        # relevant here
        path = result.sequences.path_maker(sample_id=sample_id,
                                           barcode_id=i,
                                           lane_number=1,
                                           read_number=1)

        
        
        # per q2-demux, barcode ID, lane number and read number are not
        # relevant here
        path = result.sequences.path_maker(sample_id=sample_id,
                                           barcode_id=i,
                                           lane_number=1,
                                           read_number=1)

        # we do not open a writer by default in the event that all sequences
        # for a sample are filtered out; an empty fastq file is not a valid
        # fastq file.
        writer = None
        
        seqs_iterator = skbio.io.read(str(fastqz), format='fastq')
        
        sampled_counts_series = filtered_counts_df[sample_id]
        seqs_in_sample = list(sampled_counts_series[sampled_counts_series > 0].index)
        
        log_records_totalread_counts[sample_id] = 0

        for seq in seqs_iterator:
            sha1_hash = hashlib.sha1(str(seq).encode('utf-8')).hexdigest()
            if sha1_hash in seqs_in_sample:
                log_records_totalread_counts[sample_id] += 1
                if writer is None:
                    writer = gzip.open(str(path), mode='w')
            
                seq_strIO = StringIO()
                skbio.io.write(seq, format='fastq', phred_offset=phred_offset, into=seq_strIO)
                seq_string = seq_strIO.getvalue()

                writer.write(seq_string.encode('utf-8'))
        

        if writer is not None:
            manifest_fh.write('%s,%s,%s\n' % (sample_id, path.name, 'forward'))
            
            writer.close()
            

    if set(log_records_totalread_counts.values()) == {0, }:
        raise ValueError("All sequences from all samples were filtered out. "
                         "The parameter choices may be too stringent for the "
                         "data.")

    manifest_fh.close()
    result.manifest.write_data(manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump(demux_metadata_dict))
    result.metadata.write_data(metadata, YamlFormat)


    return result



def abundance_filter(sequences:SingleLanePerSampleSingleEndFastqDirFmt) -> SingleLanePerSampleSingleEndFastqDirFmt:
    
    from qiime2 import Artifact
    from qiime2.plugins.vsearch.methods import dereplicate_sequences

    sequences_with_type = Artifact.import_data("SampleData[SequencesWithQuality]", sequences)
    
    from q2_vsearch._cluster_sequences import dereplicate_sequences as plugin_dereplicate_sequences
    
    dereplicated_table, derep_sequences = plugin_dereplicate_sequences(sequences_with_type.view(QIIME1DemuxDirFmt))

    counts_df = biomtable_to_dataframe(dereplicated_table).astype(int)
    
    
    filtered_df = pd.DataFrame()
    
    for sample in counts_df.columns:
        
        sample_col = counts_df[sample]
        
        threshold = abundance_filter_threshold_Wang_et_al(sample_col.to_frame())
    
        filtered_df[sample] = sample_col[sample_col > threshold]
        
        seq_ids_to_keep = sample_col[sample_col > threshold].index
        
        n_kept = len(list(seq_ids_to_keep))
        
        print("Abundance filtering threshold for sample {} set at > {} counts. {} *unique* sequences kept.".format(sample,
                                                                                                                   str(threshold),
                                                                                                                   str(n_kept)))
    
    filtered_df = filtered_df.fillna(0).astype(int)
    
    #exclude these sequences that were excluded from all samples
    sum_series = filtered_df.sum(axis=1)
    zero_sum_sequences = sum_series[sum_series > 0].index
    filtered_df = filtered_df.reindex(zero_sum_sequences)
    

    af_seqs = abundance_filter_seqs(demux=sequences,
                                    filtered_counts_df=filtered_df)
    
    return af_seqs


