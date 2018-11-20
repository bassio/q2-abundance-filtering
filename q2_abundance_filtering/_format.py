import qiime2.plugin.model as model


STATS_HEADER = ['sample-id', 'threshold', 'n_input_seqs', 'n_seqs_kept', 'n_seqs_unique_kept']

STATS_DESCRIPTIONS = {
    'sample-id': "The sample ID",
    'threshold': "The threshold for abundance-filtering taken for each sample (include > n_threshold reads)",
    'n_input_seqs': "The total number of input sequence presented before abundance-filteringunique reads following dereplicaton",
    'n_seqs_kept': "The total number of sequences kept in fastq file following abundance-filtering",
    'n_seqs_unique_kept': "The number of unique  sequences (i.e. dereplicated 100% similarity) kept following abundance-filtering",
    }


class AbundanceFilteringStatsFmt(model.TextFileFormat):
    def sniff(self):
        line = open(str(self)).readline()
        hdr = line.strip().split(',')
        return hdr == STATS_HEADER


AbundanceFilteringStatsDirFmt = model.SingleFileDirectoryFormat(
'AbundanceFilteringStatsDirFmt', 'stats.csv', AbundanceFilteringStatsFmt) 
