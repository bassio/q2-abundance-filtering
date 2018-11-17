
import qiime2.plugin
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    Sequences, SequencesWithQuality, PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality)

citations = qiime2.plugin.Citations.load('citations.bib', package='q2_abundance_filtering')


import q2_abundance_filtering

plugin = qiime2.plugin.Plugin(
    name='abundance-filtering',
    version=q2_abundance_filtering.__version__,
    website='https://github.com/bassio/q2-abundance-filtering',
    package='q2_abundance_filtering',
    user_support_text=None,
    short_description='Plugin for abundance filtration of sequences according to the method of Wang et al.',
    description=('This plugin filters out sequences according to the method of Wang et al.'),
    citations=[citations['Wang2018']]
)

plugin.methods.register_function(
    function=q2_abundance_filtering.abundance_filter,
    inputs={
        'sequences': SampleData[JoinedSequencesWithQuality | SequencesWithQuality]
        }
    ,
    parameters={
    },
    outputs=[
        ('filtered_sequences', SampleData[SequencesWithQuality])
    ],
    input_descriptions={
        'sequences': "The input sequences."
    },
    parameter_descriptions={
    },
    output_descriptions={
        'filtered_sequences': 'The filtered sequences.',
    },
    name='Abundance filtering',
    description=('Filter low-abundance sequences according to the method of Wang et al. Requires the dereplicated sequences as well as the dereplicated sequence counts table. '
                 'Outputs the filtered sequences.')
)


