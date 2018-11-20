
import pandas as pd

from .plugin_setup import plugin
from ._format import AbundanceFilteringStatsFmt


@plugin.register_transformer
def _1(data: pd.DataFrame) -> AbundanceFilteringStatsFmt:
    stats_obj = AbundanceFilteringStatsFmt()
    stats_obj_filepath_str = str(stats_obj)
    data.to_csv(stats_obj_filepath_str)
    return stats_obj

@plugin.register_transformer
def _2(stats_obj: AbundanceFilteringStatsFmt) -> pd.DataFrame:
    stats_obj_filepath_str = str(stats_obj)
    return pd.read_csv(stats_obj_filepath_str, index_col='sample-id')

