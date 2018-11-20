import pandas as pd


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


def abundance_filter_threshold_Wang_et_al(counts_series):
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
 
