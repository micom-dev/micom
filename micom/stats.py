"""Data warngling and statistics for MICOM."""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal, spearmanr, pearsonr
from .workflows import workflow


def fdr_adjust(p):
    """Get FDR cutoffs for p-values with Benjamini-Hochberg.

    Arguments
    ---------
    p : list[float]
        The original p-values. Can not contain naNs.

    Returns
    -------
    A numpy array of FDR cutoffs (q-values). This is commonly known as "adjusted"
    p-values.
    """
    p = np.array(p)
    if np.isnan(p).any():
        raise ValueError("`p` can not contain NAs. Please filter those before.")
    order = np.argsort(-p)
    reverse_order = np.argsort(order)
    n = len(p)
    i = np.arange(n, 0, -1)
    q = np.minimum.accumulate(n / i * p[order])[reverse_order]
    q[q > 1.0] = 1.0
    return q


def _run_test(args):
    df, col, groups = args
    met = df.metabolite.iloc[0]
    if len(groups) == 2:
        lfc = np.log2(df.flux[df[col] == groups[1]].mean() + 1e-6) - np.log2(
            df.flux[df[col] == groups[0]].mean() + 1e-6
        )
        try:
            p = mannwhitneyu(
                df.flux[df[col] == groups[0]], df.flux[df[col] == groups[1]]
            )[1]
        except Exception:
            p = np.nan

        return pd.DataFrame(
            {
                "metabolite": met,
                "log_fold_change": lfc,
                "p": p,
                "n": df.shape[0],
            },
            index=[0],
        )
    elif len(groups) > 2:
        lstd = df.groupby(col).flux.apply(lambda x: np.log2(x.mean() + 1e-6)).std()
        try:
            gs = [df.flux[df[col] == c] for c in groups]
            p = kruskal(*gs)[1]
        except Exception:
            p = np.nan

        return pd.DataFrame(
            {
                "metabolite": met,
                "log_mean_std": lstd,
                "p": p,
                "n": df.shape[0],
            },
            index=[0],
        )


def _run_corr(args):
    df, col = args
    met = df.metabolite.iloc[0]
    lpr = pearsonr(np.log2(df["flux"] + 1e-6), df[col])[0]
    try:
        p = spearmanr(df["flux"], df["time"])[1]
    except Exception:
        p = np.nan
    return pd.DataFrame(
        {
            "metabolite": met,
            "log_pearson_rho": lpr,
            "p": p,
            "n": df.shape[0],
        },
        index=[0],
    )


def compare_groups(fluxes, metadata_column, groups=None, threads=1, progress=True):
    """Compare fluxes form different sample groups.

    Note
    ----
    This uses a non-parametric test by default. By default it will use a
    Mann-Whitney test for two groups and a Kruskal-Wallis test for >2 groups.

    Arguments
    ---------
    fluxes : pandas.DataFrame
        A frame with net fluxes as returned by `production_rates` or
        `consumption_rates`.
    metatdata_column : str
        The column of the DataFrame denoting the groups.
    groups : list[str] or None
        Specify a subset of groups you want to compare or define the order (1st will
        be the reference group). If None will use the groups as they appear in
        the DataFrame.
    threads : int
        How many threads to use to run tests in parallel.
    progress : bool
        Whether to show a progress bar.

    Returns
    -------
    Returns the metabolite with their respective test statistics.
    """
    if groups is not None:
        fluxes = fluxes[fluxes[metadata_column].isin(groups)]
    else:
        groups = np.unique(fluxes[metadata_column])

    met_data = [
        (fluxes[fluxes.metabolite == m], metadata_column, groups)
        for m in fluxes.metabolite.unique()
    ]
    res = workflow(_run_test, met_data, threads=threads, progress=progress)
    res = pd.concat(res)
    res["q"] = np.nan
    res.loc[res.p.notnull(), "q"] = fdr_adjust(res.loc[res.p.notnull(), "p"])
    return res


def correlate_fluxes(fluxes, metadata_column, groups=None, threads=1, progress=True):
    """Compare fluxes form different sample groups.

    Note
    ----
    This uses a non-parametric test by default (Spearman rank correlation).

    Arguments
    ---------
    fluxes : pandas.DataFrame
        A frame with net fluxes as returned by `production_rates` or
        `consumption_rates`.
    metatdata_column : str
        The column of the DataFrame denoting the covariate.
    threads : int
        How many threads to use to run tests in parallel.
    progress : bool
        Whether to show a progress bar.

    Returns
    -------
    Returns the metabolite with their respective test statistics.
    """
    met_data = [
        (fluxes[fluxes.metabolite == m], metadata_column)
        for m in fluxes.metabolite.unique()
    ]
    res = workflow(_run_corr, met_data, threads=threads, progress=progress)
    res = pd.concat(res)
    res["q"] = np.nan
    res.loc[res.p.notnull(), "q"] = fdr_adjust(res.loc[res.p.notnull(), "p"])
    return res
