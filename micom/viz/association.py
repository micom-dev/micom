"""Visualization for phenotype prediction."""

from datetime import datetime
from micom.viz import Visualization
from micom.logger import logger
from micom.measures import production_rates, consumption_rates
from micom import stats
import json
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from sklearn.model_selection import (
    cross_val_predict,
    cross_val_score,
    LeaveOneOut,
)
from sklearn.linear_model import (
    LogisticRegressionCV,
    LassoCV,
    LogisticRegression,
    Lasso,
)
from sklearn.preprocessing import StandardScaler

PANDAS_VERSION = tuple(int(x) for x in pd.__version__.split("."))


def plot_association(
    results,
    phenotype,
    variable_type="binary",
    variable_name="phenotype",
    filename="association_%s.html" % datetime.now().strftime("%Y%m%d"),
    flux_type="production",
    fdr_threshold=0.05,
    threads=1,
    atol=1e-6,
):
    """Test for differential metabolite production.

    This will check for associations of the `phenotype` with metabolite fluxes. Individual
    tests are performed using non-parametric tests of the overall consumption or production
    fluxes for each samples versus the phenotype.

    To assess the the global association, this will fit L1-regularized linear models
    with log-fluxes as features. Will use LASSO regression for a continuous
    response and L1-regularized Logistic regression for a binary response.

    Parameters
    ----------
    results : micom.workflows.GrowthResults
        The results returned by the `grow` workflow.
    phenotype : pandas.Series
        The data to be fitted. Its index must correspond to `sample_id` in
        `exchanges`.
    variable_type : str of ["binary", "continuous"]
        The type of the variable.
    variable_name : str
        A short description of the phenotype for instance "disease status".
    filename : str
        The HTML file where the visualization will be saved.
    flux_type : str of ["import", "production"]
        Whether to fit using import or production fluxes.
    threads : int
        The number of threads to use.
    fdr_threshold : float
        The false discovery rate cutoff to use (FDR-corrected p-value cutoff). Defaults
        to 0.05.
    atol : float
        Tolerance to consider a flux different from zero. Will default
        to the solver tolerance.

    Returns
    -------
    Visualization
        A MICOM visualization. Can be served with `viz.view`.

    """
    exchanges = results.exchanges
    if flux_type == "import":
        exchanges = consumption_rates(results)
    else:
        exchanges = production_rates(results)
    exchanges = exchanges.loc[exchanges.flux > atol]
    if exchanges.shape[1] < 1:
        raise ValueError("None of the fluxes passed the tolerance threshold :(")
    if variable_type == "binary" and phenotype.nunique() != 2:
        raise ValueError(
            "Binary variables must have exactly two unique values, yours "
            "has: %s." % ", ".join(phenotype.unique())
        )
    elif variable_type == "continuous" and not is_numeric_dtype(phenotype):
        raise ValueError(
            "Continuous variables must have a numeric type, but yours is"
            " of type `%s`." % phenotype.dtype
        )
    elif variable_type not in ["binary", "continuous"]:
        raise ValueError(
            "Unsupported variable type. Must be either `binary` or " "`continuous`."
        )
    exchanges.loc[:, variable_name] = phenotype[exchanges.sample_id].values

    fluxes = exchanges.pivot_table(
        index="sample_id", columns="metabolite", values="flux", fill_value=atol
    )
    if PANDAS_VERSION >= (2, 1, 0):
        fluxes = fluxes.map(np.log)
    else:
        fluxes = fluxes.applymap(np.log)
    meta = phenotype[fluxes.index]
    stds = fluxes.std(axis=1)
    bad = stds < atol
    if bad.any():
        logger.warning("Removing %d fluxes due to zero variance." % bad.sum())
        fluxes = fluxes.loc[:, ~bad]
    scaled = StandardScaler().fit_transform(fluxes)
    if variable_type == "binary":
        model = LogisticRegressionCV(
            penalty="l1",
            scoring="accuracy",
            solver="liblinear",
            cv=2,
            Cs=np.power(10.0, np.arange(-6, 6, 0.5)),
            max_iter=50000,
        )
        fit = model.fit(scaled, meta)
        model = LogisticRegression(
            penalty="l1",
            solver="liblinear",
            C=fit.C_[0],
            max_iter=10000,
        )
        fit = model.fit(scaled, meta)
        score = cross_val_score(model, X=scaled, y=meta, cv=2)
        tests = stats.compare_groups(
            exchanges, metadata_column=variable_name, threads=threads, progress=False
        )
        statistic_name = "log fold-change"
        tests.rename(columns={"log_fold_change": "statistic"}, inplace=True)
    else:
        model = LassoCV(cv=2, max_iter=50000)
        fit = model.fit(scaled, meta)
        model = Lasso(alpha=fit.alpha_, max_iter=50000)
        fit = model.fit(scaled, meta)
        score = cross_val_score(model, X=scaled, y=meta, cv=2)
        tests = stats.correlate_fluxes(
            exchanges, metadata_column=variable_name, threads=threads, progress=False
        )
        statistic_name = "Spearman ρ"
        tests.rename(columns={"spearman_rho": "statistic"}, inplace=True)
    score = [np.mean(score), np.std(score)]
    score.append(model.score(scaled, meta))

    data = {"fluxes": exchanges, "tests": tests}
    significant = tests[tests.q < fdr_threshold]
    predicted = cross_val_predict(model, scaled, meta, cv=LeaveOneOut())
    fitted = pd.DataFrame({"real": meta, "predicted": predicted}, index=meta.index)

    exchanges = exchanges.loc[
        exchanges.metabolite.isin(significant.metabolite.values)
    ].copy()
    exchanges[variable_name] = meta[exchanges.sample_id].values
    var_type = "nominal" if variable_type == "binary" else "quantitative"
    viz = Visualization(filename, data, "tests.html")

    viz.save(
        fitted=fitted.to_json(orient="records"),
        tests=significant.to_json(orient="records"),
        exchanges=exchanges.to_json(orient="records"),
        metabolites=json.dumps([None] + significant.metabolite.tolist()),
        variable=variable_name,
        statistic=statistic_name,
        type=var_type,
        direction=flux_type,
        q_threshold=fdr_threshold,
        score=score,
        width=400,
        cwidth=max(12 * significant.shape[0], 160),
    )

    return viz
