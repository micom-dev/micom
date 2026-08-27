"""Type definitions and validations for the micom library."""

from functools import wraps
from pathlib import Path
from os import PathLike
import inspect
import logging
import pandas as pd
import typing
import types

from .constants import RANKS

logger = logging.getLogger(__name__)


def pathify(func):
    """Handle arguments that are hinted as Path or PathLike (including Unions) and convert them to pathlib.Path objects.

    Decorator that automatically converts arguments hinted as Path or PathLike
    (including Unions) into pathlib.Path objects before calling the function.

    Parameters
    ----------
    func : function
        The function to be decorated.

    Returns
    -------
    function
        The decorated function that converts Path or PathLike arguments to pathlib.Path.

    """

    sig = inspect.signature(func)
    type_hints = typing.get_type_hints(func)

    def is_path_hint(hint):
        if hint is None:
            return False
        origin = typing.get_origin(hint)
        if origin is typing.Union or origin is types.UnionType:
            return any(is_path_hint(arg) for arg in typing.get_args(hint))
        try:
            return hint is Path or issubclass(hint, PathLike)
        except TypeError:
            return False

    @wraps(func)
    def wrapper(*args, **kwargs):
        bound_args = sig.bind(*args, **kwargs)
        bound_args.apply_defaults()

        for name, value in bound_args.arguments.items():
            hint = type_hints.get(name)
            if hint and is_path_hint(hint):
                if isinstance(value, (str, bytes, PathLike)):
                    bound_args.arguments[name] = Path(value)

        return func(*bound_args.args, **bound_args.kwargs)

    return wrapper


def check_medium(medium):
    """Validate a medium for simulation.


    Parameters
    ----------
    medium : pd.DataFrame
        The medium to validate.

    Raises
    ------
    ValueError
        If the medium is invalid.

    Returns
    -------
    None

    """
    if "reaction" not in medium.columns:
        raise ValueError("The medium must have a 'reaction' column.")
    if "flux" not in medium.columns:
        raise ValueError("The medium must have a 'flux' column.")
    if "sample_id" in medium.columns:
        cn = medium.dropna().groupby("sample_id").reaction.value_counts()
        if any(cn > 1):
            raise ValueError("The medium contains duplicate reactions per sample.")
        if any(cn == 0):
            raise ValueError("The medium contains samples with no reactions.")
    else:
        cn = medium.dropna().reaction.value_counts()
        if any(cn < 1):
            raise ValueError(
                "The medium contains duplicate reactions and no 'sample_id' column."
            )
        if any(cn == 0):
            raise ValueError(
                "The medium contains no reactions and no 'sample_id' column."
            )


def check_taxonomy(taxonomy: pd.DataFrame, samples: bool = True) -> None:
    """Check if the taxonomy is valid.

    Parameters
    ----------
    taxonomy : pd.DataFrame
        The taxonomy of the community models.
    samples : bool
        Whether to allow multiple samples in the taxonomy. If False, will raise an error if multiple samples are found. Default is True.

    Raises
    ------
    ValueError
        If the taxonomy has issues.

    Returns
    -------
    Nothing.

    """
    found = taxonomy.columns.isin(["id"] + RANKS)
    if not "abundance" in taxonomy.columns:
        raise ValueError(f"Taxonomy must contain columns 'abundance' and 'sample_id'")
    if found.sum() == 0:
        raise ValueError(
            f"Taxonomy must contain at least one column from: {", ".join(RANKS)}"
        )

    if samples:
        if "sample_id" not in taxonomy.columns:
            raise ValueError(
                f"Taxonomy must contain a 'sample_id' column in multi-sample setups (like using `Batch()`)."
            )
    elif "sample_id" in taxonomy.columns and taxonomy.sample_id.nunique() > 1:
        raise ValueError(
            f"Taxonomy can only contain a single 'sample_id' in single-sample setups."
        )
    else:
        taxonomy = taxonomy.copy()
        taxonomy["sample_id"] = "sample 1"

    if "id" in taxonomy.columns:
        if taxonomy.id.str.contains(r"[^a-zA-Z0-9_]").any():
            raise ValueError(
                "Taxonomy 'id' column contains invalid characters. Only alphanumeric characters and underscores are allowed."
            )
        lowest_rank = "id"
    else:
        lowest_rank = [r for r in RANKS if r in taxonomy.columns][-1]

    # Check for some common mistakes

    # Check for duplicate entries for single samples
    lowest_rank_counts = taxonomy.groupby("sample_id")[lowest_rank].value_counts()
    if (lowest_rank_counts > 1).any():
        raise ValueError(
            f"Found duplicate entries for single samples for '{lowest_rank}' in the taxonomy."
            " Each sample should have only one collapsed abundance for the lowest rank/ID."
            " Please check your taxonomy file."
        )

    # Check if each lowest rank appears only once in the taxonomy
    lowest_rank_counts = taxonomy[lowest_rank].value_counts()
    s_counts = taxonomy.sample_id.nunique()
    if (s_counts > 1) & (lowest_rank_counts == 1).all():
        logger.warning(
            f"Each '{lowest_rank}' appears only once in the taxonomy."
            " Note that taxa names and IDs should be unique for each *organism* in the community, not for each sample."
            " This might be okay if your samples do not share any taxa, but it is worth checking your taxonomy file."
        )

    # Check for zero abundance samples
    zero_abundance_samples = taxonomy.groupby("sample_id")["abundance"].sum()
    zero_abundance_samples = zero_abundance_samples[zero_abundance_samples == 0]
    if len(zero_abundance_samples) > 0:
        raise ValueError(
            f"Found {len(zero_abundance_samples)} samples with zero abundances: "
            ", ".join(zero_abundance_samples.index)
        )
