# News and release notes for MICOM

This includes a list of major changes for each minor version starting from 0.19.0.

For information on how to use `micom` please refer to
[the documentation](https://micom-dev.github.io/micom).

### 0.37.1

Pins the OSQP version to a range compatible with optlang.

### 0.37.0

**Visualizations**

`plot_exchanges_per_taxon` and `plot_growth` now allow to color points by arbitrary
(categorical) metadata similar to `plot_mes`.

### 0.36.4

Fix a bug in `check_db_medium` where the growth rate order was random.

Do not overwrite the handlers for the root logger anymore but attach the handler to
a package logger instead.

### 0.36.3

Require a newer optlang version to stay compatible with numpy 2.0.

### 0.36.2

Limit logger default change to workflows only.

Make some attempts to clean the taxon name in interactions.

Workflows accepting per sample media will now raise an error if the medium is missing
samples.

### 0.36.1

Fixes a crash in the `minimal_media` workflow when growth rates are infeasible.

Changed the default log level to ERROR to avoid message spamming in the workflows. This
can be re-enabled with:

```python
from micom.logger import logger
logger.setLevel("WARNING")
```

### 0.36.0

**Minimal media workflow**

The minimal media workflow is now more flexible and supports all options of the
`minimal_media` base function including:

- specification of a required community growth rate, taxon growth rates, or any combination
- support for minimizing the number of media metabolites
- support for weights allowing to minimize mass uptake, carbon uptake etc.
- returning the growth results obtained by the media minimization

**Growth Results**

`GrowthResults` is now a DataClass and gains some new helpers:

- combining two GrowthResults via the addition operator (`r1 + r2`)
- concatenating many growth results via `combine_results`
- converting a `CommunitySolution` to a growth results object

**Other**

Renamed `fix_medium` to `complete_db_medium` for consistency with `complete_db_medium`
and to make it more clear what is happening.

Some black style fixes across the code base.

Updated optlang version for better numpy support.

Fix a bug where `minimal_medium(manifest, ..., summarize=False)` would crash.

### 0.35.0

`plot_fit` got renamed to `plot_association` and will not use LASSO coefficients anymore
as a proxy for univariate associations. Instead univariate non-parametric tests are
run for each metabolite. This should be more stable and avoid issues with low reproducibility
of coefficients. LASSO model performance is still reported as a global measure
of association. In case you need access to the previous `plot_fit` function you can simply
install a prior version of micom in a python or conda environment.

Visualizations are updates requiring less dependencies now and using a current version
of vega-lite.

The transition to the hybrid solver is now completed and the installation does not require
manual steps anymore.

Lowered the iteration limit for the OSQP step in the hybrid solver as more iterations
did usually not lead to an improvement or convergence.

Convert some warnings to debug messages.

Switch the the start method for multiprocessing to "spawn" which should improve
stability for workflows.

### 0.34.1

Fix a typo in the hybrid settings.

### 0.34.0

Better default settings for the hybrid solver.

Fixes an issue with taxonomy tables with species rank and without if column.

### 0.33.2

Fixes yet another annotation bug with AGORA2 due to duplicated exchanges.

### 0.33.1

Fixes the rank prefix trimming when setting it to False in Qiime2 imports.

### 0.33.0

Replaces the default OSQP solver with a new hybrid solver using OSQP for QPs and HIGHS
for everything else which improves the accuracy by several orders of magnitude. Note the new
[installation instructions for solvers](https://micom-dev.github.io/micom/installing.html).

The presence or absence of rank prefixes such as `s__` or `g__` for genus and species names,
respectively, can now be detected and resolved automatically by MICOM. This improves
compatibility with Qiime2 and GTDB.

New solver is now tested automatically on all platforms.

Small bugfixes and updates.

### 0.32.5

Backports the visualizers to be compatible with scikit-learn 0.24.1 to avoid some issues
with Qiime2.

Re-enables the Qiime model loading test in the CI.

### 0.32.4

Fixes a bug where model annotation tables could not be generated if there were multiple
exchange reactions for a single metabolite. This made MICOM not work with AGORA2 in some
instances.

Now checks for demand reactions and in database construction as well.

### 0.32.3

Bump the required cobra version to avoid numpy version errors.

### 0.32.1 - 0.32.2

Adjustments to solver settings and lost of improvements to Gurobi that make it much
faster.

### 0.32.0

Helper functions to save results were added.

### 0.31.7

Fixed a bug where minimal_medium would maximize with OSQP instead of minimize.

### 0.31.6

More adjustments to OSQP.

### 0.31.5

Fix a bug in the stabilization term used for OSQP. Try to avoid fragmentation warning.

### 0.31.4

Fix regression in T-SNE.

### 0.31.2 - 0.31.3

Some more fixes to the media completion.

### 0.31.1

Show added flux as fraction of total in `complete_db_medium`.

### 0.31.0

The medium completion functions now have a flag that allows adding additional flux
to components already in the candidate medium. The default is now changed to that
new non-strict mode which makes it much easier and more efficient to complete media.
In that case the `max_added_import` amount now specifies the maximum flux added on top
the limit specified in the candidate medium.

### 0.30.5

Fixed a bug where T-SNE would complain in low sample numbers.

### 0.30.3 - 0.30.4

Those were fixes to deployment and dependencies only.

### 0.30.2

Avoids a bug in COBRAPY that can yield to huge serialized versions of the models.

### 0.30.1

Fix compatibility with Python 3.7.

### 0.30.0

`build` will now detect existing models in the output folder and skip them. This now
allows to resume interrupted builds. Explicit rebuild requires you to delete the
folder now.

### 0.29.6

Fixes a bug where using more than 1 thread for workflows would not return results
in some cases.

### 0.29.5

Fixed some dependency errors on older Python versions.

### 0.29.4

Port tests to cobrapy 0.25 and above. Fix an import error in workflows.core.

### 0.29.3

Suppresses some warning in annotations from cobrapy.

### 0.29.2

Fixes some CI and build issues.

### 0.29.0

Renames the `n_jobs` arguments to `threads` across the code base for consistency and
updates the arguments for `workflow`.

`workflow` now does not spawn a process for single-threaded runs, reducing overhead.

Fixes an issue where `workflow` would deadlock with pytest-cov.

### 0.28.1

Fixes an issue where an accumulation in warnings in spawned processes trips up
jupyter.

### 0.28.0

Fixes a deployment issue with the previous build.

### 0.27.0

Annotations now include the molecular weight, number of carbon atoms, and the number of
nitrogen atoms for each compound.

### 0.26.0

Now adjusts active demands to be optional. This was an issue with some models
in AGORA that had forced biotin demands. This effectively prohibits zero growth
rates and is problematic for MICOM.

Fix the logger for multiprocessing.

### 0.25.1

Handle CARVEME models better in `workflows.complete_medium`.

Workflows now have a deterministic return error (same order as inputs).

### 0.25.0

`build_database` now allows to set the compression algorithm and level.

### 0.24.1

Fixed settings to get better convergence behavior with OSQP.

### 0.24.0

MICOM now works with OSQP and works out of the box without a commercial QP solver. See
the [installation docs]() for more information.

### 0.23.1

The sample heatmap in `plot_exchanges_per_sample` is now automatically rotated when
there are mores samples than reactions.

media workflows will use presolving by default since those are often numerically
problematic.

### 0.23.0

Fix the signature and add deprecation warnings for optimize_* methods.

`plot_exchanges_per_taxon` will now use normalized fluxes (per 1 gDW for each taxon)
which allows better comparability. The old behavior can be enabled with
`use_total_flux=True`.

Avoid negative growth rate constraints in `_apply_min_growth`.

Can now enable presolving/scaling in `grow` and `tradeoff`.

### 0.22.7

Fix some warnings from pandas.

Avoid a crash in `reularize_l2_norm` when all reactions for a taxon have been fixed to
zero.

Raise a better error if the minimal medium optimization fails in grow.

Use the right tolerance when setting atol and rtol automatically.

### 0.22.6

Fixed an error in `micom.workflows.build` if build was run without a model database
but with a `file` column in the taxonomy.

### 0.22.5

Fixed missing data files.

### 0.22.4

Improves version compatibility with cobrapy and QIIME 2.

### 0.22.3

Fixed a bug where incorrectly labeled biomass demand reactions were treated like an
exchange.

### 0.22.2

Lowered the required version for pandas to re-enable compatibility with Python 3.6.

Docs are now built on all pushes.

### 0.22.1

`atol` and `rtol` are now consistently exposed.

Remove six dependency.

### 0.22.0

Got a bit smarter in cleaning up compartment suffixes. This fixes the odd "_e_m" suffix
in CARVEME-derived models. This will change the names of exchange reactions compared
to version 0.21.x.

### 0.21.0

Stabilize minimal_medium a bit more.

The strategy to calculate fluxes in the grow workflow can now be set with the
`strategy` argument. For instance, you can now get fluxes via parsimonious FBA instead of
assuming minimal imports.

Fixed a bug where minimal media classification with `weights="mass"` would fail due to
invalid elements in the formula (for instance "X").

### 0.20.0

This version brings new functionality to design growth media from a skeleton medium.
This also allows for a quicker verification of media against the model databases.

Added workflows:

- `check_db_media` check if models in a model database can grow on a given medium
- `complete_db_media`completes a growth medium so all models in the db can grow on it

Together with this we now provide several new revised growth media for the gut on
github.com/micom-dev/media. In particular, we finally provide the often requested growth
medium for the carveME database.

### 0.19.0

`minimal_medium` now accepts weighting the fluxes bei molecular weight or
any elemental content. For instance, you can now calculate a minimal medium that
minimizes carbon import for instance. The used `weigths` argument propagates
to any workflow using this function including `complete_medium`, `fix_medium`
and `grow`.

### 0.1.0 - 0.18.0

`minimal_medium` now has an option to return all fluxes along with the
import and export fluxes. Useful if you want to check what every individual
taxa consumes.

Addition of `Community.cooperative_tradeoff` which brings a fast method to
get nearly unique growth rates for a particular growth tradeoff (for instance
50% maximal growth rate).

Addition of elasticity coefficients which allows you to study which import
and export reactions are sensitive to changes in diet or taxa abundances.

Several changes to make `micom` capable of handling numerical instable models.
This includes a crossover strategy to obtain an optimal or near-optimal
from an incomplete interior point optimization.
