## News and release notes for micom

This includes a list of major changes for each minor version starting from 0.19.0.

For information on how to use `micom` please see the docs at
https://micom-dev.github.io/micom.

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

### 0.9.0

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
