## News and release notes for micom

For information on how to use `micom` please see the docs at
https://resendislab.github.io/micom.

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
