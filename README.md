# hetfamily: heterogeneous treatment effects for jamovi

Menu: Heterogeneity

| Analysis | What it does |
|---|---|
| **Heterogeneous Effects (BLP, GATES, CLAN)** | Generic machine learning inference of Chernozhukov, Demirer, Duflo and Fernandez-Val: the best linear predictor of the effects (average effect and a test for predictable heterogeneity), the average effects of groups sorted by predicted effect (GATES) with the top-minus-bottom difference, and the covariate means of the least and most affected groups (CLAN) |
| **Priority Curve (RATE)** | Rank-weighted average treatment effects (Yadlowsky et al.): the TOC curve (how much more the units ranked highest by the predicted effect gain than the average unit), its area (AUTOC) and the Qini, with bootstrap standard errors |
| **Effect Profile (CATE regression)** | The doubly robust scores regressed on the moderators: the effect at the average and its change per standard deviation of each moderator (best linear projection of the conditional effect) |
| **Treatment Rule** | Value of treating the units whose predicted effect exceeds a cost: gain over treating nobody and over treating everybody, with the share treated |

Learners: lasso, ridge, elastic net (`glmnet`), random forest (`ranger`), linear regression; several can be compared by the lambda and lambda-bar measures. Each of many random sample splits trains the learners on one part of the sample and estimates on the other; the splits are combined by medians with the p-value and interval adjustments of the paper. The treatment must be randomized (or have a known propensity score).

The engine is base R (`R/engine_het.R`). `tests/reference.R` checks the second stage against the R package `GenericML` on the same predictions: BLP, GATES (with sorting and the top-minus-bottom difference), lambda and lambda-bar for classical, HC1 and HC3 standard errors, and the CLAN means, differences and standard errors, to rounding error (it needs GenericML in a library named by `HET_REFLIB`). `tests/profile.R` checks the effect profile: the regression of doubly robust scores equals `grf::best_linear_projection` (coefficients and HC3 standard errors) and recovers known slopes. `tests/rate.R` checks the priority curve: AUTOC, Qini and the TOC curve equal `grf::rank_average_treatment_effect.fit` on the same scores and priorities (to 1e-10, also between observations), the bootstrap standard error is about the same as grf's, and by simulation the AUTOC is near zero for random priorities, large for perfect ones, and not significant under a constant effect (grf in `HET_REFLIB`). `tests/properties.R` checks the whole procedure by simulation: detection of heterogeneity with a slope near 1, no rejection when the effect is constant, correct ordering and range of the group effects, the value of the rule against the known value. `tests/smoke.R` runs both analysis classes on the sample data and edge cases, also in jamovi's own R.

Example data (`data/`): `sim_het_trial.csv` (2000 people, effect 0.3 + 0.5 female - 0.4 baseline) and `sim_homogeneous_trial.csv` (600 people, constant effect); generator `dev/make_data.R`.

## Ribbon icon

`icon/analysis-het.svg` (74 x 44): five group effects with their intervals rising from the least to the most affected group, the average effect in orange and the null as a dotted line. `dev/icon_preview.ps1` adds one CSS rule for the menu group Heterogeneity to jamovi's client stylesheet (administrator PowerShell, restart jamovi afterwards).
