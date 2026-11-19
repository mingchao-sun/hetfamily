# hetfamily — heterogeneous treatment effects in jamovi

A jamovi module for finding out who gains from a treatment in a randomized experiment: generic machine learning inference on heterogeneous effects
(Chernozhukov, Demirer, Duflo and Fernandez-Val), priority curves, effect profiles and treatment rules. Chinese summary: [README.zh.md](README.zh.md).

## Analyses (jamovi menu: Heterogeneity)

| Analysis | What it does |
|---|---|
| **Heterogeneous Effects (BLP, GATES, CLAN)** | The best linear predictor of the effects (average effect and a test for predictable heterogeneity), the average effects of groups sorted by predicted effect (GATES) with the top-minus-bottom difference, and the covariate means of the least and most affected groups (CLAN) |
| **Priority Curve (RATE)** | Rank-weighted average treatment effects (Yadlowsky et al.): the TOC curve (how much more the units ranked highest by the predicted effect gain than the average unit), its area (AUTOC) and the Qini, with bootstrap standard errors |
| **Effect Profile (CATE regression)** | The doubly robust scores regressed on the moderators: the effect at the average and its change per standard deviation of each moderator (best linear projection of the conditional effect) |
| **Treatment Rule** | Value of treating the units whose predicted effect exceeds a cost: gain over treating nobody and over treating everybody, with the share treated |

Learners: lasso, ridge, elastic net (`glmnet`), random forest (`ranger`), linear regression; several can be compared by the lambda and lambda-bar measures.
Each of many random sample splits trains the learners on one part of the sample and estimates on the other; the splits are combined by medians with the
p-value and interval adjustments of the paper. The treatment must be randomized (or have a known propensity score). Every analysis has *Method notes* and
*Data and setup guide* in its Output section.

## Install
- **jamovi library**: in jamovi open the modules button (top right), *jamovi library*, find *Heterogeneous Treatment Effects Family* and install (once the module is listed).
- **Sideload**: download the `.jmo` file and use *Install from file* in the same menu. A `.jmo` works only on the operating system, processor type and jamovi version it was built for.
- **From source** (needs R and the jamovi compiler `jmvtools`; adjust the path to your jamovi folder):
```r
jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")
```

## Use from R (Rj editor, scripts)
Every analysis is also an R function, for example
```r
hetfamily::hetgeneric(data = d, outcome = "y", treat = "arm", covs = c("age", "female", "income", "baseline"), nsplit = 100)
```

## Data format
One row per person. A binary treatment (randomized; the first level is the control), a numeric outcome and numeric or categorical covariates that may
modify the effect. Rows with missing values are dropped and reported.

## Example data
Shown in jamovi under **Open > Data Library > Heterogeneous Treatment Effects Family**; each entry says which analysis to use and which variable goes where.
Files in `data/` (generator `dev/make_data.R`):

| File | Design |
|---|---|
| `sim_het_trial.csv` | 2,000 people randomized with probability 0.5; effect 0.3 + 0.5 female - 0.4 baseline |
| `sim_homogeneous_trial.csv` | 600 people; constant effect (no heterogeneity to find) |

## References
Papers and R packages behind each analysis appear in the References section of the jamovi output (defined in `jamovi/00refs.yaml`).

## Validation
The engine is base R (`R/engine_het.R`).
- `tests/reference.R` checks the second stage against the R package `GenericML` on the same predictions: BLP, GATES (with sorting and the top-minus-bottom difference), lambda and lambda-bar for classical, HC1 and HC3 standard errors, and the CLAN means, differences and standard errors, to rounding error (needs `GenericML` in a library named by `HET_REFLIB`).
- `tests/profile.R`: the regression of doubly robust scores equals `grf::best_linear_projection` (coefficients and HC3 standard errors) and recovers known slopes.
- `tests/rate.R`: AUTOC, Qini and the TOC curve equal `grf::rank_average_treatment_effect.fit` on the same scores and priorities (to 1e-10); by simulation the AUTOC is near zero for random priorities, large for perfect ones, and not significant under a constant effect.
- `tests/properties.R` checks the whole procedure by simulation: detection of heterogeneity with a slope near 1, no rejection when the effect is constant, correct ordering and range of the group effects, and the value of the rule against the known value.
- `tests/smoke.R` runs every analysis headless, including bad input.

## Known limitations
- Randomized treatment (or a known propensity score) only: with confounding, use `matchfamily` or `dmlfamily` first.
- One binary treatment; no multi-arm or continuous treatments.
- Results depend on the random splits (seed option, number of splits); more splits make them stable at the cost of time.
- Heterogeneity that exists but is small relative to the sample size will not be detected: the tests have limited power in samples of a few hundred.
- The effect profile is a linear summary of the conditional effect, not a non-parametric estimate.

## Notes for developers
- Layout: `R/engine_het.R` is plain R with no jamovi dependency; `R/*.b.R` and `jamovi/*.yaml` are the jamovi front end (`dev/gen_het.py` generates the analysis definitions).
- Tests: `Rscript tests/smoke.R`, `tests/reference.R`, `tests/profile.R`, `tests/rate.R`, `tests/properties.R`; `tests/ui_checklist.md` lists what to check by hand in jamovi.
- Ribbon icon: `icon/analysis-het.svg` (74 x 44): five group effects rising from the least to the most affected group. jmvtools 28.17 cannot attach an icon to a module, so `dev/icon_preview.ps1` adds one CSS rule for the menu group to jamovi's client stylesheet (administrator PowerShell, restart jamovi afterwards).

## How to cite
Sun, M. *hetfamily: Heterogeneous Treatment Effects Family for jamovi* (jamovi module). https://github.com/mingchao-sun/hetfamily. Please also cite Chernozhukov, Demirer, Duflo and Fernandez-Val and the papers behind the analysis you use.

## License
GPL-3, as stated in `DESCRIPTION`; third-party code and data are listed in `NOTICE.md`. Bug reports and suggestions: https://github.com/mingchao-sun/hetfamily/issues
