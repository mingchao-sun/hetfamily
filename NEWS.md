# hetfamily 1.0.0

* First public release (jamovi library submission). Contains everything listed below; the earlier 0.x entries are development history.

# hetfamily 0.2.0

* New: **Effect Profile (CATE regression)**: the doubly robust scores of the estimation sample regressed on the moderators (standardised), giving the effect at the average and its change per standard deviation of each moderator, with the same sample-splitting and median aggregation as the other analyses (Semenova and Chernozhukov, 2021). The regression equals `grf::best_linear_projection` on the same scores.
* New (0.1.x, listed here): **Priority Curve (RATE)** with AUTOC and Qini, equal to `grf`.

# hetfamily 0.1.0

First version. Three analyses in the menu group Heterogeneity:

* **Heterogeneous Effects (BLP, GATES, CLAN)**: generic machine learning inference of Chernozhukov, Demirer, Duflo and Fernandez-Val: the best linear predictor of the effects (average effect and a heterogeneity test), the average effects of groups sorted by predicted effect (GATES, with the difference between the most and least affected group), and the covariate means of those groups (CLAN). Learners: lasso, ridge, elastic net, random forest, linear regression; many random sample splits combined by medians; learners compared by the lambda measures.
* **Priority Curve (RATE)**: rank-weighted average treatment effects of Yadlowsky et al. (2021): the TOC curve for the share of units treated first, the AUTOC and the Qini, to judge whether ranking by predicted effect puts the largest gains first.
* **Treatment Rule Value**: the gain from treating the units whose predicted effect exceeds a cost, compared with treating nobody and with treating everybody.

The second-stage regressions, group effects and classification analysis equal the R package `GenericML` on the same predictions, and the priority curve equals `grf`. Method notes and a data and setup guide are in the Output section of each analysis.
