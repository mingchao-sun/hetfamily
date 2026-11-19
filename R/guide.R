# "Data and setup guide" and "Method notes" for each analysis.

het_settings <- function(rows) {
    paste0("<table style='border-collapse:collapse;margin:4px 0 8px 0'>",
           "<tr><th style='text-align:left;padding:2px 14px 2px 0;border-bottom:1px solid #999'>Field</th>",
           "<th style='text-align:left;padding:2px 0;border-bottom:1px solid #999'>Set to</th></tr>",
           paste0("<tr><td style='padding:2px 14px 2px 0'>", names(rows), "</td><td style='padding:2px 0'>",
                  unlist(rows), "</td></tr>", collapse = ""), "</table>")
}

het_guide_raw <- function(key) {
    switch(key,
    hetgeneric = paste0(
        "<h3>Heterogeneous Effects: data and setup</h3>",
        "<p>Asks whether the effect of a randomized treatment differs between people and who gains most. It needs one row per unit, a treatment variable with a control and a treated group, an outcome, and the baseline covariates that may change the effect (measured before the treatment).</p>",
        "<h4>Sample data: heterogeneity found</h4>",
        het_settings(list("Dataset" = "<b>Simulated heterogeneous trial</b>", "Outcome" = "y", "Treatment (arms)" = "arm (control group: control, treatment group: treated; leave both blank)",
                          "Covariates" = "age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2", "Learners" = "Lasso")),
        "<p>2000 people randomized with probability 0.5. The true effect is 0.3 + 0.5 female - 0.4 baseline: average 0.55, larger for women and for people with a low baseline. Expected: average effect 0.533 (s.e. 0.064, interval 0.392 to 0.677); slope on the predicted effect 0.972 (s.e. 0.140, interval 0.667 to 1.278, p &lt; .001); group effects rising from -0.120 (group 1) to 1.153 (group 5), the difference 1.271 (s.e. 0.200, interval 0.823 to 1.718); in the classification table female is 0.20 in the least and 0.85 in the most affected group (difference 0.655, p &lt; .001) and baseline 1.212 against -1.160 (difference -2.366, p &lt; .001), while the noise variables do not differ. With the random forest ticked as well, lasso is still the learner with the larger lambda (0.200 against 0.129).</p>",
        "<h4>Sample data: no heterogeneity</h4>",
        het_settings(list("Dataset" = "<b>Simulated homogeneous trial</b>", "Outcome" = "y", "Treatment (arms)" = "arm", "Covariates" = "x1, x2, x3", "Number of groups" = "3")),
        "<p>The effect is 0.4 for everybody. Expected: average effect 0.479 (interval 0.219 to 0.736), slope on the predicted effect -0.336 (interval -2.162 to 1.610, p 0.979) and top minus bottom group difference 0.316 (p 0.523): no heterogeneity is found. In the classification table x1 appears to differ between the groups, but the groups are formed from noise in the predictions: read that table only when heterogeneity has been found.</p>",
        "<h4>Reading the tables</h4><p>The <b>slope</b> in the first table tests whether the learner finds real heterogeneity; the group effects show how much the effect rises from the lowest to the highest predicted group; the last table shows who the most affected are. Heavy computations run many random sample splits: reduce the number of splits to try things out.</p>"),
    hetrate = paste0(
        "<h3>Priority Curve: data and setup</h3>",
        "<p>Asks whether treating first the units with the largest predicted effect really puts the largest gains first. The data are the same as for <i>Heterogeneous Effects</i>.</p>",
        "<h4>Sample data: a useful ranking</h4>",
        het_settings(list("Dataset" = "<b>Simulated heterogeneous trial</b>", "Outcome" = "y", "Treatment (arms)" = "arm", "Covariates" = "age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2", "Learners" = "Lasso")),
        "<p>Expected: AUTOC 0.400 (s.e. 0.065, interval 0.252 to 0.547, p &lt; .001) and Qini 0.124 (s.e. 0.018, interval 0.084 to 0.165). The 10% of the units with the largest predicted effect gain 0.829 more than the average unit (interval 0.392 to 1.270), the top half 0.340 more, and the curve falls to zero when everybody is treated.</p>",
        "<h4>Sample data: no heterogeneity</h4>",
        het_settings(list("Dataset" = "<b>Simulated homogeneous trial</b>", "Outcome" = "y", "Treatment (arms)" = "arm", "Covariates" = "x1, x2, x3")),
        "<p>Expected: AUTOC -0.052 (interval -0.323 to 0.202, p 1.000) and Qini -0.020 (p 0.890): the ranking is no better than a random order.</p>"),
    hetprofile = paste0(
        "<h3>Effect Profile: data and setup</h3>",
        "<p>Describes how the treatment effect changes with each covariate: the doubly robust score of every unit (an unbiased, noisy measure of its own effect) is regressed on the covariates, which gives the best linear approximation of the conditional effect. It answers questions such as <i>by how much does the effect rise when baseline falls by one standard deviation</i>. The data are the same as for <i>Heterogeneous Effects</i>.</p>",
        "<h4>Sample data</h4>",
        het_settings(list("Dataset" = "<b>Simulated heterogeneous trial</b>", "Outcome" = "y", "Treatment (arms)" = "arm", "Covariates" = "age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2", "Moderators" = "leave empty (all covariates)", "Learners" = "Lasso")),
        "<p>The true effect is 0.3 + 0.5 female - 0.4 baseline; per standard deviation of the covariate that is +0.25 for female and -0.40 for baseline. Expected: effect at the average of the moderators 0.538 (s.e. 0.063, interval 0.395 to 0.678); female +0.234 (s.e. 0.063, p &lt; .001), baseline -0.410 (s.e. 0.064, p &lt; .001); the other eight moderators between -0.045 and +0.047, none significant. ",
        "With only <i>female</i> and <i>baseline</i> as moderators (and covariates age, female, baseline, noise1): 0.544, +0.239 and -0.436.</p>",
        "<h4>Data layout</h4><p>One row per unit. Moderators must be among the covariates (they are what the learners and the regression use); leave the box empty to use all covariates. Nominal moderators enter as 0/1 indicators and are standardised like the others.</p>"),
    hetpolicy = paste0(
        "<h3>Treatment Rule Value: data and setup</h3>",
        "<p>Estimates how much better off the population would be if only the units whose predicted effect exceeds a cost were treated. The data are the same as for <i>Heterogeneous Effects</i>.</p>",
        "<h4>Sample data</h4>",
        het_settings(list("Dataset" = "<b>Simulated heterogeneous trial</b>", "Outcome" = "y", "Treatment (arms)" = "arm", "Covariates" = "age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2", "Cost of treating a unit" = "0")),
        "<p>Expected: the rule treats 87.6% of the units; the gain over treating nobody is 0.554 (s.e. 0.060, interval 0.418 to 0.691), the gain over treating everybody is 0.019 (interval -0.030 to 0.067, p 0.758): nearly everybody gains, so selecting hardly helps. The average effect is 0.531. Set the cost to 0.5: the rule then treats 54.5% of the units and the gain over treating everybody becomes 0.153 (interval 0.055 to 0.253).</p>"),
    stop("unknown guide: ", key))
}

het_notes_raw <- function(key) {
    switch(key,
    hetgeneric = paste0(
        "<h3>Heterogeneous Effects: method notes</h3>",
        "<h4>The idea</h4><p>Machine learning can predict who gains most, but its predictions are noisy and, used naively on the same data, give misleading inference. The generic approach of Chernozhukov, Demirer, Duflo and Fernandez-Val (2018) avoids this: split the sample at random into two halves; train the learner on one half, which gives a prediction of the outcome without treatment (B) and of the effect (S), then estimate on the other half. Because the treatment is randomized, the estimates are valid whatever the quality of the learner.</p>",
        "<h4>What is estimated</h4><ul><li><b>Best linear predictor (BLP)</b>: a regression of the outcome on B and on the treatment and the treatment times the demeaned prediction S, weighted by 1/(p(1-p)). The first coefficient is the average effect; the second is the slope on the prediction. A slope of 1 means the predictions are right on average, 0 that they carry no information; a significant positive slope shows heterogeneity.</li>",
        "<li><b>Sorted group average treatment effects (GATES)</b>: the average effect in groups formed by quantiles of S. A large and significant difference between the top and bottom groups shows who gains.</li>",
        "<li><b>Classification analysis (CLAN)</b>: the average covariates of the least and most affected groups, to describe who they are.</li></ul>",
        "<h4>Splits and learners</h4><p>One split gives a result that depends on the luck of the split. The analysis is repeated for many random splits and combined: estimates by the median over splits, the p-value by twice the median p-value, and the confidence interval by the medians of the bounds of intervals with level 1 - alpha/2. This costs a little power but makes the inference valid for the whole procedure. With several learners, the one with the largest lambda (BLP) or lambda-bar (GATES, CLAN) is shown.</p>",
        "<h4>Limits</h4><p>Randomized experiments (or a known propensity score), a binary comparison and a continuous outcome. Covariates are used as given: add the transformations and interactions you think matter. Absence of significant heterogeneity does not mean a constant effect; it means the learner could not find a predictable difference at this sample size. The group effects are conditional on the predictions, which depend on the split. Second-stage results equal the R package GenericML on the same predictions; the learners and random numbers differ between packages, so complete results are close but not identical.</p>"),
    hetrate = paste0(
        "<h3>Priority Curve: method notes</h3>",
        "<h4>What it does</h4><p>Rank the units by the predicted effect, from the largest down. The <b>TOC curve</b> (treatment operating characteristic) shows, for each share q of the units treated first, how much more the top q gain on average than the average unit. If the ranking is informative the curve starts high and falls to zero at q = 1; if it is no better than random it stays near zero. The <b>AUTOC</b> is the area under the curve (the average of the TOC over all ranks) and the <b>Qini</b> weights the curve by q; both are zero for a random ranking and positive for a useful one (Yadlowsky et al., 2021).</p>",
        "<h4>Estimation</h4><p>As in the other analyses of this module the learner is trained on one part of the sample and the curve is evaluated on the other part, with doubly robust scores (the difference of the predicted outcomes plus a correction using the realized outcome and the known assignment probability), so the estimates are valid whatever the quality of the learner. The point estimates equal grf::rank_average_treatment_effect. Standard errors come from a bootstrap of the evaluation units; results of many random splits are combined by medians with the adjustments of Chernozhukov et al. (2018).</p>",
        "<h4>Limits</h4><p>AUTOC is a global measure: it can be positive with a ranking that is useful only for the top few percent, and the Qini gives more weight to those. The curve says nothing about costs; use <i>Treatment Rule Value</i> for a rule with a cost. Randomized experiments (or a known propensity score), a binary comparison and a continuous outcome.</p>"),
    hetprofile = paste0(
        "<h3>Effect Profile: method notes</h3>",
        "<h4>What it does</h4><p>As in the other analyses of this module the learners are trained on one part of the sample and the estimation is done on the other. For each unit of the estimation part the <b>doubly robust score</b> (the difference of the predicted outcomes plus a correction using the realized outcome and the assignment probability) has the conditional effect as its expectation given the covariates, so regressing the scores on the moderators gives the best linear approximation of the conditional effect (Semenova and Chernozhukov, 2021), valid whatever the quality of the learners. The standard errors are robust (HC3).</p>",
        "<h4>Reading the table</h4><ul><li>The moderators are standardised with the sample mean and standard deviation: the first row is the effect for a unit at the average of all moderators and every other row is the change of the effect per standard deviation of that moderator, holding the others fixed in the linear approximation.</li>",
        "<li>A significant row says the effect varies with that variable; it does not say the variable causes the variation (it may stand for another one). Compare with <i>Heterogeneous Effects</i>: CLAN describes who is most affected without regression, this table gives marginal changes.</li>",
        "<li>With many moderators some will look significant by chance; results of many random splits are combined with the conservative adjustments of the module (p-values twice the median).</li></ul>",
        "<h4>Limits</h4><p>A linear approximation: it cannot show non-linear or interaction patterns (add squares or products as variables to explore them). Randomized experiments (or a known propensity score), a binary comparison and a continuous outcome. The regression of scores on the moderators equals grf::best_linear_projection on the same scores.</p>"),
    hetpolicy = paste0(
        "<h3>Treatment Rule Value: method notes</h3>",
        "<h4>What it does</h4><p>A treatment rule treats a unit when its predicted effect exceeds the cost of treating it. The value of the rule is its average effect on the outcome net of cost. The analysis learns the predicted effects on one half of the sample and evaluates the rule on the other half with doubly robust scores (the difference of the predicted outcomes plus a correction using the realized outcome and the known assignment probability), so the estimate is valid whatever the quality of the learner. It is repeated over random splits and combined by medians as in <i>Heterogeneous Effects</i>.</p>",
        "<h4>Reading the table</h4><p>The gain over treating nobody says whether the rule helps at all; the gain over treating everybody says whether being selective pays. If the treatment helps almost everybody, the second gain is near zero. The cost is in outcome units: with a positive cost, only units whose predicted effect exceeds it are treated.</p>",
        "<h4>Limits</h4><p>The rule is not the optimal rule but one based on the learner's predictions of the effect; the value is estimated for the rule learned in each split. Treatment capacity limits and fairness constraints are not modelled. See also Kitagawa and Tetenov (2018) on empirical welfare maximization.</p>"),
    stop("unknown notes: ", key))
}

het_guide_html <- function(key) {
    x <- het_guide_raw(key)
    sub("</p>", paste0("</p><p><i>Open the sample data from File &gt; Open &gt; Data Library &gt; Heterogeneous Treatment Effects Family; the library lists each data set by the bold name below.</i></p>"), x, fixed = TRUE)
}
het_notes_html <- function(key) het_notes_raw(key)
