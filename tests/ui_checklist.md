# hetfamily: checklist for testing in jamovi

Menu: Heterogeneity. Install: `jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`.
Open sample data from File > Open > Data Library > Heterogeneous Treatment Effects Family. Defaults: lasso, 5 groups, 50 splits, seed 1234; with these the numbers below must appear exactly. The analysis takes a few seconds (a progress bar shows).

## 1. Heterogeneous Effects (BLP, GATES, CLAN)
Dataset: Simulated heterogeneous trial. Outcome `y`; Treatment `arm`; Covariates: all of age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2 (drag all ten).
Expect: 2000 units (1012 treated, 988 control); average effect 0.533 (s.e. 0.064, interval 0.392 to 0.677); heterogeneity slope 0.972 (s.e. 0.140, interval 0.667 to 1.278, p < .001); group effects -0.120, 0.303, 0.552, 0.770, 1.153 and the difference 1.271 (s.e. 0.200); CLAN: female 0.20 against 0.85, baseline 1.212 against -1.160, p < .001, noise1 and noise2 not significant. The plot shows five points rising, the orange band is the average effect.
Then: tick Random forest as well -> a comparison table with lambda 0.200 (lasso) and 0.129 (forest), "yes" in the lasso row; a few more seconds. Number of groups 3 or 10; splits 20 (faster, slightly different numbers); share used for estimation 0.3; standard errors Classical / HC1; untick sorting; untick Classification analysis (the table disappears) and the plot. Known propensity score: add a variable that is constant 0.5 (the numbers hardly change). Method notes and guide on/off.
Dataset: Simulated homogeneous trial. Outcome `y`; Treatment `arm`; Covariates x1, x2, x3; 3 groups. Expect: slope -0.336 (p 0.979), top minus bottom 0.316 (p 0.523); the CLAN note warns that the groups may reflect noise.
Messages: untick all learners -> "Select at least one learner"; put `y` also in the covariates -> "A variable can sit in one box only"; treatment group value `zzz` -> message listing the groups; a known propensity score equal to 1 -> message; no covariates -> hint.

## 2. Treatment Rule
Dataset: Simulated heterogeneous trial, same variables. Cost 0.
Expect: the rule treats 87.6% of the units; gain over treating nobody 0.554 (interval 0.418 to 0.691); gain over treating everybody 0.019 (interval -0.030 to 0.067, p 0.758); average effect 0.531.
Then: cost 0.5 -> share treated 54.5%, gain over treating everybody 0.153 (interval 0.055 to 0.253); a negative cost or a very large cost (no unit treated, share 0) still gives a table. Two learners -> comparison table.

## 3. Priority Curve (RATE)
Dataset: Simulated heterogeneous trial, same variables as above (ten covariates, lasso).
Expect: AUTOC 0.400 (s.e. 0.065, interval 0.252 to 0.547, p < .001), Qini 0.124 (s.e. 0.018); the TOC table has ten rows: 0.829 at share 0.1 (interval 0.392 to 1.270), 0.340 at 0.5, 0.087 at 0.9 and 0 at 1 (no standard error in that row). The plot shows a curve falling from the left to zero at the right with a shaded band; the caption gives the AUTOC.
Then: tick Random forest -> comparison table with "yes" in the lasso row; bootstrap samples 20 (faster); untick the plot (it disappears); level 90; method notes and guide on/off.
Dataset: Simulated homogeneous trial, covariates x1, x2, x3. Expect: AUTOC -0.052 (p 1.000), Qini -0.020 (p 0.890); the curve stays around zero.
Messages: untick all learners -> "Select at least one learner"; no covariates -> hint.

## 4. Effect Profile (CATE regression)
Dataset: Simulated heterogeneous trial; Outcome `y`; Treatment `arm`; Covariates: all ten (age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2); Moderators: leave empty; lasso.
Expect: effect at the average of the moderators 0.538 (s.e. 0.063, interval 0.395 to 0.678); change per standard deviation: female +0.234 (s.e. 0.063, p < .001), baseline -0.410 (s.e. 0.064, p < .001), the other eight between -0.045 and +0.047 and not significant. The coefficient plot shows ten bars with intervals; female to the right and baseline to the left of the dashed zero line.
Then: put `female` and `baseline` in the Moderators box (covariates unchanged): the table has three rows (0.544, +0.239, -0.436 when the covariates are only age, female, baseline, noise1); tick Random forest (comparison table); standard errors HC1 / Classical; untick the plot. Method notes and guide on/off.
Messages: a moderator that is not among the covariates -> message; no covariates -> hint; untick all learners -> "Select at least one learner".
