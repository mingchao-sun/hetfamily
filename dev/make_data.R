# Sample data of hetfamily. Run from the module root:  Rscript dev/make_data.R
set.seed(20261117)
n <- 2000
age <- round(rnorm(n, 40, 11)); female <- rbinom(n, 1, 0.5); income <- round(exp(rnorm(n, 3.2, 0.5)), 1)
baseline <- round(rnorm(n), 3); edu <- sample(1:4, n, TRUE, c(.2, .4, .3, .1)); urban <- rbinom(n, 1, 0.6)
hhsize <- rpois(n, 2.5) + 1; health <- round(rnorm(n), 3)
noise1 <- round(rnorm(n), 3); noise2 <- round(rnorm(n), 3)
arm <- ifelse(rbinom(n, 1, 0.5) == 1, "treated", "control")
# effect: larger for women and for people with a low baseline score; zero for most noise variables
tau <- 0.3 + 0.5 * female - 0.4 * baseline
y0 <- 1 + 0.6 * baseline + 0.02 * (age - 40) + 0.3 * health + 0.15 * urban + rnorm(n)
y <- round(y0 + (arm == "treated") * tau, 3)
d <- data.frame(id = seq_len(n), arm = arm, age, female, income, baseline, edu, urban, hhsize, health, noise1, noise2, y)
write.csv(d, "data/sim_het_trial.csv", row.names = FALSE)
cat("true average effect", round(mean(tau), 3), " sd of the effect", round(sd(tau), 3), " share with positive effect", round(mean(tau > 0), 3), "\n")
# a small noisy example without heterogeneity
set.seed(20261118); m <- 600
x1 <- round(rnorm(m), 3); x2 <- round(rnorm(m), 3); x3 <- rbinom(m, 1, .5); a2 <- ifelse(rbinom(m, 1, .5) == 1, "treated", "control")
d2 <- data.frame(id = seq_len(m), arm = a2, x1, x2, x3, y = round(0.5 + 0.5 * x1 + (a2 == "treated") * 0.4 + rnorm(m), 3))
write.csv(d2, "data/sim_homogeneous_trial.csv", row.names = FALSE)
