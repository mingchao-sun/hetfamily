# Effect profile (regression of the doubly robust scores on moderators): against grf::best_linear_projection on the same scores, and recovery of known effects.
# Needs grf in a library named by HET_REFLIB. Run from the module root:  Rscript tests/profile.R
lib <- Sys.getenv("HET_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
source("R/engine_het.R")
nfail <- 0
cmp <- function(label, a, b, tol) { ok <- all(is.finite(a)) && max(abs(a - b)) < tol; cat(sprintf("[%s] %-66s max diff %.1e\n", if (ok) "PASS" else "FAIL", label, max(abs(a - b)))); if (!ok) nfail <<- nfail + 1 }
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }
set.seed(3); n <- 1500
X <- matrix(rnorm(n * 4), n, 4, dimnames = list(NULL, paste0("x", 1:4))); d <- rbinom(n, 1, .5)
tau <- 0.4 + 0.5 * X[, 1] - 0.3 * X[, 2]; y <- 1 + X[, 3] + d * tau + rnorm(n)

if (requireNamespace("grf", quietly = TRUE)) {
    cf <- grf::causal_forest(X, y, d, W.hat = rep(0.5, n), num.trees = 300, seed = 1)
    G <- grf::get_scores(cf)                                  # doubly robust scores of the forest
    bl <- grf::best_linear_projection(cf, X)
    me <- het_wls(G, cbind(1, X), rep(1, n), "hc3")
    cmp("regression of the scores on moderators: coefficients equal grf BLP", me$b, bl[, "Estimate"], 1e-8)
    cmp("... and HC3 standard errors", me$se, bl[, "Std. Error"], 1e-8)
}
# whole procedure: slopes per standard deviation of the moderators; the true slopes are 0.5 sd(x1) and -0.3 sd(x2)
r <- het_generic(y, d, X, NULL, learners = "lasso", K = 3, nsplit = 21, seed = 5, moderators = X)
pf <- r$profile$lasso
ok(all(c("(Intercept)", colnames(X)) %in% rownames(pf)), "the profile has an intercept and one row per moderator")
ok(pf["(Intercept)", "lo"] < 0.4 && pf["(Intercept)", "hi"] > 0.4, sprintf("intercept (effect at the mean) %.2f, interval %.2f to %.2f contains 0.4", pf["(Intercept)", "est"], pf["(Intercept)", "lo"], pf["(Intercept)", "hi"]))
t1 <- 0.5 * sd(X[, 1]); t2 <- -0.3 * sd(X[, 2])
ok(pf["x1", "lo"] < t1 && pf["x1", "hi"] > t1 && pf["x1", "p"] < 0.001, sprintf("slope on x1 %.2f (true %.2f), p %.1e", pf["x1", "est"], t1, pf["x1", "p"]))
ok(pf["x2", "lo"] < t2 && pf["x2", "hi"] > t2, sprintf("slope on x2 %.2f (true %.2f)", pf["x2", "est"], t2))
ok(pf["x3", "p"] > 0.01 && pf["x4", "p"] > 0.01, "irrelevant moderators are not significant")
r2 <- het_generic(y, d, X, NULL, learners = "lasso", K = 3, nsplit = 21, seed = 5, moderators = NULL)
ok(isTRUE(all.equal(r$blp$lasso, r2$blp$lasso)) && length(r2$profile) == 0, "asking for the profile does not change the other results")
# a moderator that is not a learner covariate (a subset of the covariates is enough, any variable works)
r3 <- het_generic(y, d, X, NULL, learners = "ols", K = 3, nsplit = 11, seed = 5, moderators = X[, 1, drop = FALSE])
ok(rownames(r3$profile$ols)[2] == "x1" && nrow(r3$profile$ols) == 2, "a single moderator")
cat(if (nfail == 0) "\nAll effect profile checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
