# Properties of the whole procedure on simulated experiments: power of the heterogeneity test, its size under a constant effect,
# calibration of the BLP slope, ordering of the group effects and the value of the treatment rule. Run from the module root.
source("R/engine_het.R")
nfail <- 0
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }
sim <- function(n, hetero, seed) {
    set.seed(seed); x <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("x", 1:5))); d <- rbinom(n, 1, .5)
    tau <- if (hetero) 0.5 + x[, 1] else rep(0.5, n)
    list(x = x, d = d, y = 1 + x[, 2] + 0.5 * x[, 3] + d * tau + rnorm(n), tau = tau)
}
s <- sim(2400, TRUE, 3)
r <- het_generic(s$y, s$d, s$x, s$x, learners = c("ols", "lasso"), K = 5, nsplit = 21, seed = 11, cost = 0)
b <- r$blp$lasso
ok(b["HET", "p"] < 0.001 && b["HET", "lo"] > 0.4 && b["HET", "hi"] < 1.6, sprintf("heterogeneity detected: slope %.2f (interval %.2f to %.2f), p %.1e", b["HET", "est"], b["HET", "lo"], b["HET", "hi"], b["HET", "p"]))
ok(b["ATE", "lo"] < 0.5 && b["ATE", "hi"] > 0.5, sprintf("ATE interval (%.2f, %.2f) contains the true average effect 0.5", b["ATE", "lo"], b["ATE", "hi"]))
g <- r$gates$lasso
ok(all(diff(g[, "est"]) >= 0), "group effects are sorted (monotonised)")
ok(g["G5", "est"] > 1.2 && g["G1", "est"] < -0.2, sprintf("group effects span the true range: %.2f to %.2f (true about -0.7 to 1.7)", g["G1", "est"], g["G5", "est"]))
gd <- r$gates_diff$lasso; ok(gd[, "p"] < 0.001 && gd[, "est"] > 1.5, sprintf("top minus bottom group effect %.2f, p %.1e", gd[, "est"], gd[, "p"]))
cl <- r$clan$lasso; ok(cl["x1", "p"] < 0.001 && cl["x1", "diff"] > 1 && all(cl[-1, "p"] > 0.001), "CLAN: only x1 differs between the most and least affected groups")
ok(identical(rownames(r$lambda), c("ols", "lasso")) && all(r$lambda > 0), "learner comparison: lambda and lambda.bar are positive for both learners")
po <- r$policy$lasso
ok(po["gain_none", "est"] > 0.25 && po["gain_none", "lo"] > 0, sprintf("treating the units with a positive predicted effect gains %.2f over treating nobody (interval from %.2f)", po["gain_none", "est"], po["gain_none", "lo"]))
ok(po["gain_all", "est"] > 0, sprintf("... and %.2f over treating everybody", po["gain_all", "est"]))
ok(abs(po["share", "est"] - mean(s$tau > 0)) < 0.12, sprintf("share treated by the rule %.2f (true share with a positive effect %.2f)", po["share", "est"], mean(s$tau > 0)))
# Size of the heterogeneity test under a constant effect, and behaviour of the aggregated p-value
rej <- replicate(60, { s0 <- sim(800, FALSE, sample.int(1e6, 1)); a <- het_generic(s0$y, s0$d, s0$x, NULL, learners = "ols", K = 4, nsplit = 11, seed = sample.int(1e5, 1))
                       c(a$blp$ols["HET", "p"] < 0.05, a$gates_diff$ols[, "p"] < 0.05) })
ok(mean(rej[1, ]) <= 0.10, sprintf("constant effect: the BLP heterogeneity test rejects in %.0f%% of 60 experiments (nominal 5%%, conservative aggregation)", 100 * mean(rej[1, ])))
ok(mean(rej[2, ]) <= 0.12, sprintf("constant effect: the top-bottom GATES test rejects in %.0f%%", 100 * mean(rej[2, ])))
# Reproducibility, and errors
r2 <- het_generic(s$y, s$d, s$x, s$x, learners = "lasso", K = 5, nsplit = 21, seed = 11)
ok(all.equal(r2$blp$lasso, r$blp$lasso) == TRUE, "the same seed gives the same result")
ok(inherits(try(het_generic(s$y[1:40], s$d[1:40], s$x[1:40, ], NULL, K = 10, nsplit = 3), silent = TRUE), "try-error"), "a sample too small for the number of groups is refused")
ok(inherits(try(het_generic(s$y, s$d, s$x, NULL, p = rep(1, length(s$y)), nsplit = 3), silent = TRUE), "try-error"), "a propensity score of 1 is refused")
cat(if (nfail == 0) "\nAll property checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
