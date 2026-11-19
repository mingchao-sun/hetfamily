# Second-stage regressions (BLP, GATES, CLAN) and aggregation against the R package GenericML, on the same proxies.
# GenericML needs to be installed in a library named by HET_REFLIB. Run from the module root:  Rscript tests/reference.R
lib <- Sys.getenv("HET_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
source("R/engine_het.R")
if (!requireNamespace("GenericML", quietly = TRUE)) { cat("GenericML is not installed: skipped\n"); quit(status = 0) }
suppressMessages(library(GenericML)); ns <- asNamespace("GenericML")
nfail <- 0
cmp <- function(label, a, b, tol = 1e-8) { ok <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)); cat(sprintf("[%s] %s  (max diff %.2e)\n", if (ok) "PASS" else "FAIL", label, max(abs(as.numeric(a) - as.numeric(b))))); if (!ok) nfail <<- nfail + 1 }

set.seed(11); n <- 1200; p <- 0.5
x <- matrix(rnorm(n * 4), n, 4, dimnames = list(NULL, paste0("x", 1:4))); d <- rbinom(n, 1, p)
tau <- 0.5 + 0.8 * x[, 1]; y <- 1 + x[, 2] + x[, 3] * 0.5 + d * tau + rnorm(n)
B <- 1 + x[, 2] + 0.5 * x[, 3] + rnorm(n, 0, 0.3); S <- tau + rnorm(n, 0, 0.5)      # proxies of some quality
K <- 5

for (vc in c("const", "hc1", "hc3")) {
    ours <- het_main(y, d, rep(p, n), B, S, x, K, vcov = vc, monotonize = TRUE)
    vctl <- GenericML::setup_vcov(estimator = "vcovHC", arguments = list(type = if (vc == "const") "const" else toupper(vc)))
    ref_b <- ns$BLP.classic(D = d, Y = y, propensity_scores = rep(p, n), proxy_BCA = B, proxy_CATE = S, vcov_control = vctl)$generic_targets
    cmp(sprintf("BLP estimates and s.e., vcov %s", vc), c(ours$blp[, "est"], ours$blp[, "se"]), c(ref_b[, "Estimate"], ref_b[, "Std. Error"]))
    grp <- ns$quantile_group_NoChecks(S, seq_len(K - 1) / K)
    ref_g <- ns$GATES.classic(D = d, Y = y, propensity_scores = rep(p, n), proxy_BCA = B, proxy_CATE = S, membership = grp, vcov_control = vctl, monotonize = TRUE)$generic_targets
    cmp(sprintf("GATES estimates and s.e., vcov %s", vc), c(ours$gates[, "est"], ours$gates[, "se"]), c(ref_g[1:K, "Estimate"], ref_g[1:K, "Std. Error"]))
    cmp(sprintf("GATES difference top - bottom, vcov %s", vc), ours$gates_diff[c("est", "se")], c(ref_g[K + 1, "Estimate"], ref_g[K + 1, "Std. Error"]))
}
# lambda parameters
bl <- ns$BLP.classic(D = d, Y = y, propensity_scores = rep(p, n), proxy_BCA = B, proxy_CATE = S)
gt <- ns$GATES.classic(D = d, Y = y, propensity_scores = rep(p, n), proxy_BCA = B, proxy_CATE = S, membership = ns$quantile_group_NoChecks(S, seq_len(K - 1) / K))
lp <- ns$lambda_parameters_NoChecks(BLP = bl, GATES = gt, proxy_CATE = S, membership = ns$quantile_group_NoChecks(S, seq_len(K - 1) / K))
o <- het_main(y, d, rep(p, n), B, S, NULL, K, vcov = "const")
cmp("lambda and lambda.bar", c(o$lambda, o$lambda_bar), c(lp$lambda, lp$lambda.bar))
# CLAN: means by group; the difference of the groups with a Welch standard error
g <- het_groups(S, K); cl <- het_main(y, d, rep(p, n), B, S, x, K)$clan
ref_c <- ns$CLAN_NoChecks(Z_CLAN = x, membership = ns$quantile_group_NoChecks(S, seq_len(K - 1) / K), external_weights = rep(1, n))
cmp("CLAN means of the least and most affected groups", c(cl[, "mean_low"], cl[, "mean_high"]), c(sapply(1:4, function(j) mean(x[g == 1, j])), sapply(1:4, function(j) mean(x[g == K, j]))))
rd <- sapply(1:4, function(j) ref_c$generic_targets[[j]]["delta.5-delta.1", c("Estimate", "Std. Error")])
cmp("CLAN difference between the most and least affected groups, estimate and s.e.", rbind(cl[, "diff"], cl[, "se"]), rd, 1e-6)
cat(if (nfail == 0) "\nAll checks passed.\n" else sprintf("\n%d check(s) failed.\n", nfail))
if (nfail > 0) quit(status = 1)
