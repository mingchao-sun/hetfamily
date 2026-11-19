# Rank-weighted average treatment effects (AUTOC, Qini, TOC curve) against grf::rank_average_treatment_effect.fit on the same scores and
# priorities, and properties on simulated experiments. Needs grf in a library named by HET_REFLIB. Run from the module root:  Rscript tests/rate.R
lib <- Sys.getenv("HET_REFLIB"); if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
source("R/engine_het.R")
nfail <- 0
ok <- function(cond, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); if (!isTRUE(cond)) nfail <<- nfail + 1 }

# ---- point estimates and the TOC curve against grf ----
if (requireNamespace("grf", quietly = TRUE)) {
    set.seed(5); n <- 700
    S <- rnorm(n); Gam <- 0.4 + 0.8 * S + rnorm(n, 0, 1.5)
    q <- seq(0.1, 1, by = 0.1)
    mine <- het_rate_point(Gam, S, q)
    for (tg in c("AUTOC", "QINI")) {
        ref <- grf::rank_average_treatment_effect.fit(Gam, S, target = tg, q = q, R = 2)
        ok(abs(mine[if (tg == "AUTOC") "autoc" else "qini"] - ref$estimate) < 1e-10, sprintf("%s point estimate equals grf (%.6f)", tg, ref$estimate))
        if (tg == "AUTOC") ok(max(abs(mine[-(1:2)] - ref$TOC$estimate)) < 1e-10, "TOC curve at 10 grid points equals grf")
    }
    q2 <- c(0.05, 0.13, 0.5, 0.77, 1)                      # grid points that fall between observations
    ref2 <- grf::rank_average_treatment_effect.fit(Gam, S, target = "AUTOC", q = q2, R = 2)
    ok(max(abs(het_rate_point(Gam, S, q2)[-(1:2)] - ref2$TOC$estimate)) < 1e-10, "TOC curve at points between observations equals grf")
    # standard errors: grf uses half-sample bootstrap, here the ordinary bootstrap of the units, so only the size is compared
    sd_grf <- grf::rank_average_treatment_effect.fit(Gam, S, target = "AUTOC", q = q, R = 400)$std.err
    sd_me <- het_rate(Gam, S, q, nboot = 400, seed = 3)["AUTOC", "se"]
    ok(sd_me / sd_grf > 0.7 && sd_me / sd_grf < 1.4, sprintf("AUTOC standard error %.3f against grf %.3f (different bootstrap schemes)", sd_me, sd_grf))
} else cat("grf is not installed: comparison skipped\n")

# ---- properties ----
set.seed(8); n <- 4000
ok(abs(het_rate_point(rnorm(n), rnorm(n), 1)["autoc"]) < 0.06, "random priorities: AUTOC near zero")
Sg <- rnorm(n); tau <- Sg; Gam <- tau + rnorm(n, 0, 1)
r1 <- het_rate(Gam, Sg, nboot = 100)
ok(r1["AUTOC", "est"] > 0.5 && r1["AUTOC", "est"] / r1["AUTOC", "se"] > 5, sprintf("perfect priorities: AUTOC %.2f clearly positive", r1["AUTOC", "est"]))
r0 <- het_rate(Gam, rnorm(n), nboot = 100)
ok(abs(r0["AUTOC", "est"]) < 3 * r0["AUTOC", "se"], "unrelated priorities: AUTOC not significant")
ok(abs(r1["q1", "est"]) < 1e-12 && r1["q0.05", "est"] > r1["q0.5", "est"], "TOC is zero at q = 1 and largest for the top few percent")
had <- exists(".Random.seed", envir = globalenv()); before <- get(".Random.seed", envir = globalenv()); invisible(het_rate(Gam, Sg, nboot = 10))
ok(identical(before, get(".Random.seed", envir = globalenv())), "the bootstrap leaves the random number state of the caller untouched")
# whole procedure
sim <- function(n, seed) { set.seed(seed); x <- matrix(rnorm(n * 5), n, 5); d <- rbinom(n, 1, .5); list(x = x, d = d, y = 1 + x[, 2] + d * (0.5 + x[, 1]) + rnorm(n)) }
s <- sim(2400, 3)
r <- het_generic(s$y, s$d, s$x, NULL, learners = "lasso", K = 3, nsplit = 15, seed = 11, rate = TRUE, nboot = 50)
a <- r$rate$lasso
ok(a["AUTOC", "lo"] > 0 && a["QINI", "lo"] > 0 && a["AUTOC", "p"] < 0.001, sprintf("heterogeneous effect: AUTOC %.2f (interval %.2f to %.2f), Qini %.2f", a["AUTOC", "est"], a["AUTOC", "lo"], a["AUTOC", "hi"], a["QINI", "est"]))
r2 <- het_generic(s$y, s$d, s$x, NULL, learners = "lasso", K = 3, nsplit = 15, seed = 11, rate = FALSE)
ok(isTRUE(all.equal(r$blp$lasso, r2$blp$lasso)), "asking for the priority curve does not change the other results")
rej <- replicate(40, { s0 <- sim(800, sample.int(1e6, 1)); s0$y <- 1 + s0$x[, 2] + s0$d * 0.5 + rnorm(800)
                       het_generic(s0$y, s0$d, s0$x, NULL, learners = "ols", K = 3, nsplit = 11, seed = sample.int(1e5, 1), rate = TRUE, nboot = 30)$rate$ols["AUTOC", "p"] < 0.05 })
ok(mean(rej) <= 0.10, sprintf("constant effect: AUTOC significant in %.0f%% of 40 experiments (nominal 5%%, conservative aggregation)", 100 * mean(rej)))
cat(if (nfail == 0) "\nAll priority curve checks passed.\n" else sprintf("\n%d check(s) FAILED.\n", nfail))
if (nfail > 0) quit(status = 1)
