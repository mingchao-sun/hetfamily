# Generic machine learning inference on heterogeneous treatment effects in a randomized experiment
# (Chernozhukov, Demirer, Duflo and Fernandez-Val, 2018/2023): best linear predictor (BLP), sorted group average treatment
# effects (GATES), classification analysis (CLAN), and the value of a treatment rule built on the same proxy.
# Each split of the sample trains the learners on an auxiliary half and estimates on the main half; the results of many
# splits are combined by medians, with the p-value and interval adjustments of the paper.
# The second-stage regressions and CLAN equal the R package GenericML on the same proxies (tests/reference.R).

het_learners <- c(ols = "Linear regression", lasso = "Lasso", ridge = "Ridge", enet = "Elastic net (alpha 0.5)", rf = "Random forest")

# A learner is function(Xtr, ytr, Xte) returning predictions of the conditional mean at Xte. X is a numeric matrix without intercept.
het_learner <- function(type, lambda_rule = "lambda.min", nfolds_cv = 5, trees = 300) {
    force(type); force(lambda_rule); force(nfolds_cv); force(trees)
    function(Xtr, ytr, Xte) {
        if (type == "ols") {
            fit <- stats::lm.fit(cbind(1, Xtr), ytr); b <- fit$coefficients; b[is.na(b)] <- 0
            return(as.numeric(cbind(1, Xte) %*% b))
        }
        if (type %in% c("lasso", "ridge", "enet")) {
            alpha <- switch(type, lasso = 1, ridge = 0, enet = 0.5)
            if (ncol(Xtr) < 2) { Xtr <- cbind(Xtr, 0); Xte <- cbind(Xte, 0) }
            cv <- glmnet::cv.glmnet(Xtr, ytr, alpha = alpha, nfolds = nfolds_cv)
            return(as.numeric(stats::predict(cv, Xte, s = lambda_rule)))
        }
        if (type == "rf") {
            if (!requireNamespace("ranger", quietly = TRUE)) stop("The random forest needs the ranger package.", call. = FALSE)
            fit <- ranger::ranger(x = Xtr, y = ytr, num.trees = trees)
            return(as.numeric(stats::predict(fit, Xte)$predictions))
        }
        stop("Unknown learner.", call. = FALSE)
    }
}

# Weighted least squares with a classical or heteroskedasticity-robust covariance matrix, computed as in sandwich::vcovHC
# applied to lm(weights = w): the regression is run on sqrt(w) y and sqrt(w) X.
het_wls <- function(y, X, w, vcov = c("hc3", "hc2", "hc1", "hc0", "const")) {
    vcov <- match.arg(vcov)
    sw <- sqrt(w); Xw <- X * sw; yw <- y * sw
    qrx <- qr(Xw); n <- length(y); k <- qrx$rank
    if (k < ncol(X)) stop("The regressors of the second stage are collinear.", call. = FALSE)
    XtXi <- chol2inv(qr.R(qrx)); b <- as.numeric(XtXi %*% crossprod(Xw, yw)); e <- as.numeric(yw - Xw %*% b)
    V <- if (vcov == "const") XtXi * sum(e^2) / (n - k) else {
        h <- rowSums((Xw %*% XtXi) * Xw)
        om <- switch(vcov, hc0 = e^2, hc1 = e^2 * n / (n - k), hc2 = e^2 / (1 - h), hc3 = e^2 / (1 - h)^2)
        XtXi %*% crossprod(Xw * sqrt(om)) %*% XtXi
    }
    list(b = b, se = sqrt(diag(V)), V = V)
}

# Groups of the units by the quantiles of the proxy (K groups, lowest to highest). Ties are broken by small noise if a group would be too small.
het_groups <- function(S, K) {
    cut <- seq_len(K - 1) / K
    for (tries in 1:4) {
        q <- stats::quantile(S, cut, names = FALSE)
        g <- 1L + rowSums(outer(S, q, ">="))
        if (all(tabulate(g, K) >= 2)) return(g)
        S <- S + stats::rnorm(length(S), 0, 0.001)
    }
    stop("The proxy for the treatment effect has too little variation for this number of groups: use fewer groups.", call. = FALSE)
}

# Rank-weighted average treatment effects (Yadlowsky et al., 2021) for the priority score S and the doubly robust scores Gam of a
# sample: TOC(q) = mean of the scores of the top share q by priority minus the overall mean, AUTOC = average of the TOC curve over all
# ranks, Qini = the same with weights q. Point estimates equal grf::rank_average_treatment_effect.fit (no ties); the standard errors
# come from a bootstrap of the units (the ranking is held fixed), with its own seed so that the random numbers of the learners are untouched.
het_rate_point <- function(Gam, S, q) {
    n <- length(Gam); g <- Gam[order(S, decreasing = TRUE)]
    ate <- mean(g); cm <- cumsum(g) / seq_len(n); toc <- cm - ate
    nw <- q * n; idx <- floor(nw + 1e-9)
    base <- ifelse(idx > 0, cumsum(g)[pmax(idx, 1)], 0)
    nxt <- g[pmin(idx + 1, n)]
    tq <- (base + (nw - idx) * nxt) / nw - ate
    c(autoc = mean(toc), qini = sum(seq_len(n) / n * toc) / n, tq)
}
het_rate <- function(Gam, S, q = seq(0.05, 1, by = 0.05), nboot = 100, seed = 1) {
    est <- het_rate_point(Gam, S, q); n <- length(Gam)
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had) old <- get(".Random.seed", envir = globalenv())
    on.exit(if (had) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    bs <- t(vapply(seq_len(nboot), function(b) { i <- sample.int(n, replace = TRUE); het_rate_point(Gam[i], S[i], q) }, numeric(length(est))))
    se <- apply(bs, 2, stats::sd)
    out <- cbind(est = est, se = se, p = 2 * stats::pnorm(-abs(est / pmax(se, 1e-12))))
    rownames(out) <- c("AUTOC", "QINI", paste0("q", q))
    out
}

# One analysis of the main sample. y, d, p, B (proxy for the control mean), S (proxy for the effect), Z (columns for CLAN, may be NULL),
# all of the main sample. Returns BLP, GATES, CLAN, the parameters lambda and lambda.bar and the value of the rule "treat if S > cost".
het_main <- function(y, d, p, B, S, Z, K, vcov = "hc3", monotonize = TRUE, funsZ = "B", cost = 0, m1 = NULL, m0 = NULL, rate = FALSE, nboot = 100, qgrid = seq(0.05, 1, by = 0.05), rate_seed = 1, prof_Z = NULL) {
    n <- length(y); w <- 1 / (p * (1 - p))
    X1 <- cbind(1, if ("B" %in% funsZ) B, if ("S" %in% funsZ) S, if ("p" %in% funsZ && stats::sd(p) > 0) p)
    dp <- d - p
    # BLP
    fit <- het_wls(y, cbind(X1, beta1 = dp, beta2 = dp * (S - mean(S))), w, vcov)
    k <- ncol(X1); ix <- k + 1:2
    blp <- cbind(est = fit$b[ix], se = fit$se[ix]); blp <- cbind(blp, p = 2 * stats::pnorm(-abs(blp[, 1] / blp[, 2])))
    # GATES
    g <- het_groups(S, K); G <- outer(g, seq_len(K), "==") * 1
    fg <- het_wls(y, cbind(X1, dp * G), w, vcov)
    ig <- k + seq_len(K); ge <- fg$b[ig]; gV <- fg$V[ig, ig, drop = FALSE]; gse <- sqrt(diag(gV))
    if (monotonize) { o <- order(ge); ge <- ge[o]; gV <- gV[o, o, drop = FALSE]; gse <- sqrt(diag(gV)) }          # sorted estimates, as in GenericML
    gdiff <- ge[K] - ge[1]; gdse <- sqrt(gV[K, K] + gV[1, 1] - 2 * gV[K, 1])
    gates <- cbind(est = ge, se = gse, p = 2 * stats::pnorm(-abs(ge / gse)))
    gates_diff <- c(est = gdiff, se = gdse, p = 2 * stats::pnorm(-abs(gdiff / gdse)))
    # lambda parameters for the choice of the learner
    lambda <- as.numeric(blp[2, 1]^2 * stats::var(S)); lambda_bar <- as.numeric(colMeans(G) %*% (fg$b[ig])^2)
    # CLAN: means of the covariates in the least and the most affected group, and their difference
    clan <- NULL
    if (!is.null(Z) && ncol(Z) > 0) {
        clan <- t(vapply(seq_len(ncol(Z)), function(j) {
            lo <- Z[g == 1, j]; hi <- Z[g == K, j]
            sl <- if (length(lo) > 1) stats::var(lo) / length(lo) else 0; sh <- if (length(hi) > 1) stats::var(hi) / length(hi) else 0
            dd <- mean(hi) - mean(lo); se <- sqrt(sl + sh)
            c(mean_low = mean(lo), mean_high = mean(hi), diff = dd, se = se, p = if (se > 0) 2 * stats::pnorm(-abs(dd / se)) else 1)
        }, numeric(5)))
        rownames(clan) <- colnames(Z)
    }
    # value of the rule "treat if the proxy exceeds the cost": doubly robust scores built on the auxiliary-sample predictions
    pol <- NULL
    if (!is.null(m1)) {
        Gam <- (m1 - m0) + d * (y - m1) / p - (1 - d) * (y - m0) / (1 - p)
        tr <- as.numeric(S > cost)
        sc <- function(v) c(est = mean(v), se = stats::sd(v) / sqrt(n))
        pol <- rbind(sc(tr * (Gam - cost)), sc((tr - 1) * (Gam - cost)), sc(Gam), c(est = mean(tr), se = NA_real_))
        rownames(pol) <- c("gain_none", "gain_all", "ate", "share")
    }
    rt <- NULL
    if (rate && !is.null(m1)) {
        Gam <- (m1 - m0) + d * (y - m1) / p - (1 - d) * (y - m0) / (1 - p)
        rt <- het_rate(Gam, S, qgrid, nboot, rate_seed)
    }
    # effect profile: linear projection of the conditional effect on the moderators = regression of the doubly robust scores on them (Semenova and Chernozhukov, 2021)
    prof <- NULL
    if (!is.null(prof_Z) && !is.null(m1)) {
        Gam <- (m1 - m0) + d * (y - m1) / p - (1 - d) * (y - m0) / (1 - p)
        pf <- het_wls(Gam, cbind(1, prof_Z), rep(1, n), vcov)
        prof <- cbind(est = pf$b, se = pf$se, p = 2 * stats::pnorm(-abs(pf$b / pf$se)))
        rownames(prof) <- c("(Intercept)", colnames(prof_Z))
    }
    list(blp = blp, gates = gates, gates_diff = gates_diff, lambda = lambda, lambda_bar = lambda_bar, clan = clan, policy = pol, rate = rt, profile = prof)
}

# Median aggregation across splits (the paper's quantile aggregation): estimates by the median; p-values by min(1, 2 x median);
# the interval bounds are the medians of intervals built with level 1 - alpha/2, so that the interval has level at least 1 - alpha.
het_median <- function(x) { v <- sort(x[is.finite(x)]); if (!length(v)) return(NA_real_); mean(c(stats::quantile(v, .5, type = 1, names = FALSE), -stats::quantile(-v, .5, type = 1, names = FALSE))) }
het_agg <- function(est, se, p = NULL, level = 0.95) {
    z <- stats::qnorm(1 - (1 - level) / 4)                      # level 1 - alpha/2 within each split, alpha = 1 - level
    out <- c(est = het_median(est), se = het_median(se), lo = het_median(est - z * se), hi = het_median(est + z * se))
    c(out, p = if (!is.null(p)) min(1, 2 * het_median(p)) else NA_real_)
}

# The whole procedure. y outcome, d 0/1 treatment, X covariates for the learners (numeric matrix), Z covariates for CLAN (numeric matrix),
# p propensity score (scalar, or vector for the full sample). learners: names from het_learners.
het_generic <- function(y, d, X, Z = NULL, p = NULL, learners = "lasso", K = 5, nsplit = 50, share_main = 0.5, level = 0.95, vcov = "hc3",
                        monotonize = TRUE, funsZ = "B", cost = 0, seed = 1234, lambda_rule = "lambda.min", trees = 300, checkpoint = NULL, rate = FALSE, nboot = 100, qgrid = seq(0.05, 1, by = 0.05), moderators = NULL) {
    n <- length(y)
    PZ <- NULL                                           # moderators standardised with the full-sample mean and standard deviation: the intercept is the effect at the mean
    if (!is.null(moderators)) { PZ <- scale(moderators); PZ[is.na(PZ)] <- 0; attr(PZ, "scaled:center") <- NULL; attr(PZ, "scaled:scale") <- NULL; colnames(PZ) <- colnames(moderators) }
    if (is.null(p)) p <- mean(d)
    if (length(p) == 1) p <- rep(p, n)
    if (any(p <= 0 | p >= 1)) stop("The propensity score must lie strictly between 0 and 1.", call. = FALSE)
    if (K < 2 || K > 10) stop("The number of groups must be between 2 and 10.", call. = FALSE)
    if (min(sum(d == 1), sum(d == 0)) * min(share_main, 1 - share_main) < 4 * K) stop("The sample is too small for this number of groups; use fewer groups.", call. = FALSE)
    lr <- stats::setNames(lapply(learners, het_learner, lambda_rule = lambda_rule, trees = trees), learners)
    res <- stats::setNames(lapply(learners, function(l) vector("list", nsplit)), learners)
    for (s in seq_len(nsplit)) {
        if (!is.null(checkpoint)) checkpoint()
        set.seed(seed + s)
        main <- logical(n)
        for (v in 0:1) { i <- which(d == v); main[sample(i, round(length(i) * share_main))] <- TRUE }
        for (l in learners) {
            tr <- !main
            B <- lr[[l]](X[tr & d == 0, , drop = FALSE], y[tr & d == 0], X[main, , drop = FALSE])
            M1 <- lr[[l]](X[tr & d == 1, , drop = FALSE], y[tr & d == 1], X[main, , drop = FALSE])
            res[[l]][[s]] <- het_main(y[main], d[main], p[main], B, M1 - B, if (!is.null(Z)) Z[main, , drop = FALSE], K, vcov, monotonize, funsZ, cost, m1 = M1, m0 = B, rate = rate, nboot = nboot, qgrid = qgrid, rate_seed = seed + s, prof_Z = if (!is.null(PZ)) PZ[main, , drop = FALSE])
        }
    }
    pick <- function(f) vapply(res, function(r) het_median(vapply(r, f, numeric(1))), numeric(1))
    lam <- cbind(lambda = pick(function(r) r$lambda), lambda_bar = pick(function(r) r$lambda_bar))
    best <- list(blp = learners[which.max(lam[, "lambda"])], gates = learners[which.max(lam[, "lambda_bar"])])
    agg <- function(l, f, rows = NULL) {      # f extracts a matrix (rows x est/se/p) from a split
        arr <- simplify2array(lapply(res[[l]], f))
        if (length(dim(arr)) == 2) arr <- array(arr, c(1, dim(arr)), dimnames = list(NULL, rownames(arr), NULL))
        t(vapply(seq_len(dim(arr)[1]), function(i) het_agg(arr[i, 1, ], arr[i, 2, ], if (dim(arr)[2] >= 3) arr[i, 3, ], level), numeric(5)))
    }
    out <- list(learners = learners, lambda = lam, best = best, nsplit = nsplit, K = K, level = level, n = n, nmain = round(n * share_main), blp = list(), gates = list(), gates_diff = list(), clan = list(), policy = list(), rate = list(), profile = list())
    for (l in learners) {
        out$blp[[l]] <- agg(l, function(r) r$blp); rownames(out$blp[[l]]) <- c("ATE", "HET")
        out$gates[[l]] <- agg(l, function(r) r$gates); rownames(out$gates[[l]]) <- paste0("G", seq_len(K))
        out$gates_diff[[l]] <- agg(l, function(r) matrix(r$gates_diff, 1)); colnames(out$gates_diff[[l]]) <- colnames(out$gates[[l]])
        if (!is.null(Z) && ncol(Z) > 0) {
            nm <- colnames(Z)
            out$clan[[l]] <- do.call(rbind, lapply(seq_along(nm), function(j) {
                one <- function(f) het_median(vapply(res[[l]], function(r) r$clan[j, f], numeric(1)))
                a <- het_agg(vapply(res[[l]], function(r) r$clan[j, "diff"], numeric(1)), vapply(res[[l]], function(r) r$clan[j, "se"], numeric(1)),
                             vapply(res[[l]], function(r) r$clan[j, "p"], numeric(1)), level)
                c(mean_low = one("mean_low"), mean_high = one("mean_high"), diff = a[["est"]], se = a[["se"]], lo = a[["lo"]], hi = a[["hi"]], p = a[["p"]])
            }))
            rownames(out$clan[[l]]) <- nm
        }
        out$policy[[l]] <- do.call(rbind, lapply(c("gain_none", "gain_all", "ate", "share"), function(nmr) {
            if (nmr == "share") return(c(est = het_median(vapply(res[[l]], function(r) r$policy["share", "est"], numeric(1))), se = NA, lo = NA, hi = NA, p = NA))
            het_agg(vapply(res[[l]], function(r) r$policy[nmr, "est"], numeric(1)), vapply(res[[l]], function(r) r$policy[nmr, "se"], numeric(1)),
                    vapply(res[[l]], function(r) 2 * stats::pnorm(-abs(r$policy[nmr, "est"] / r$policy[nmr, "se"])), numeric(1)), level)
        }))
        rownames(out$policy[[l]]) <- c("gain_none", "gain_all", "ate", "share")
    }
    if (!is.null(PZ)) {
        for (l in learners) {
            rn <- rownames(res[[l]][[1]]$profile)
            out$profile[[l]] <- t(vapply(rn, function(r) het_agg(vapply(res[[l]], function(x) x$profile[r, "est"], numeric(1)), vapply(res[[l]], function(x) x$profile[r, "se"], numeric(1)),
                                                                 vapply(res[[l]], function(x) x$profile[r, "p"], numeric(1)), level), numeric(5)))
            rownames(out$profile[[l]]) <- rn
        }
    }
    if (rate) {
        out$qgrid <- qgrid
        for (l in learners) {
            rn <- rownames(res[[l]][[1]]$rate)
            out$rate[[l]] <- t(vapply(rn, function(r) het_agg(vapply(res[[l]], function(x) x$rate[r, "est"], numeric(1)), vapply(res[[l]], function(x) x$rate[r, "se"], numeric(1)),
                                                              vapply(res[[l]], function(x) x$rate[r, "p"], numeric(1)), level), numeric(5)))
            rownames(out$rate[[l]]) <- rn
        }
    }
    out$splits <- res
    out
}
