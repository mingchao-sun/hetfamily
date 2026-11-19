# Helpers shared by the heterogeneity analyses.

het_fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, drop0trailing = TRUE)
het_guard <- function(expr) tryCatch(expr, error = function(e) jmvcore::reject(conditionMessage(e), code = ""))
`%||%` <- function(a, b) if (is.null(a)) b else a

het_show_text <- function(results, opts, key) {
    if (isTRUE(opts$showGuide)) results$guide$setContent(het_guide_html(key))
    if (isTRUE(opts$showHelp)) results$helpText$setContent(het_notes_html(key))
}

# Treatment variable as a factor whose first level is the control group.
het_arm <- function(x, control = "") {
    f <- if (is.factor(x)) droplevels(x) else factor(x, levels = sort(unique(x[!is.na(x)])))
    lv <- levels(f)
    if (length(lv) < 2) stop("The treatment needs at least two groups.", call. = FALSE)
    if (length(lv) > 12) stop("The treatment has more than 12 groups.", call. = FALSE)
    control <- trimws(control %||% "")
    if (nzchar(control)) {
        if (!control %in% lv) stop(sprintf("The control group value '%s' is not among the treatment values (%s).", control, paste(lv, collapse = ", ")), call. = FALSE)
        f <- stats::relevel(f, ref = control)
    }
    f
}

het_num <- function(v, label) {
    if (is.factor(v)) stop(sprintf("The %s must be a numeric variable (set its type to Continuous in the Data tab).", label), call. = FALSE)
    out <- as.numeric(v)
    if (any(is.infinite(out), na.rm = TRUE)) stop(sprintf("The %s has infinite values.", label), call. = FALSE)
    out
}

# Covariates as a numeric matrix (nominal variables become dummies, constant columns are dropped); missing values stay NA.
het_covmat <- function(data, covs) {
    if (length(covs) == 0) return(NULL)
    df <- data[, covs, drop = FALSE]
    for (v in covs) if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]]) else df[[v]] <- het_num(df[[v]], paste("covariate", v))
    f <- stats::reformulate(paste0("`", covs, "`"))
    mm <- stats::model.matrix(f, stats::model.frame(f, df, na.action = stats::na.pass))[, -1, drop = FALSE]
    mm[, apply(mm, 2, function(z) stats::sd(z, na.rm = TRUE)) > 0, drop = FALSE]
}

# Everything the analyses need: outcome, 0/1 treatment, covariate matrix, optional known propensity score, for the complete rows.
het_prepare <- function(data, outcome, treat, control, armTreat, covs, prop = NULL) {
    arm <- het_arm(data[[treat]], control)
    ctl <- levels(arm)[1]
    trt <- trimws(armTreat %||% "")
    if (!nzchar(trt)) trt <- levels(arm)[2]
    if (!trt %in% levels(arm) || identical(trt, ctl)) stop(sprintf("The treatment group value must be one of the non-control values (%s).", paste(levels(arm)[-1], collapse = ", ")), call. = FALSE)
    y <- het_num(data[[outcome]], "outcome")
    X <- het_covmat(data, covs)
    if (is.null(X) || ncol(X) == 0) stop("Select at least one covariate.", call. = FALSE)
    pr <- if (!is.null(prop)) het_num(data[[prop]], "propensity score") else NULL
    in2 <- arm %in% c(ctl, trt) & !is.na(arm)
    keep <- in2 & stats::complete.cases(y, X, if (!is.null(pr)) pr)
    list(y = y[keep], d = as.numeric(arm[keep] == trt), X = X[keep, , drop = FALSE], p = if (!is.null(pr)) pr[keep] else NULL, ctl = ctl, trt = trt,
         n_drop = sum(in2) - sum(keep), n = sum(keep), n1 = sum(arm[keep] == trt), n0 = sum(arm[keep] == ctl))
}

het_learner_names <- function(o) {
    map <- c(lrnLasso = "lasso", lrnRidge = "ridge", lrnEnet = "enet", lrnRf = "rf", lrnOls = "ols")
    unname(map[vapply(names(map), function(k) isTRUE(o[[k]]), TRUE)])
}

# Prepares the data and runs the procedure; returns the data summary and the result. K is the number of groups.
het_collect <- function(data, o, K, vcov = "hc3", monotonize = TRUE, cost = 0, checkpoint = NULL, rate = FALSE, nboot = 100, mods = NULL) {
    lrn <- het_learner_names(o)
    if (length(lrn) == 0) stop("Select at least one learner.", call. = FALSE)
    if (length(intersect(o$covs, c(o$treat, o$outcome, o$prop))) || identical(o$treat, o$outcome)) stop("A variable can sit in one box only.", call. = FALSE)
    pp <- het_prepare(data, o$outcome, o$treat, o$control, o$armTreat, o$covs, o$prop)
    if (pp$n1 < 20 || pp$n0 < 20) stop("Each group needs at least 20 units.", call. = FALSE)
    PM <- NULL
    if (!is.null(mods)) {
        bad <- setdiff(mods, o$covs)
        if (length(bad)) stop(sprintf("The moderators must be among the covariates (not in the covariates: %s).", paste(bad, collapse = ", ")), call. = FALSE)
        cols <- which(vapply(colnames(pp$X), function(cn) any(cn == mods | startsWith(cn, mods)), TRUE))
        PM <- pp$X[, cols, drop = FALSE]
        if (ncol(PM) == 0) stop("No moderator column was found.", call. = FALSE)
        PM <- PM[, apply(PM, 2, stats::sd) > 0, drop = FALSE]
    }
    res <- het_generic(pp$y, pp$d, pp$X, pp$X, pp$p, learners = lrn, K = K, nsplit = o$nsplit, share_main = o$shareMain, level = o$level / 100, vcov = vcov,
                       monotonize = monotonize, cost = cost, seed = o$seed, lambda_rule = if (identical(o$lambdaRule, "1se")) "lambda.1se" else "lambda.min", checkpoint = checkpoint, rate = rate, nboot = nboot, moderators = PM)
    list(pp = pp, res = res, learners = lrn)
}

het_row <- function(v, names) as.list(stats::setNames(as.numeric(v[names]), names))
