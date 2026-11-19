# jamovi analysis class: value of a treatment rule based on the predicted effects.

hetpolicyClass <- R6::R6Class(
    "hetpolicyClass",
    inherit = hetpolicyBase,
    private = list(

        .run = function() {
            o <- self$options
            het_show_text(self$results, o, "hetpolicy")
            if (is.null(o$outcome) || is.null(o$treat) || length(o$covs) == 0) {
                self$results$notes$setContent("<p>Select the outcome, the treatment variable and the covariates that may change the effect to begin.</p>")
                for (nm in c("fit", "val")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            out <- het_guard(het_collect(self$data, o, 2, "hc3", TRUE, o$cost, function() private$.checkpoint()))
            r <- out$res; pp <- out$pp; lrn <- out$learners
            lab <- het_learners[lrn]
            use <- r$best$gates
            fb <- self$results$fit
            for (i in seq_along(lrn))
                fb$addRow(rowKey = i, values = list(learner = unname(lab[i]), lam = r$lambda[lrn[i], "lambda"], lamb = r$lambda[lrn[i], "lambda_bar"], bg = if (identical(use, lrn[i])) "yes" else ""))
            fb$setVisible(length(lrn) > 1)
            pv <- r$policy[[use]]
            stat <- c("est", "se", "lo", "hi", "p")
            vb <- self$results$val
            what <- c(gain_none = "Gain from the rule compared with treating nobody", gain_all = "Gain from the rule compared with treating everybody", ate = "Average effect of the treatment (for reference)")
            for (i in seq_along(what)) vb$addRow(rowKey = i, values = c(list(what = unname(what[i])), het_row(pv[names(what)[i], ], stat)))
            vb$addRow(rowKey = 4, values = list(what = "Share of units the rule treats", est = pv["share", "est"], se = NaN, lo = NaN, hi = NaN, p = NaN))
            vb$setNote("a", sprintf("Learner: %s. Value in outcome units per unit of the population, net of the cost of %s per treated unit. The rule is learned on one part of the sample and evaluated with doubly robust scores on the other part; medians over %d splits.", unname(lab[use]), het_fmt(o$cost, 3), o$nsplit))
            notes <- sprintf("%d units (%d treated with '%s', %d control with '%s').", pp$n, pp$n1, pp$trt, pp$n0, pp$ctl)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (pv["gain_all", "est"] < 0 && pv["gain_none", "est"] > 0) notes <- c(notes, "Treating everybody would be worse than treating only the units with a large predicted effect.")
            self$results$notes$setContent(paste0("<p>", paste(notes, collapse = " "), "</p>"))
        }
    )
)
