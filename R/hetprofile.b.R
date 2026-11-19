# jamovi analysis class: effect profile, the linear projection of the conditional effect on the moderators (regression of doubly robust scores).

hetprofileClass <- R6::R6Class(
    "hetprofileClass",
    inherit = hetprofileBase,
    private = list(

        .run = function() {
            o <- self$options
            het_show_text(self$results, o, "hetprofile")
            if (is.null(o$outcome) || is.null(o$treat) || length(o$covs) == 0) {
                self$results$notes$setContent("<p>Select the outcome, the treatment variable and the covariates that may change the effect to begin. Moderators default to all covariates.</p>")
                for (nm in c("fit", "prof")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            mods <- if (length(o$moderators) > 0) o$moderators else o$covs
            out <- het_guard(het_collect(self$data, o, 2, o$vcov, TRUE, 0, function() private$.checkpoint(), mods = mods))
            r <- out$res; pp <- out$pp; lrn <- out$learners
            lab <- het_learners[lrn]
            use <- r$best$gates
            fb <- self$results$fit
            for (i in seq_along(lrn))
                fb$addRow(rowKey = i, values = list(learner = unname(lab[i]), lam = r$lambda[lrn[i], "lambda"], lamb = r$lambda[lrn[i], "lambda_bar"], bg = if (identical(use, lrn[i])) "yes" else ""))
            fb$setVisible(length(lrn) > 1)
            pf <- r$profile[[use]]; stat <- c("est", "se", "lo", "hi", "p")
            pb <- self$results$prof
            for (i in seq_len(nrow(pf)))
                pb$addRow(rowKey = i, values = c(list(what = if (i == 1) "Effect at the average of the moderators" else paste("Change per standard deviation of", rownames(pf)[i])), het_row(pf[i, ], stat)))
            pb$setNote("a", sprintf("Learner: %s. The moderators are standardised (mean 0, standard deviation 1), so the first row is the average effect and each other row is how much the effect changes when that moderator rises by one standard deviation, the others held fixed in this linear approximation. It is a description of how the effect varies, not a causal effect of the moderator.", unname(lab[use])))
            if (isTRUE(o$plot) && nrow(pf) > 1) self$results$profPlot$setState(list(label = rownames(pf)[-1], est = unname(pf[-1, "est"]), lo = unname(pf[-1, "lo"]), hi = unname(pf[-1, "hi"])))
            notes <- sprintf("%d units (%d treated with '%s', %d control with '%s'); %d random splits, each estimating on %d%% of the sample and learning on the rest; the regression of the doubly robust scores on the moderators is done on the estimation part and combined by medians, with the adjustments of Chernozhukov et al. (2018).", pp$n, pp$n1, pp$trt, pp$n0, pp$ctl, o$nsplit, round(100 * o$shareMain))
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (nrow(pf) > 6) notes <- c(notes, "Many moderators are compared, so single p-values need care.")
            self$results$notes$setContent(paste0("<p>", paste(notes, collapse = " "), "</p>"))
        },

        .plotProf = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(label = factor(st$label, levels = rev(st$label)), est = st$est, lo = st$lo, hi = st$hi)
            g <- ggplot2::ggplot(df, ggplot2::aes(x = est, y = label)) + ggtheme +
                ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
                ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2, colour = "#6b9de8", linewidth = 0.8) +
                ggplot2::geom_point(size = 2.8, colour = "#6b9de8") +
                ggplot2::labs(x = "Change of the effect per standard deviation of the moderator", y = NULL)
            print(g)
            TRUE
        }
    )
)
