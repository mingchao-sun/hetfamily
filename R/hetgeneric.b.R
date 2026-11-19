# jamovi analysis class: generic machine learning inference on heterogeneous effects (BLP, GATES, CLAN).

hetgenericClass <- R6::R6Class(
    "hetgenericClass",
    inherit = hetgenericBase,
    private = list(

        .run = function() {
            o <- self$options
            het_show_text(self$results, o, "hetgeneric")
            if (is.null(o$outcome) || is.null(o$treat) || length(o$covs) == 0) {
                self$results$notes$setContent("<p>Select the outcome, the treatment variable and the covariates that may change the effect to begin.</p>")
                for (nm in c("fit", "blp", "gates", "clan")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            out <- het_guard(het_collect(self$data, o, o$groups, o$vcov, isTRUE(o$monotonize), 0, function() private$.checkpoint()))
            r <- out$res; pp <- out$pp; lrn <- out$learners; K <- r$K
            lab <- het_learners[lrn]
            fb <- self$results$fit
            for (i in seq_along(lrn))
                fb$addRow(rowKey = i, values = list(learner = unname(lab[i]), lam = r$lambda[lrn[i], "lambda"], lamb = r$lambda[lrn[i], "lambda_bar"],
                                                    bb = if (identical(r$best$blp, lrn[i])) "yes" else "", bg = if (identical(r$best$gates, lrn[i])) "yes" else ""))
            fb$setNote("a", "Lambda = (slope)^2 x variance of the predicted effects: how much of the heterogeneity the learner captures in the best linear predictor; lambda-bar = the same for the group effects. The learner with the largest value is used below (medians over the sample splits).")
            fb$setVisible(length(lrn) > 1)
            stat <- c("est", "se", "lo", "hi", "p")
            # BLP
            lb <- r$best$blp; bt <- r$blp[[lb]]
            bb <- self$results$blp
            bb$addRow(rowKey = 1, values = c(list(what = "Average effect of the treatment"), het_row(bt["ATE", ], stat)))
            bb$addRow(rowKey = 2, values = c(list(what = "Heterogeneity (slope on the predicted effect)"), het_row(bt["HET", ], stat)))
            bb$setNote("a", sprintf("Learner: %s. The slope is 1 if the predicted effects are right on average and 0 if they carry no information about who gains; the p-value tests no heterogeneity that the learner can predict. A large p-value does not prove that the effect is constant.", lab[[lb]]))
            # GATES
            lg <- r$best$gates; gt <- r$gates[[lg]]; gd <- r$gates_diff[[lg]]
            gb <- self$results$gates
            for (k in seq_len(K))
                gb$addRow(rowKey = k, values = c(list(what = sprintf("Group %d%s", k, if (k == 1) " (least affected)" else if (k == K) " (most affected)" else "")), het_row(gt[k, ], stat)))
            gb$addRow(rowKey = K + 1, values = c(list(what = sprintf("Group %d minus group 1", K)), het_row(gd[1, ], stat)))
            gb$setNote("a", sprintf("Learner: %s. Units are sorted into %d groups by the predicted effect.%s", lab[[lg]], K, if (isTRUE(o$monotonize)) " The estimates are sorted so that they increase from group 1 to group K." else ""))
            # CLAN
            if (isTRUE(o$clan)) {
                cl <- r$clan[[lg]]; cb <- self$results$clan
                for (j in seq_len(nrow(cl)))
                    cb$addRow(rowKey = j, values = list(what = rownames(cl)[j], m1 = cl[j, "mean_low"], m2 = cl[j, "mean_high"], est = cl[j, "diff"], se = cl[j, "se"], lo = cl[j, "lo"], hi = cl[j, "hi"], p = cl[j, "p"]))
                cb$setNote("a", sprintf("Means of the covariates in the least affected (group 1) and most affected (group %d) group, learner %s. Nominal variables appear as 0/1 indicators. Many covariates are compared, so single p-values need care.%s", K, lab[[lg]], if (bt["HET", "p"] > 0.05) " The heterogeneity test does not reject, so the groups may reflect noise in the predictions: read this table with caution." else ""))
            }
            if (isTRUE(o$plot)) self$results$gatesPlot$setState(list(est = unname(gt[, "est"]), lo = unname(gt[, "lo"]), hi = unname(gt[, "hi"]), ate = unname(bt["ATE", c("est", "lo", "hi")]), K = K))
            notes <- sprintf("%d units (%d treated with '%s', %d control with '%s'); %d random splits, each estimating on %d%% of the sample and learning on the rest; estimates are medians over the splits, and p-values and intervals follow the adjustments of Chernozhukov et al. (2018).",
                             pp$n, pp$n1, pp$trt, pp$n0, pp$ctl, o$nsplit, round(100 * o$shareMain))
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (is.null(pp$p)) notes <- c(notes, sprintf("The probability of treatment is taken to be the share treated, %s.", het_fmt(mean(pp$d), 3)))
            self$results$notes$setContent(paste0("<p>", paste(notes, collapse = " "), "</p>"))
        },

        .plotGates = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(g = seq_len(st$K), est = st$est, lo = st$lo, hi = st$hi)
            g <- ggplot2::ggplot(df, ggplot2::aes(x = g, y = est)) + ggtheme +
                ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = st$ate[2], ymax = st$ate[3], fill = "#e6ac40", alpha = 0.18) +
                ggplot2::geom_hline(yintercept = st$ate[1], colour = "#e6ac40", linewidth = 0.8) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
                ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.15, colour = "#6b9de8", linewidth = 0.8) +
                ggplot2::geom_point(size = 3, colour = "#6b9de8") +
                ggplot2::scale_x_continuous(breaks = seq_len(st$K)) +
                ggplot2::labs(x = "Group of predicted effect (1 = lowest)", y = "Average treatment effect", caption = "Orange: average effect with its interval; dashed: no effect")
            print(g)
            TRUE
        }
    )
)
