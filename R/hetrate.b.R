# jamovi analysis class: priority curve (rank-weighted average treatment effects, AUTOC, Qini) of a ranking by predicted effect.

hetrateClass <- R6::R6Class(
    "hetrateClass",
    inherit = hetrateBase,
    private = list(

        .run = function() {
            o <- self$options
            het_show_text(self$results, o, "hetrate")
            if (is.null(o$outcome) || is.null(o$treat) || length(o$covs) == 0) {
                self$results$notes$setContent("<p>Select the outcome, the treatment variable and the covariates that may change the effect to begin.</p>")
                for (nm in c("fit", "rate", "toc")) self$results[[nm]]$setVisible(FALSE)
                return()
            }
            out <- het_guard(het_collect(self$data, o, 2, "hc3", TRUE, 0, function() private$.checkpoint(), rate = TRUE, nboot = o$nboot))
            r <- out$res; pp <- out$pp; lrn <- out$learners
            lab <- het_learners[lrn]
            use <- r$best$gates
            fb <- self$results$fit
            for (i in seq_along(lrn))
                fb$addRow(rowKey = i, values = list(learner = unname(lab[i]), lam = r$lambda[lrn[i], "lambda"], lamb = r$lambda[lrn[i], "lambda_bar"], bg = if (identical(use, lrn[i])) "yes" else ""))
            fb$setVisible(length(lrn) > 1)
            rt <- r$rate[[use]]; stat <- c("est", "se", "lo", "hi", "p")
            rb <- self$results$rate
            rb$addRow(rowKey = 1, values = c(list(what = "AUTOC (area under the TOC curve)"), het_row(rt["AUTOC", ], stat)))
            rb$addRow(rowKey = 2, values = c(list(what = "Qini (weights the top ranks more)"), het_row(rt["QINI", ], stat)))
            rb$setNote("a", sprintf("Learner: %s. A positive value with a small p-value means that units ranked high by the predicted effect gain more than the average unit, so the ranking is useful for deciding whom to treat first; zero means the ranking is no better than a random order. The units are those of the estimation part of each split; medians over %d splits.", unname(lab[use]), o$nsplit))
            q <- r$qgrid; show <- which(abs(q * 10 - round(q * 10)) < 1e-9)
            tb <- self$results$toc
            for (j in seq_along(show)) {
                qq <- q[show[j]]
                tb$addRow(rowKey = j, values = if (abs(qq - 1) < 1e-9) list(q = qq, est = 0, se = NaN, lo = NaN, hi = NaN, p = NaN) else c(list(q = qq), het_row(rt[paste0("q", qq), ], stat)))
            }
            tb$setNote("a", "TOC(q): the average effect of the top share q of the units (ranked by predicted effect) minus the average effect of all units. It is zero when everybody is included (q = 1).")
            if (isTRUE(o$plot)) {
                rows <- paste0("q", q)
                self$results$tocPlot$setState(list(q = q, est = unname(rt[rows, "est"]), lo = unname(rt[rows, "lo"]), hi = unname(rt[rows, "hi"]), autoc = unname(rt["AUTOC", "est"])))
            }
            notes <- sprintf("%d units (%d treated with '%s', %d control with '%s'); %d random splits, %d bootstrap samples each for the standard errors; estimates are medians over the splits, and p-values and intervals follow the adjustments of Chernozhukov et al. (2018).", pp$n, pp$n1, pp$trt, pp$n0, pp$ctl, o$nsplit, o$nboot)
            if (pp$n_drop > 0) notes <- c(notes, sprintf("%d observation(s) with missing values were dropped.", pp$n_drop))
            if (is.null(pp$p)) notes <- c(notes, sprintf("The probability of treatment is taken to be the share treated, %s.", het_fmt(mean(pp$d), 3)))
            self$results$notes$setContent(paste0("<p>", paste(notes, collapse = " "), "</p>"))
        },

        .plotToc = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            df <- data.frame(q = st$q, est = st$est, lo = st$lo, hi = st$hi)
            g <- ggplot2::ggplot(df, ggplot2::aes(x = q, y = est)) + ggtheme +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "#6b9de8", alpha = 0.2) +
                ggplot2::geom_line(colour = "#6b9de8", linewidth = 1) +
                ggplot2::geom_point(colour = "#6b9de8", size = 1.8) +
                ggplot2::labs(x = "Share of units treated first (highest predicted effect)", y = "TOC: gain over the average effect", caption = sprintf("AUTOC %s; dashed: ranking no better than random", signif(st$autoc, 3)))
            print(g)
            TRUE
        }
    )
)
