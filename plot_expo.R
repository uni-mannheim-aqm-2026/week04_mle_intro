loglik_exp_plot <-
  function (X,
            poss = NULL,
            plot.likfunc = TRUE,
            plot.density = TRUE,
            xlab = NULL,
            ylab = NULL,
            conv = 0.005,
            anim = TRUE,
            est.col = 2,
            density.leg = TRUE,
            cex.leg = 0.9,
            interval = 1,
            ...)
  {
    old.par <- par(no.readonly = TRUE)
    possibilities <- poss
    if (is.null(poss)) {
      possibilities <- seq(conv, 1, min(conv,
                                        0.01))
      logl <- possibilities
    }
    if (!is.null(poss))
      possibilities <- poss[poss != 0]
    logl <- possibilities
    for (i in 1:length(possibilities)) {
      logl[i] <- log(prod(dexp(X, possibilities[i])))
    }
    max.p <- ifelse(is.null(poss), mean(X), possibilities[logl ==
                                                            max(logl)][1])
    MLE <- 1 / mean(X)
    Max.lik <- log(prod(dexp(X, 1/mean(X))))
    if (anim == FALSE) {
      if (plot.likfunc == TRUE & plot.density == TRUE) {
        dev.new(height = 4, width = 8)
        par(
          mfrow = c(1, 2),
          mar = c(4.4, 4.5, 1, 0.5),
          cex = 0.9,
          las = 1,
          bty = 'n'
        )
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
      }
      if (plot.likfunc == TRUE) {
        plot(
          possibilities,
          logl,
          type = "l",
          ylab = ifelse(
            is.null(ylab),
            "Exponential Log-likelihood function",
            ylab
          ),
          xlab = ifelse(is.null(xlab), expression(paste(
            "Estimates for ",
            lambda
          )), xlab),
          las = 1,
          bty="n",
          ...
        )
        legend(
          "topright",
          legend = bquote(paste("ML est. = ",
                                .(MLE))),
          cex = cex.leg,
          bty = "n"
        )
        abline(v = mean(X),
               lty = 2,
               col = est.col)
      }
      if (plot.density == TRUE) {
        x <- NULL
        rm(x)
        curve(
          dexp(x, 1 / mean(X)),
          from = min(X),
          to = max(X),
          las = 1,
          bty="n",
          xlab = expression(italic(x)),
          ylab = expression(paste(italic(f),
                                  "(", italic(x), ")", sep = ""))
        )
        legend(
          "topright",
          legend = c(
            paste("X~EXP(1/", bquote(.(
              round(mean(X),
                    3)
            )), ")", sep = ""),
            paste("Loglik = ", bquote(.(
              round(Max.lik,
                    2)
            )))
          ),
          bty = "n",
          cex = cex.leg
        )
        if (density.leg == TRUE)
        {
          legend(
            "topleft",
            col = est.col,
            lty = 2,
            legend = "Obs. density",
            bty = "n",
            cex = cex.leg
          )
        }
        segments(X,
                 rep(0, length(X)),
                 X,
                 dexp(X, 1 / mean(X)),
                 lty = 2,
                 col = est.col)
      }
    }
    if (anim == TRUE) {
      nm <- which(logl == max(logl))[1]
      if (plot.likfunc == TRUE & plot.density == TRUE)
        dev.new(height = 4, width = 8)
      par(
        mfrow = c(1, 2),
        mar = c(4.4, 4.5, 1, 0.5),
        cex = 0.9
      )
      for (i in 1:(nm - 1)) {
        dev.hold()
        if (plot.likfunc == TRUE) {
          plot(
            possibilities,
            logl,
            type = "n",
            ylab = ifelse(
              is.null(ylab),
              
              "Exponential Log-likelihood function",
              ylab
            ),
            xlab = ifelse(is.null(xlab), expression(paste(
              "Estimates for ",
              lambda
            )), xlab),
            las = 1,
            bty="n",
            ...
          )
          arrows(
            possibilities[i],
            logl[i],
            possibilities[i +
                            1],
            logl[i + 1],
            col = est.col,
            length = 0.15,
            lwd = 1
          )
          points(
            possibilities[1:i],
            logl[1:i],
            lty = 2,
            col = est.col,
            lwd = 1,
            type = "l"
          )
          if (i == (nm - 1)) {
            points(possibilities, logl, type = "l")
            abline(v = 1 / mean(X),
                   lty = 2,
                   col = est.col)
            legend(
              "topright",
              legend = bquote(paste("ML est. = ",
                                    .(paste(round(MLE, 2))))),
              cex = cex.leg,
              bty = "n"
            )
          }
        }
        if (plot.density == TRUE) {
          x <- NULL
          rm(x)
          curve(
            dexp(x, possibilities[i]),
            from = min(X),
            to = max(X),
            xlab = expression(italic(x)),
            ylab = expression(paste(
              italic(f), "(", italic(x),
              ")", sep = ""
            )),
            las = 1,
            bty="n"
          )
          segments(
            x0 = X,
            y0 = rep(0, length(X)),
            x1 = X,
            y1 = dexp(X, possibilities[i]),
            col = est.col,
            lty = 2
          )
          if (i != (nm - 1)) {
            legend(
              "topright",
              legend = c(as.expression(bquote(
                paste(lambda,
                      " = ", .(round(
                        possibilities[i], 2
                      )))
              )),
              as.expression(bquote(
                paste("Loglik = ", .(round(logl[i],
                                           2)))
              ))),
              bty = "n",
              cex = cex.leg
            )
          }
          if (i == (nm - 1)) {
            legend(
              "topright",
              legend = c(as.expression(bquote(
                paste(lambda,
                      " = ", .(paste(round(MLE, 2))))
              )), as.expression(bquote(
                paste("Loglik = ",
                      .(round(
                        Max.lik, 2
                      )))
              ))),
              bty = "n",
              cex = cex.leg
            )
          }
        }
        dev.flush()
        Sys.sleep(interval)
      }
    }
    on.exit(par(old.par))
  }
