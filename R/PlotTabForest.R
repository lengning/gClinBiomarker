#' @keywords internal

PlotTabForest <- function(label.text,
                          mean,
                          lower,
                          upper,
                          headings,
                          cols=4,
                          xticks=NULL,
                          box.size=NULL,
                          main="Forest Plot",
                          sub.main=c("Treatment better", "Control better"),
                          sub="Subtitle",
                          hline,
                          group.hline,
                          vline,
                          note,
                          clip=c(-Inf, Inf),
                          widths,
                          heights,
                          cex.headings=1,
                          cex.note=0.8,
                          xlog=FALSE, cex=1, cex.main=2, cex.sub=1.3, cex.axis=1, pch=19,
                          par.param=NULL){


    n.row = dim(label.text)[1]
    n.col = dim(label.text)[2] + 1

    if (length(cols) == 1) {
        cols <- rep(cols, n.row)
    }

    if (missing(headings)) {
        headings = rep("", n.col)
    } else {
        headings = rep(headings, length.out=n.col)
    }

    if (missing(widths)) {
        widths = c(lcm(3), rep(lcm(2.5), n.col-2), lcm(7))
    }

    if (missing(heights)) {
        heights = c(lcm(2.5), lcm(1.5), rep(lcm(1), n.row), lcm(1.5))
    } else {
        heights = c(lcm(2.5), heights, lcm(1.5))
    }

    tmp = matrix(seq(n.row*n.col)+n.col+1, ncol=n.col, byrow=F)
    tmp2 = rbind(rep(1, n.col), seq(n.col)+1, tmp, rep((n.row+1)*n.col+2, n.col))

    layout(tmp2, widths=widths, heights=heights)

    par(mar=c(0, 0, 0, 0))

    # Write main title
    plot.new()
    plot.window(xlim=c(0, 1), ylim=c(0, 1))
    text(0.5, 0.5, adj=c(0.5, 0.5), lab=main, font=2, cex=cex.main) #centered, bold-face

    # Write headings
    for(jj in c(1:n.col)) {
        plot.new()
        plot.window(xlim=c(0, 1), ylim=c(0, 1))

        if (jj == 1) {
            text(0, 0.9, adj=c(0, 1), lab=paste(headings[jj]), font=2, cex=cex.headings) #left-justified, bold-face
        } else {
            text(0.5, 0.9, adj=c(0.5, 1), lab=paste(headings[jj]), font=2, cex=cex.headings) #centered, bold-face
        }

        lines(c(-0.5, 1.5), c(1, 1))
        lines(c(-0.5, 1.5), c(0, 0))
        lines(c(-0.5, 1.5), c(0+0.05, 0+0.05))
    }

    # Display the table
    for(jj in c(1:(n.col-1))) {
        for(ii in c(1:n.row)) {
            plot.new()
            plot.window(xlim=c(0, 1), ylim=c(0, 1))

            if (jj==1) {
                text(0, 0.5, adj=c(0, 0.5), lab=paste(label.text[ii, 1]), cex=cex) #left-justified
            } else {
                text(0.5, 0.5, adj=c(0.5, 0.5), lab=paste(label.text[ii, jj]), cex=cex) #centered, bold-face
            }

            if (ii %in% hline) {
                lines(c(-0.5, 1.5), c(0, 0), col=8)
            }

            if (ii %in% group.hline) {
                lines(c(-0.5, 1.5), c(0, 0), col=1)
            }
        }
    }

    l <- vector(mode="numeric", length=length(lower))
    u <- vector(mode="numeric", length=length(upper))

    # Display the forest plot
    if (xlog) {
        mean = log(mean)
        lower = log(lower)
        upper = log(upper)

        for (i in 1:length(lower)) {
            if (!is.na(lower[i])) {
                if (abs(lower[i]) == Inf) {
                    lower[i] <- min(lower[is.finite(lower)])
                    l[i] <- 1
                }
            }
        }

        for (i in 1:length(upper)) {
            if (!is.na(upper[i])) {
                if (abs(upper[i]) == Inf) {
                    upper[i] <- max(upper[is.finite(upper)])
                    u[i] <- 1
                }
            }
        }

        clip = ifelse(is.finite(clip), log(clip), clip)
        vline = log(vline)
    }

    par(mar=c(0, 0, 0, 0))

    if (is.null(xticks)) {
        xmax = max(upper, na.rm=T)
        xmin = min(lower, na.rm=T)

        if(is.finite(clip[2])) {
            xmax = clip[2]
        }

        if(is.finite(clip[1])) {
            xmin = clip[1]
        }

    } else {
        if (xlog) {
            xmax = max(log(xticks))
            xmin= min(log(xticks))
        } else {
            xmax = max(xticks)
            xmin = min(xticks)
        }
    }

    for (ii in c(1:length(mean))) {
        plot.new()
        plot.window(xlim=c(xmin, xmax), ylim=c(0,1))

        if (ii == 1 & !is.null(sub.main)) {
            mtext(sub.main[1], side=3, at=xmin+(xmax-xmin)/4, line=1, font=2, cex=0.7)
            mtext(sub.main[2], side=3, at=xmin+3*(xmax-xmin)/4, line=1, font=2, cex=0.7)
        }

        if (!is.na(lower[ii]) & !is.na(upper[ii])) {
            if (lower[ii] < xmin & upper[ii] > xmax) {
                arrows(xmin, 0.5, xmax, 0.5, code=3, length=0.1, col=cols[ii], pch=pch)
            }

            if (lower[ii] < xmin & upper[ii] <= xmax) {
                arrows(xmin, 0.5, upper[ii], 0.5, code=1, length=0.1, col=cols[ii], pch=pch)
            }

            if (lower[ii] >= xmin & upper[ii] > xmax) {
                arrows(lower[ii], 0.5, xmax, 0.5, code=2, length=0.1, col=cols[ii], pch=pch)
            }

            if (lower[ii] >= xmin & upper[ii] <= xmax) {
                lines(c(lower[ii], upper[ii]), c(0.5, 0.5), col=cols[ii], pch=pch)
            }

            if (l[ii] == 1 & u[ii] == 1) {
                arrows(xmin, 0.5, xmax, 0.5, code=3, length=0.1, col=cols[ii], pch=pch)
            }

            if (l[ii] == 1 & u[ii] != 1) {
                arrows(xmin, 0.5, upper[ii], 0.5, code=1, length=0.1, col=cols[ii], pch=pch)
            }

            if (u[ii] == 1 & l[ii] != 1) {
                arrows(lower[ii], 0.5, xmax, 0.5, code=2, length=0.1, col=cols[ii], pch=pch)
            }
        }

        if (!is.na(mean[ii])) {
            points(mean[ii], 0.5, col=cols[ii], cex=box.size[ii], pch=pch)
        }

        if (ii %in% hline) {
            lines(c(xmin-10, xmax+10), c(0, 0), col=8)
        }

        if (ii %in% group.hline) {
            lines(c(xmin-10, xmax+10), c(0, 0), col=1)
        }

        abline(v=vline, lty=2)
    }

    # Draw axis
    if (xlog) {
        if (is.null(xticks)) {
            xticks = round(exp(c(0, xmin, xmax, pretty(range(xmin, xmax)))), 1)
        }
        axis(side=1, at=log(xticks), labels=xticks, font=2, lwd=2, cex.axis=cex.axis)
    } else {
        axis(side=1, at=xticks, labels=xticks, font=2, lwd=2, cex.axis=cex.axis)
    }

    # Write subtitle
    plot.new()
    plot.window(xlim=c(0, 1), ylim=c(0, 1))
    text(0.5, 0.5, adj=c(1, 0), lab=sub, font=2, cex=cex.sub) #centered, bold-face
    text(-0.04, 1, adj=c(0, 1), lab=note, cex=cex.note) #left-justified
}
