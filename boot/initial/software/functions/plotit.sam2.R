#-*- coding: utf-8 -*-

### File: cdd_survey_scripts.R
### Time-stamp: <2024-04-25 12:52:51 a23579>
###
### Created: 2016??
### Author: ?? + Jennifer Devine, Yves Reecht
###
####################################################################################################
### Description:
###
### Misc plot functions
####################################################################################################

plotit.sam2 <- function(fit, what, x = fit$data$years, ylab = what, xlab = "Years",
                        ex = numeric(0), trans = function(x)x, add = FALSE, grid = TRUE, gridcol = 1,
                        valcol = 1, vallty = 1, vallwd = 1,
                        ci = TRUE, ciborder = gray(.5, alpha = .5), cicol = gray(.5, alpha = .5), cilty = 1, cilwd = 1,
                        las = las, addCI = NA, drop = 0, unnamed.basename = "current", xlim = NULL, ...)
{
    idx <- names(fit$sdrep$value) == what
    y <- fit$sdrep$value[idx]
    lowhig <- y + fit$sdrep$sd[idx]%o%c(-2, 2)
    didx <- 1:(length(x) - drop)
    if(missing(xlim)){
        xr <- range(x)
    }else{
        xr <- xlim
    }
    x <- x[didx]
    y <- y[didx]
    lowhig <- lowhig[didx, ]
    if(!add){
        plot(x, trans(y), xlab = xlab, ylab = ylab, type = "n", xlim = xr,
             ylim = range(c(trans(lowhig), 0,ex)), ...)
        if(grid) grid(col = gridcol)
    }
    if(ci){
        polygon(c(x, rev(x)), y = c(trans(lowhig[, 1]), rev(trans(lowhig[, 2]))),
                border = NA, col = cicol)
        lines(x = x, y = trans(lowhig[, 1]), col = ciborder, lty = cilty, lwd = cilwd)
        lines(x = x, y = trans(lowhig[, 2]), col = ciborder, lty = cilty, lwd = cilwd)
    }
    lines(x, trans(y), col = valcol, lty = vallty, lwd = vallwd, ...)
}





addforecast2 <- function(fit, what,
                         dotcol = "black", dotpch = 19, dotcex = 1,
                         arrowcol = gray(.5, alpha = .5), arrowlwd = 1,
                         arrowangle = 90, arrowcode = 3, arrowlength = 0.1, arrowlty = 1,
                         ...)
{

    x <- attr(fit, "tab")
    y <- as.numeric(rownames(x))
    dummy <- sapply(1:length(y),
                    function(i)
             {
                 xx <- c(x[i, paste(what, "low", sep = ":")], x[i, paste(what, "high", sep = ":")])
                 units = par(c('usr', 'pin'))
                 xx_to_inches = with(units, pin[2L] / diff(usr[3:4]))
                 if(abs(xx_to_inches * diff(xx)) > 0.01){
                     arrows(y[i], xx[1], y[i], xx[2], lwd = arrowlwd, col = arrowcol,
                            angle = arrowangle, code = arrowcode, length = arrowlength,
                            lty = arrowlty, ...)
                 }
             })
    points(y, x[, paste(what, "median", sep = ":")], pch = dotpch, cex = dotcex, col = dotcol)
}

### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
