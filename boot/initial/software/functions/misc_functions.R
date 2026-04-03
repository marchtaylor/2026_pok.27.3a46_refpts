#-*- coding: utf-8 -*-

### File: misc_functions.R
### Time-stamp: <2025-04-28 10:32:32 a23579>
###
### Created: 25/04/2022	09:50:17
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################



## Geometric mean:
gmean <- function(x, na.rm = TRUE, zero.propagate = FALSE)
{
    ##
    if(any(x < 0, na.rm = TRUE)){
        return(NaN)
    }
    if(zero.propagate){
        if(any(x == 0, na.rm = TRUE)){
            return(0)
        }
        exp(mean(log(x), na.rm = na.rm))
    } else {
        exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
    }
}

## Functions for SRR:
extractSRinfo <- function (fit)
{
    require(tidyr)
    require(dplyr)
    require(tibble)
    X <- summary(fit)
    n <- nrow(X)
    head(X)
    lag <- fit$conf$minAge
    idxR <- (lag + 1):n
    idxS <- 1:(n - lag)
    R <- X[idxR, 1]
    S <- X[idxS, 4]
    Rnam <- colnames(X)[1]
    Snam <- colnames(X)[4]
    y <- rownames(X)

    return(tibble(Year = y[idxR],
                  R, S) %>%
           setNames(c("Year", Rnam, Snam)))
}

extractFleetInfo <- function (fit, log = TRUE,
                              fleets = unique(fit$data$aux[, "fleet"]),
                              ...)
{
    require(tidyr)
    require(dplyr)
    require(tibble)
    idx <- fit$data$aux[, "fleet"] %in% fleets
    trans <- function(x)
    {
        if (log) {
            x
        }else{
            exp(x)
        }
    }
    p <- trans(fit$obj$report(c(fit$sdrep$par.fixed,
                                fit$sdrep$par.random))$predObs[idx])
    o <- trans(fit$data$logobs[idx])
    aa <- fit$data$aux[idx, "age"]
    neg.age <- (aa < -1e-06)
    aa[neg.age] <- NA
    a <- paste0("a = ", aa, " ")
    f <- paste0("", strtrim(attr(fit$data,
                                 "fleetNames")[fit$data$aux[idx, "fleet"]],
                            50))
    Year <- fit$data$aux[idx, "year"]

    res <- tibble(fleet = f,
                  year = Year,
                  age = aa,
                  pred = p,
                  obs = o)
    return(res)
}

plot.SAMres <- function(x, fleet = "all")
{
    ## Purpose: Plot OSA residuals from a sam residual object.
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 25 Apr 2022, 18:43
    require(tidyr)
    require(dplyr)
    require(tibble)
    require(ggplot2)

    fn <- attr(x, "fleetNames")

    if (is.character(fleet))
    {
        if (fleet[1] != "all")
        {
            fl <- seq_along(fn)
        }else{
            fl <- which(fn %in% fleet)
        }
    }

    if (is.numeric(fleet))
    {
        fl <- fleet[fleet %in% seq_along(fn)]
    }

    if ((! is.numeric(fleet) && ! is.character(fleet)) ||
        length(fl) == 0)
    {
        fl <- seq_along(fn)
    }

    df <- x %>% as.list.data.frame() %>% as_tibble() %>%
        mutate(fleet = factor(fn[fleet], levels = fn)) %>%
        filter(fleet %in% fn[fl])

    ggplot(data = df,
           aes(x = year, y = residual, color = residual > 0)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_point() +
        scale_color_discrete(name = "Sign", labels = c("-", "+")) +
        facet_grid(age~fleet, labeller = labeller(age = label_both, fleet = label_value))
}


get_watage <- function(x, type = c("stock", "catch", "landing", "discard"),
                       nAveYears = 3)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 24 Apr 2023, 17:57

    ## Extract fit object:
    fittmp <- attr(x, "fit")

    fitYears <- fittmp$data$years
    fcYears <- sapply(x, function(x) x$year)

    ## Type of estimate:
    estType <- data.frame(Year = as.character(c(fitYears, tail(fcYears, -1))),
                          EstType = c(rep("Estimate", length(fitYears)),
                                      "STF, intermediate year",
                                      "STF, advice year",
                                      "STF"))

    type <- match.arg(tolower(type),
                      c("stock", "catch", "landing", "discard"))

    var <- case_when(type %in% "stock" ~ "stockMeanWeight",
                     type %in% "catch" ~ "catchMeanWeight",
                     type %in% "landing" ~ "landMeanWeight",
                     type %in% "discard" ~ "disMeanWeight")

    ## Try to guess the ages in the model:
    minAge <- fittmp$conf$minAge
    maxAge <- fittmp$conf$maxAge

    maxAgePlus <- max(fittmp$conf$maxAgePlusGroup *
                      fittmp$data$maxAgePerFleet) == maxAge

    if(maxAgePlus)
    {
        ages <- c(as.character(seq(minAge, maxAge - 1)),
                  paste0(maxAge, "+"))
    }else{
        ages <- as.character(seq(minAge, maxAge))
    }

    allYears <- as.numeric(row.names(fittmp$data[[var]]))

    ## Weights at age and mean used for forecast:
    wgtmt <- fittmp$data[[var]]
    if (length(dim(wgtmt)) == 3)
    {
        wgtmt <- as.matrix(wgtmt[ , , 1])
    }
    wgt <- do.call(rbind,
                   c(list(wgtmt),
                     rep(list(apply(tail(wgtmt, nAveYears), 2, mean)),
                         times = sum(! fcYears %in% allYears)))) # Do not add averages for years provided.

    ## Fancier and more informative outputs:
    row.names(wgt) <- c(head(row.names(wgt), -(length(fcYears) - 1)),
                        tail(fcYears, -1))

    colnames(wgt) <- ages

    wgt <- cbind(data.frame(Age = ages),
                 t(wgt))

    attr(wgt, "weightType") <- type
    attr(wgt, "estType") <- estType

    return(wgt)
}

get_moatage <- function(x,
                        nAveYears = 3)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 24 Apr 2023, 17:57

    ## Extract fit object:
    fittmp <- attr(x, "fit")

    fitYears <- fittmp$data$years
    fcYears <- sapply(x, function(x) x$year)

    ## Type of estimate:
    estType <- data.frame(Year = as.character(c(fitYears, tail(fcYears, -1))),
                          EstType = c(rep("Estimate", length(fitYears)),
                                      "STF, intermediate year",
                                      "STF, advice year",
                                      "STF"))

    ## Try to guess the ages in the model:
    minAge <- fittmp$conf$minAge
    maxAge <- fittmp$conf$maxAge

    maxAgePlus <- max(fittmp$conf$maxAgePlusGroup *
                      fittmp$data$maxAgePerFleet) == maxAge

    if(maxAgePlus)
    {
        ages <- c(as.character(seq(minAge, maxAge - 1)),
                  paste0(maxAge, "+"))
    }else{
        ages <- as.character(seq(minAge, maxAge))
    }

    names(fittmp$data)

    allYears <- as.numeric(row.names(fittmp$data[["propMat"]]))

    ## Weights at age and mean used for forecast:
    mo <- do.call(rbind,
                  c(list(fittmp$data[["propMat"]]),
                    rep(list(apply(tail(fittmp$data[["propMat"]], nAveYears), 2, mean)),
                        times = sum(! fcYears %in% allYears)))) # Do not add averages for years provided.

    ## Fancier and more informative outputs:
    row.names(mo) <- c(head(row.names(mo), -(length(fcYears) - 1)),
                        tail(fcYears, -1))

    colnames(mo) <- ages

    mo <- cbind(data.frame(Age = ages),
                 t(mo))

    attr(mo, "estType") <- estType

    return(mo)
}

get_natage <- function(x, funFC = median, funFCrec = gmean)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 24 Apr 2023, 17:12

    ## Extract fit object and info:
    fittmp <- attr(x, "fit")

    fitYears <- fittmp$data$years
    fcYears <- sapply(x, function(x) x$year)

    ## Type of estimate:
    estType <- data.frame(Year = as.character(c(fitYears, tail(fcYears, -1))),
                          EstType = c(rep("Estimate", length(fitYears)),
                                      "STF, intermediate year",
                                      "STF, advice year",
                                      "STF"))

    ## Try to guess ages in the model:
    minAge <- fittmp$conf$minAge
    maxAge <- fittmp$conf$maxAge

    maxAgePlus <- max(fittmp$conf$maxAgePlusGroup *
                      fittmp$data$maxAgePerFleet) == maxAge

    if(maxAgePlus)
    {
        ages <- c(as.character(seq(minAge, maxAge - 1)),
                  paste0(maxAge, "+"))
    }else{
        ages <- as.character(seq(minAge, maxAge))
    }

    ## Assembling (median) estimates and forecasts:
    res <- cbind(data.frame(Age = ages),
                 exp(fittmp$pl$logN),
                 ## What follows can be made more flexible later...
                 FC2 = c(funFCrec(exp(x[[2]]$sim[, 1])),
                         apply(x[[2]]$sim[, tail(seq_along(ages), -1)], 2,
                               function(x) funFC(exp(x)))),
                 FC3 = c(funFCrec(exp(x[[3]]$sim[, 1])),
                         apply(x[[3]]$sim[, tail(seq_along(ages), -1)], 2,
                               function(x) funFC(exp(x)))),
                 FC4 = c(funFCrec(exp(x[[4]]$sim[, 1])),
                         apply(x[[4]]$sim[, tail(seq_along(ages), -1)], 2,
                               function(x) funFC(exp(x)))))

    colnames(res) <- c("Age", fitYears, tail(fcYears, -1))

    attr(res, "estType") <- estType

    return(res)
}

get_fatage <- function(x, funFC = median, nAveYears = 3)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 25 Apr 2023, 12:10

    ## Extract fit object and info:
    fittmp <- attr(x, "fit")

    fitYears <- fittmp$data$years
    fcYears <- sapply(x, function(x) x$year)

    ## Type of estimate:
    estType <- data.frame(Year = as.character(c(fitYears, tail(fcYears, -1))),
                          EstType = c(rep("Estimate", length(fitYears)),
                                      "STF, intermediate year",
                                      "STF, advice year",
                                      "STF"))

    ## Try to guess ages in the model:
    minAge <- fittmp$conf$minAge
    maxAge <- fittmp$conf$maxAge

    maxAgePlus <- max(fittmp$conf$maxAgePlusGroup *
                      fittmp$data$maxAgePerFleet) == maxAge

    if(maxAgePlus)
    {
        ages <- c(as.character(seq(minAge, maxAge - 1)),
                  paste0(maxAge, "+"))
    }else{
        ages <- as.character(seq(minAge, maxAge))
    }

    agesT <- ages[ ! duplicated(fittmp$conf$keyLogFsta[1, ])]

    ## Average exploitation pattern ():
    Fave <- apply(exp(fittmp$pl$logF[, (length(fitYears) - nAveYears + 1):length(fitYears)]),
                  1, mean)

    ## Row indices corresponding to Fbar:
    fbarAi <- which(dplyr::between(as.numeric(agesT),
                                   fittmp$conf$fbarRange[1],
                                   fittmp$conf$fbarRange[2]))

    ## Fbar for average exploitation pattern:
    FbarAve <- mean(Fave[fbarAi])


    ## Equivalent to scale based on the me(di)an of simulated Fbar:
    ## Fave * median(x[[2]]$fbar) / FbarAve
    ## apply(Fave %*% t(x[[2]]$fbar / FbarAve),
    ##       1, median)

    ## Assembling (median) estimates and forecasts:
    res <- cbind(data.frame(Age = agesT),
                 exp(fittmp$pl$logF),
                 ## What follows can be made more flexible later...
                 FC2 = Fave * funFC(x[[2]]$fbar) / FbarAve,
                 FC3 = Fave * funFC(x[[3]]$fbar) / FbarAve,
                 FC4 = Fave * funFC(x[[4]]$fbar) / FbarAve)

    colnames(res) <- c("Age", fitYears, tail(fcYears, -1))

    attr(res, "estType") <- estType

    return(res)
}

## Graphic functions:

label_parse <- function(breaks) {
   parse(text = breaks)
}

##
sdState <- function(fit, y = max(fit$data$years) - 1:0)
{
    idx <- names(fit$sdrep$value) == "logR"
    sdLogR <- fit$sdrep$sd[idx][fit$data$years%in%y]
    idx <- names(fit$sdrep$value) == "logssb"
    sdLogSSB <- fit$sdrep$sd[idx][fit$data$years %in% y]
    idx <- names(fit$sdrep$value) == "logfbar"
    sdLogF <- fit$sdrep$sd[idx][fit$data$years%in%y]
    ret <- cbind(sdLogR, sdLogSSB, sdLogF)
    rownames(ret) <- y
    colnames(ret) <- c("sd(log(R))", "sd(log(SSB))", "sd(log(Fbar))")
    return(ret)
}


### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 100
### End:
