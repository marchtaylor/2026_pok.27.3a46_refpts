# utilities.R
# Code adapted from code/functions in:
# https://github.com/ices-taf/wk_WKNSMSE_pok.27.3a46/blob/master/a4a_mse_WKNSMSE_funs.R
# General functions


# Function to collect objects from current environment that match function arguments;
# will replace default arguments, otherwise use any default arguments
get.args.list <- function(fun, preargs=NULL){
  mget(formalArgs(fun),
       envir = parent.frame(),
       ifnotfound=formals(fun), 
       inherits = TRUE) |>
    (\(x){
      if(!is.null(preargs)){modifyList(x=x,preargs,keep.null=TRUE)}else{x}
    })()
}


# Convert FLStock object to list of FLQuant objects(?) ####
to.flq <- function(stk,quant){
  require(FLCore)
  flq <- quant |> do.call(args=list(object=stk))
  flq
}


# Create FLStock object from SAM fit object ####
as.FLStock.sam <- function(fit, unit.w = "kg", name = "", desc = "", predicted = FALSE){
  # list2env(list(unit.w = "kg", name = "", desc = "", predicted = FALSE), globalenv())
  #require(FLCore)
  toFLQ <- function(x, unit = "NA"){
    FLCore::FLQuant(t(x),dimnames=list(age=as.numeric(colnames(x)),
                               year=as.numeric(rownames(x)),
                               unit="unique",
                               season="all",
                               area="unique",
                               iter=1),
            units = unit)
  }
  na2zero <- function(x){
    x[is.na(x)] <- 0
    x
  }
  resize1 <- function(x, ages, years, replicate = FALSE){
    x1 <- matrix(NA,length(years),length(ages))
    colnames(x1) <- ages
    rownames(x1) <- years
    
    # Extra check on if names are missing on x dims
    if(is.null(rownames(x))){
      # Will spit error if lengths of x dims and ages/years don't match
      rownames(x) <- years
    }
    if(is.null(colnames(x))){
      colnames(x) <- ages
    }
    
    x1[match(rownames(x),rownames(x1)),match(colnames(x),colnames(x1))] <- x
    # FIX THIS BELOW!!! It is not giveing right answers
    if(replicate){
      ## Years
      ii <- which(!match(rownames(x1),rownames(x), FALSE))
      if(length(ii) > 0){
        isAfter <- as.numeric(rownames(x1)[ii]) > max(as.numeric(rownames(x)))
        x1[ii,] <- x[c(1,nrow(x))[as.numeric(isAfter)+1],]
      }
      ## Ages
      ii <- which(!match(colnames(x1),colnames(x), FALSE))
      if(length(ii) > 0){
        isAfter <- as.numeric(colnames(x1)[ii]) > max(as.numeric(colnames(x)))
        x1[,ii] <- x[,c(1,ncol(x))[as.numeric(isAfter)+1]]
      }
    }
    x1
  }
  
  # Check if fit has estimated biopars (MO, NM, SW, CW)
  fitnms <- names(fit[["sdrep"]][["value"]])
  
  biopar2dat <- function(fit,datnm = "stockMeanWeight", parnm = "logSW"){
    rnam <- rownames(fit[["data"]][[datnm]])
    cnam <- colnames(fit[["data"]][[datnm]])
    fit[["data"]][[datnm]][] <- 
      fit[["sdrep"]][["value"]] |> 
      _[fitnms==parnm] |> 
      exp() |> 
      matrix(ncol=ncol(fit[["data"]][[datnm]])) |>
      _[1:nrow(fit[["data"]][[datnm]]),]
    rownames(fit[["data"]][[datnm]]) <- rnam
    colnames(fit[["data"]][[datnm]]) <- cnam
    return(fit[["data"]][[datnm]])
  }
  
  if(any("logSW" %in% fitnms)) fit[["data"]][["stockMeanWeight"]] <- biopar2dat(fit,"stockMeanWeight", "logSW")
  if(any("logitMO" %in% fitnms)) fit[["data"]][["propMat"]] <- biopar2dat(fit,"propMat", "logitMO") |> log() |> boot::inv.logit()
  if(any("logNM" %in% fitnms)) fit[["data"]][["natMor"]] <- biopar2dat(fit,"natMor","logNM")
  if(any("logCW" %in% fitnms)){
    fit[["data"]][["catchMeanWeight"]][,,1] <- biopar2dat(fit,"catchMeanWeight", "logCW")
    fit[["data"]][["disMeanWeight"]] <- fit[["data"]][["landMeanWeight"]] <- fit[["data"]][["catchMeanWeight"]]
  }else if(any("logSW" %in% fitnms)){
    fit[["data"]][["catchMeanWeight"]][,,1] <- biopar2dat(fit,"stockMeanWeight", "logSW")[1:nrow(fit[["data"]][["catchMeanWeight"]]),]
    fit[["data"]][["disMeanWeight"]] <- fit[["data"]][["landMeanWeight"]] <- fit[["data"]][["catchMeanWeight"]]
  }
  
  # Make sure dimnames are specified
  if(any(sapply(dimnames(fit$data$propMat),is.null))) rownames(fit$data$propMat) <- fit$data$years
  if(is.null(rownames(fit$data$stockMeanWeight))) rownames(fit$data$stockMeanWeight) <- fit$data$years
  if(is.null(rownames(fit$data$catchMeanWeight))) rownames(fit$data$catchMeanWeight) <- fit$data$years
  if(is.null(rownames(fit$data$disMeanWeight))) rownames(fit$data$disMeanWeight) <- fit$data$years
  if(is.null(rownames(fit$data$landMeanWeight))) rownames(fit$data$landMeanWeight) <- fit$data$years
  if(is.null(rownames(fit$data$natMor))) rownames(fit$data$natMor) <- fit$data$years
  
  ## catch.n
  stopifnot(dim(fit[["data"]][["catchMeanWeight"]])[3]==1) # Only works when catches are for one fishing fleet (not multiple)
  ages <- fit$conf$minAge:fit$conf$maxAge
  years <- fit$data$years
  CN <- toFLQ(resize1(na2zero(getFleet(fit,1, predicted)),ages,years))
  LF <- toFLQ(resize1(fit$data$landFrac[,,1],ages,years,TRUE))
  ## mat
  MA <- toFLQ(resize1(fit$data$propMat,ages,years,TRUE))
  ## stock.wt
  SW <- toFLQ(resize1(fit$data$stockMeanWeight,ages,years,TRUE), unit.w)
  ## stock.n
  SN <- toFLQ(resize1(ntable(fit),ages,years))
  ## catch.wt
  CW <- toFLQ(resize1(fit$data$catchMeanWeight[,,1],ages,years,TRUE), unit.w)
  ## discards.wt
  DW <- toFLQ(resize1(fit$data$disMeanWeight[,,1],ages,years,TRUE), unit.w)
  ## discards.n
  DN <- CN * (1 - LF)
  ## landings.wt
  LW <- toFLQ(resize1(fit$data$landMeanWeight[,,1],ages,years,TRUE), unit.w)
  ## landings.n
  LN <- CN * LF
  ## m
  M <- toFLQ(resize1(fit$data$natMor,ages,years,TRUE))
  ## harvest.spwn
  PF <- toFLQ(resize1(fit$data$propF,ages,years,TRUE))
  ## m.spwn
  PM <- toFLQ(resize1(fit$data$propM,ages,years,TRUE))
  ## harvest
  F <- toFLQ(resize1(faytable(fit),ages,years,FALSE), "f")
  FLCore::FLStock(catch = FLCore::apply(CN * CW,2,sum), # Total catch weight (‘FLQuant’)
          catch.n = CN, # Catch numbers (‘FLQuant’)
          catch.wt = CW, # Mean catch weights (‘FLQuant’)
          discards = FLCore::apply(DW * DN,2,sum), # Total discards weight (‘FLQuant’)
          discards.n = DN, # Discard numbers (‘FLQuant’)
          discards.wt = DW, # Mean discard weights (‘FLQuant’)
          landings = FLCore::apply(LW * LN,2,sum), # Total landings weight (‘FLQuant’)
          landings.n = LN, # Landing numbers (‘FLQuant’)
          landings.wt = LW, # Landing weights (‘FLQuant’)
          stock = FLCore::apply(SW * SN,2,sum), # Total stock weight (‘FLQuant’)
          stock.n = SN, # Stock numbers (‘FLQuant’)
          stock.wt = SW, # Mean stock weights (‘FLQuant’)
          m = M, # Natural mortality (‘FLQuant’)
          mat = MA, # Proportion mature (‘FLQuant’)
          harvest = F, # Harvest rate or fishing mortality. The units of this slot should be set to 'hr' or 'f' accordingly (‘FLQuant’)
          harvest.spwn = PF, # Proportion of harvest/fishing mortality before spawning (‘FLQuant’)
          m.spwn = PM, # Proportion of natural mortality before spawning (‘FLQuant’)
          name = name, # Name of the stock (‘character’)
          desc = desc, # Description of the stock (‘character’)
          range = c(min = fit$conf$minAge,
                    max = fit$conf$maxAge,
                    plusgroup = ifelse(fit$conf$maxAgePlusGroup[1],fit$conf$maxAge,NA),
                    minyear = min(fit$data$years),
                    maxyear = max(fit$data$years),
                    minfbar = fit$conf$fbarRange[1],
                    maxfbar = fit$conf$fbarRange[2]
          ) # Named numeric vector containing the quant and year ranges, the plusgroup and the quant range that the average fishing mortality should be calculated over (‘numeric’)
  )
}


# Calculates weighted Fbar from FLStock object ####
wfbar <- function (object, ...) 
{
  mina = consis_args$samfit$conf$fbarRange[1]
  maxa = consis_args$samfit$conf$fbarRange[2]
  # mina = range(object, "minfbar") 
  # maxa = range(object, "maxfbar")
  if (units(harvest(object)) == "f" || units(harvest(object)) == 
      "hr") {
    fval <- quantSums((harvest(object) * stock.n(object))[as.character(seq(mina,maxa)),]) /
      quantSums(stock.n(object)[as.character(seq(mina,maxa)),])
    return(fval)
  }
  else {
    stop("Correct units (f or hr) not specified in the harvest slot")
  }
}


# Simulate from truncated normal distribution ####
rtnorm <- function(n,mean = 0, sd = 1, a = -6, b = 6, ...){
  # Generate random uniform
  qsim = runif(n, min=0, max=1)
  
  # Put inverse distribution
  qinv = (pnorm(a,mean,sd,...) + qsim*(pnorm(b,mean,sd,...) - pnorm(a,mean,sd,...)))
  
  return(qnorm(qinv, mean, sd,...))
}


# Moving block bootstrap ####
mov.blk.boot <- function(ts, n=NULL, blk){
  # blk specification change type of resampling done:
  # blk=1         Simple resample with replacement
  # blk={>1}      Simple block bootstrap with replacement (and overlap)
  # blk={vector}  Moving block bootstrap, with set of potential block sizes specified in blk
  
  n0 <- length(ts) # Original time series length (does NOT need to match n)
  if(is.null(n)) n <- n0
  
  if(length(blk)>1){# blk is vector possible lengths from which to sample
    # Randomly sample lengths from blk
    # Vector of block lengths
    blki <- sample(blk,1)
    while(sum(blki)<n){blki <- c(blki,sample(blk,1))}
    #k <- length(blki) # Total number of blocks
    rblks <- sapply(blki,function(j) sample.int((n0-j+1),1,replace=TRUE)) # Starting index of each block of varying length
    rblks_end <-  rblks + blki -1  # Ending index ....
    newts <- sapply(seq_along(rblks), function(i){ts[rblks[i]:rblks_end[i]]}) |> unlist() # Extract and concatenate blocks
  }else{
    j <- n0 - blk + 1 # Total number of overlapping blocks (starting indices) to be sampled from 
    k <- ceiling(n/blk) # Total number of blocks to draw
    rblks <- sample.int(j,k,replace=TRUE)
    rblks_end <-  rblks + blk -1
    newts <- sapply(seq_along(rblks), function(i){ts[rblks[i]:rblks_end[i]]}) |> c()
  }
  
  newts[1:n]
}


# Conditional moving block boostrap ####
# Resampling of sequential block is dependent on the preceding block
cond.mov.blk.boot <- function(ts, n=NULL, blk, cond_sets, use_lastobs=TRUE){
  # cond.sets should be long formatted matrix or data.frame with 3 columns identifying 
  # all valid sequential combinations between observations. 
  # In each row, column 1 is ts value in current time step, column 2 is possible 
  # ts value in following time step, and column 3 is a logical (TRUE/FALSE) whether
  # this sequence of values is acceptable
  
  # Original time series length (does NOT need to match n)
  n0 <- length(ts)
  if(is.null(n)) n <- n0
  
  # First sample  sequence of block sizes
  blkszs <- sample(blk, ceiling(n/min(blk)),replace=T)
  
  # Truncate sequence at cumsum nearest n (to the right)
  blkszs <- blkszs[1:(max(which(cumsum(blkszs)<n))+1)]
  
  # Find indices of transitions between blocks
  blkends <- cumsum(blkszs)
  
  # Loop through blocks and sample starting indices of each block
  rblks <- rep(NA,length(blkszs))
  for(i in seq_along(blkszs)){
    # ri is starting index of each block of varying length
    j <- blkszs[i]
    if(i==1){
      if(use_lastobs){
        imin1 <- n0
      }else{
        imin1 <- sample.int(n0,1)
      }
    }else{
      imin1 <- rblks[i-1] + blkszs[i-1] -1
    }
    subts <- cond_sets[cond_sets[,1]==ts[imin1],] |> 
      (\(x) x[,2][as.logical(x[,3])])()
    subts <- which(ts%in%subts)
    rblks[i] <- sample(subts[subts<=(n0-j+1)],1)
  }
  rblks_end <-  rblks + blkszs -1
  newts <- sapply(seq_along(rblks), function(i){ts[rblks[i]:rblks_end[i]]}) |> unlist() # Extract and concatenate blocks
  return(newts[1:n])
}


# Get historical process residuals (on logN) ####
# Extracted from stockassessment:::procres()
get.procres.sdrep <- function(fit,...){
  pp <- fit$pl
  attr(pp,"what") <- NULL
  pp$missing <- NULL
  #fit.co <- stockassessment:::refit(fit,startingValues=pp, run=FALSE)
  fit.co <- stockassessment:::sam.fit(fit$data,fit$conf,pp, run=FALSE, map=fit$obj$env$map)
  fit.co$obj$env$data$resFlag <- 1
  fit.co$obj$retape()
  fit$sdrep <- TMB::sdreport(fit.co$obj,fit$opt$par,...)
  return(fit)
}


# Simulate SAM parameter vectors from estimation uncertainty ####
# Returns list alternative (random) values of SAM parameters based on SAM 
# uncertainty (variance-covariance with multivariate normal)
# TODO: Argument to exclude (or reduce size?) 'obj' from list of samfit objects because it is quite large
SAM.uncertainty <- function(fit, 
                            n = 10, 
                            print_screen = FALSE, 
                            use_fast_mvrnorm=TRUE) {
  if (!is(fit, "sam")) stop("fit has to be class \"sam\"")
  
  ##Resample estimated values to get N, F and q 
  fixed_nms <- unique(names(fit$opt$par))
  fixed_plsd <- unlist(fit$plsd[names(fit$plsd)%in%fixed_nms])
  chk <- which(fixed_plsd>6 & !is.na(fixed_plsd))
  if(length(chk)>0){
    invisible(lapply(chk,function(i){
      p1 <- names(fixed_plsd)[i]
      if(!grepl("Sd",p1)){
        stop("Very large SD (>6) estimated for ",p1,"; fix parameter and re-fit.")
      }else{
        warning("Very large SD (>6) estimated for ",p1," (",round(fixed_plsd[i],3),"); will simulate parameter vectors, but suggest fixing parameter and re-fitting.")
      }
    }))
  }
  
  ### calculate standard deviations of model parameters
  if("sdrep"%in%names(fit)){
    if("jointPrecision"%in%names(fit$sdrep)) sds <- fit$sdrep
  }
  
  if(!exists("sds")){
    . <- capture.output(sds <- TMB::sdreport(obj = fit$obj, 
                                             par.fixed = fit$opt$par,
                                             getJointPrecision = TRUE))
  }
  
  ### extract values for parameters
  est <- c(sds$par.fixed, sds$par.random)
  ### get covariance matrix of all model parameters
  cov <- Matrix::solve(sds$jointPrecision)
  parnms <- names(est)
  
  ### create random values based on estimation and covariance
  # set.seed(62390233)
  if(use_fast_mvrnorm){
    sim.pars <- mvnfast::rmvn(n, est, cov) # 24.10.2024:  39 s for 10 samples and est vector size of 5756
    #sim.states <- mvnfast::rmvn(n, sds$value, sds$cov) 
  }else{
    sim.pars <- mvrnorm(n, est, cov) # 24.10.2024:  453 s for 10 samples and est vector size of 5756
    #sim.states <- mvrnorm(n=1, sds$value, sds$cov) 
  }
  
  ### combine SAM estimate and random samples
  lastpar0 <- fit[["obj"]][["env"]][["last.par"]]
  data2 <- fit$obj$env$data
  #data2$reportingLevel <- as.integer(TRUE)
  data2$resFlag <- 1
  args <- list(data       = data2, 
               #parameters = fit$obj$env$parameters, 
               parameters = fit$pl, 
               map        = fit[["obj"]][["env"]][["map"]],
               type       = "ADFun",
               #random     = fit[["obj"]][["env"]][[".random"]],
               ADreport   = TRUE, 
               DLL        = fit$obj$env$DLL, 
               silent     = fit$obj$env$silent)
  
  obj2 <- do.call(TMB::MakeADFun,args)
  
  #unc <- sim.pars
  if(n==1) sim.pars <- matrix(sim.pars,nrow=1)
  #print("Iter ")
  unc <- lapply(1:nrow(sim.pars), function(x){
    
    #attributes(pl) <- attributes(fit$obj$env$parameters)
    #attr(pl, "check.passed") <- TRUE
    last.par <- lastpar0
    last.par[] <- sim.pars[x,]
    obj2[["env"]][["last.par"]] <- last.par
    pl <- obj2[["env"]][["parList"]](last.par)
    # pn <- try(obj2$report(last.par), silent = TRUE)[["pn"]]
    # resN <- pl[["logN"]] - pn
    
    phi <- try(obj2$fn(last.par), silent = TRUE)
    #phi["resN"] <- c(t(resN[,-1]))
    sdrep <- list(value=phi, sd=phi*NA)
    
    pseudofit <- list(data=fit$data,
                      conf=fit$conf,
                      pl=pl,
                      last.par=last.par,
                      obj=obj2,
                      rep=obj2$report(),
                      sdrep=sdrep)
    class(pseudofit) <- "sam"
    pseudofit
  })
  
  ## Below extracted from stockassessment::sdreport()
  # par <- obj$env$last.par.best
  # phi <- try(obj2$fn(par), silent = TRUE)
  # ans <- list(value = phi, sd = sd, cov = cov, par.fixed = par.fixed, 
  #             cov.fixed = Vtheta, pdHess = pdHess, gradient.fixed = gradient.fixed)
  
  return(unc)
}

