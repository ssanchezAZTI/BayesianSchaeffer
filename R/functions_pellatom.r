# Complementary functions for Pella-Tomlinson (and Schaeffer as specific case) model in R
# BayesianSchaeffer/R/functions_pellatom.R 

# Copyright 2015 AZTI Team. Distributed under the GPL 2 or later 
# Maintainer: Leire Ibaibarriaga (libaibarriaga@azti.es) 
# $Id: functions_pellatom.R 14/9/2015 libaibarriaga $ 


####################################################################
####################################################################
#
# Functions:
#
# - calc.biomass :
#   function to calculate P (biomass/k) with or without process error 
# - calc.nuisance: 
#   function to derive q and var.logI max lik estimates (assuming no process error) 
#   only for 1 index
# - ff:
#   log likelihood function to be maximised
# - jags2coda:
#   generic function to transform the jags object into coda 
# - pellatom.inits
#   function to randomly draw initial values
# - plot.chain: 
#   function to get the standard figures from coda by parameter
# - proj.biomass:
#   function for short term projections (only one iteration)
#
####################################################################
####################################################################

# function to calculate P (biomass/k) with or without process error 

calc.biomass <- function(r,k,p,P0,var.logP,C){
  n <- length(C)
  P <- rep(0, n)
  prerror <- rnorm(n, 0, sqrt(var.logP)) # process error in log scale
  P[1] <- P0 * exp(prerror[1])
  for (i in 2:n){
    P[i] <- ( P[i-1]+r/p*P[i-1]*(1-(P[i-1]^p))-C[i-1]/k ) * exp(prerror[i])
  }
  ok <- all(P>0)
  out <- c(P, ok)
  names(out) <- c(paste("P[", 1:n, "]", sep = ""), "ok")
  return(out)
}

####################################################################
####################################################################

# function to derive q and var.logI max lik estimates (assuming no process error)
# only for 1 index

calc.nuisance <- function(r, k, p=NULL, P0=NULL, Index, C){
  var.logP <- 0 # assuming var.logP=0, i.e. no process error
  if(is.null(P0)){
    P0 <- 1
  }
  if(is.null(p)){
    p <- 1
  }
  n <- length(C)
  aux <- calc.biomass(r=r,k=k,p=p,P0=P0,var.logP=var.logP,C=C) 
  ok <- aux["ok"]
  if(ok){
    P <- aux[paste("P[", 1:n, "]", sep = "")]  # same length as the C vector
    q <- exp(1/n*sum(log(Index)-log(k)-log(P)))
    sig <- sqrt(1/n*sum( (log(Index)-log(q)-log(k)-log(P))^2 ))
  }else{
    q <- NA
    sig <- NA
  }
  return(list(q=q, var.logI=sig^2))
}

####################################################################
####################################################################

# function to compute the log-likelihood for the case with only observation error
# it assumes the same number of years in the indices and C
# it allows Index to be a vector for 1 index of a list of vectors for 2 or more indices
# note that q and var.logI will be derived analitically

# log likelihood function to be maximised

ff <- function(vec, p=NULL, P0=NULL, Index, C){
  r <- vec[1]
  k <- vec[2]
  if(is.null(P0)){
    P0 <- 1
  }
  if(is.null(p)){
    p <- 1
  }
  aux <- calc.biomass(r=r,k=k,p=p,P0=P0,var.logP=0,C=C) # assuming var.logP=0, i.e. no process error
  ok <- aux["ok"]
  n <- length(C)
  ll <- 0
  if (ok){ 
    P <- aux[paste("P[", 1:n, "]", sep = "")]  # same length as the C vector
    if (is.list(Index)){
      for (j in 1:length(Index)){  # for each index
        qsig <- calc.nuisance(r=r, k=k, p=p, P0=P0, Index=Index[[j]], C=C)
        q <- qsig$q
        var.logI <- qsig$var.logI
        for (i in 1:n){
          ll <- ll + dnorm(log(Index[[j]][i]), mean=log(q)+log(k)+log(P[i]), sd=sqrt(var.logI), log=TRUE)    
        }
      }
    }else{
      q <- exp(1/n*sum(log(Index)-log(k)-log(P)))
      sig <- sqrt(1/n*sum( (log(Index)-log(q)-log(k)-log(P))^2 ))
      for (i in 1:n){
        ll <- ll + dnorm(log(Index[i]), mean=log(q)+log(k)+log(P[i]), sd=sig, log=TRUE)    
      }      
    }  
  }else{
    ll <- -1.0e+6 # we give it very low values at the incorrect area not to affect the maximization
  }  
  return(ll)
}

####################################################################
####################################################################

# function to randomly draw initial values

pellatom.inits <- function(C,
                           mu.logr=-1.4, psi.logr=3.8,
                           mu.logk=5, psi.logk=3.7,
                           mu.logq=0, psi.logq=1,
                           mu.logp=0, psi.logp=1,
                           a.P0=0, b.P0=1,                     
                           a.psi.logI=3.8, b.psi.logI=0.01,
                           a.psi.logP=1.7, b.psi.logP=0.009,
                           fx.logp=T, fx.P0=T,
                           n.chains = 1, maxit = 100){
  out <- vector("list", length = n.chains)
  for (i in 1:n.chains) {
    ok <- F
    iter <- 0
    cont <- ((!ok) & (iter < maxit))
    while (cont) {
      logr <- rnorm(1, mu.logr, sqrt(1/psi.logr))
      logk <- rnorm(1, mu.logk, sqrt(1/psi.logk))
      logq <- rnorm(1, mu.logq, sqrt(1/psi.logq))
      logp <- ifelse(fx.logp, mu.logp, rnorm(1, mu.logp, sqrt(1/psi.logp))) 
      P0 <- ifelse(fx.P0, (a.P0+b.P0)/2, runif(1, a.P0, b.P0))
      psi.logI <- rgamma(1, a.psi.logI, b.psi.logI)
      psi.logP <- rgamma(1, a.psi.logP, b.psi.logP)
      
      tmp <- calc.biomass(r=exp(logr),k=exp(logk),p=exp(logp),P0=P0,var.logP=1/psi.logP,C=C)
      ok <- tmp["ok"]
      P <- tmp[paste("P[", 1:n, "]", sep = "")]
      iter <- iter + 1
      cont <- ((!ok) & (iter < maxit))
    }
    if (!ok) {
      print(paste("The maximum number of iterations has been reached for chain", 
                  i))
    }
    logP <- as.numeric(log(P))
    if(fx.logp==T & fx.P0==T){
      out[[i]] <- list(logr=logr, logk=logk, logq=logq, psi.logI=psi.logI, psi.logP=psi.logP, logP=logP)        
    }else if(fx.logp==T & fx.P0==F){
      out[[i]] <- list(logr=logr, logk=logk, logq=logq, P0=P0, psi.logI=psi.logI, psi.logP=psi.logP, logP=logP)        
    }else if(fx.logp==F & fx.P0==T){
      out[[i]] <- list(logr=logr, logk=logk, logq=logq, logp=logp, psi.logI=psi.logI, psi.logP=psi.logP, logP=logP)        
    }else if(fx.logp==F & fx.P0==F){
      out[[i]] <- list(logr=logr, logk=logk, logq=logq, logp=logp, P0=P0, psi.logI=psi.logI, psi.logP=psi.logP, logP=logP)        
    }  
    
  }
  return(out)     
}

################################################################################
################################################################################

# generic function to transform the jags object into coda 
# it is taken from the function coda.samples(), which is a wrapper of the functions jags.samples() in the package rjags

jags2coda <- function(out){
  nchain <- as.numeric(dim(out[[1]])["chain"])
  if(is.na(nchain)){
    nchain <- 1
  }
  ans <- vector("list", nchain)
  for (ch in 1:nchain) {
    ans.ch <- vector("list", length(out))
    vnames.ch <- NULL
    for (i in seq(along = out)) {
      varname <- names(out)[[i]]
      d <- dim(out[[i]])
      if (length(d) < 3) {
        stop("Invalid dimensions for sampled output")
      }
      vardim <- d[1:(length(d) - 2)]
      nvar <- prod(vardim)
      niter <- d[length(d) - 1]
      nchain <- d[length(d)]
      values <- as.vector(out[[i]])
      var.i <- matrix(NA, nrow = niter, ncol = nvar)
      for (j in 1:nvar) {
        var.i[, j] <- values[j + (0:(niter - 1)) * nvar + 
                               (ch - 1) * niter * nvar]
      }
      vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
      ans.ch[[i]] <- var.i
    }
    ans.ch <- do.call("cbind", ans.ch)
    colnames(ans.ch) <- vnames.ch
    ans[[ch]] <- mcmc(ans.ch)
  }
  mcmc.list(ans)
}

################################################################################
################################################################################

# function to get the standard figures from coda by parameter

plot.chain <- function(mcmc, probs=c(0.05, 0.5, 0.95), lag.max=100, param.names=NULL, filename="plot.pdf"){
  
  # mcmc is an mcmc type object (a unique chain; sublist of mcmc.list after applying the burn-in, the thinning, etc)
  # probs is a vector of length 3 with the probs to plot the low, median and upper quantiles
  
  if (length(probs)!=3)
    stop("The length of vector probs must be 3")
  
  if (is.null(param.names)){
    param.names <- varnames(mcmc)
  }
  
  pdf(filename)
  for (param in param.names){
    par(mfrow=c(2,2))
    chain <- mcmc[,param]
    quant <- quantile(chain, probs=probs)
    plot(as.vector(chain), type="l", xlab="", ylab="", main="Trace")
    abline(h=quant, lty=c(2,1,2), col=2)
    cumuplot(chain, probs=probs, auto.layout=F, xlab="", ylab="", main="Cumulative quantile")
    abline(h=quant, lty=c(2,1,2), col=2)
    try(acf(chain, lag.max=lag.max, xlab="", ylab="", main="Auto-correlation"), silent=T)
    hist(chain, freq=F, xlab="", ylab="", main="Density") 
    lines(density(chain, from=0))
    abline(v=quant, lty=c(2,1,2), col=2)
    title(param, outer=T, line=-2)
  }
  par(mfrow=c(1,1))
  crosscorr.plot(mcmc[, param.names])
  dev.off()
}

################################################################################
################################################################################

# function for short term projections (only one iteration)
# I don't use the function calc.biomass because it adds process error in the starting point

proj.biomass <- function(r, k, p=1, Pini, var.logP, C=NULL, F=NULL, type="C", nproj=2){
  if (type=="C"){
    if(is.null(C)){
      C <- rep(0, nproj)
    }else if(length(C)==1 & nproj>1){
      C <- rep(C, nproj)
    }
    F <- rep(0, nproj)
  } 
  if (type=="F"){
    if(is.null(F)){
      F <- rep(0, nproj)      
    }else if(length(F)==1 & nproj>1){
      F <- rep(F, nproj)
    }
    C <- rep(0, nproj)
  } 
  proj <- vector("numeric", length=nproj)
  prerror <- rnorm(nproj, 0, sqrt(var.logP)) # process error in log scale
  for (j in 1:nproj){
    Plast <- ifelse(j==1, Pini, proj[j-1])
    if (type=="F"){
      C[j] <- F[j]*Plast*k
    }else if (type=="C"){
      F[j] <- C[j]/(Plast*k)
    }
    proj[j] <- ( Plast+r/p*Plast*(1-(Plast^p))-C[j]/k ) * exp(prerror[j])
  }
  Bproj <- proj * as.numeric(k)
  out <- list(Bproj=Bproj, Cproj=C, Fproj=F) # be careful because B and (C or F) are not of the same year, lag of 1 year
  return(out)
}


