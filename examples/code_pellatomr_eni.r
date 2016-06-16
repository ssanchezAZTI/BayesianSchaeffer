# Pella-Tomlinson (and Schaeffer as specific case) model in R
# BayesianSchaeffer/examples/code_pellatomr_eni.R 

# Copyright 2015 AZTI Team. Distributed under the GPL 2 or later 
# Maintainer: Leire Ibaibarriaga (libaibarriaga@azti.es) 
# $Id: code_pellatomr_eni.R 26/11/2015 libaibarriaga $ 

####################################################################
####################################################################
#
# Orange spotted grouper (Epinephelus coioides)      - ENI - AZTI
#
####################################################################
####################################################################

# species

sp <- "ENI"

####################################################################
####################################################################

# set working directory

wd <- paste(".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/",tolower(sp),"/surplus/",sep="")
  
setwd(wd)

####################################################################
####################################################################

# load data file

sa.dat <- read.table(paste(".../IM12KFUPM/WP5.2_PopDyn&StockAss/data_analysis/data_",sp,".txt",sep=""), header=T)

# figures of time series of catches, ndays and cpue

jpeg("series_catch_sa_recofi.jpg", width=700, quality=100)
matplot(sa.dat$year, cbind(sa.dat$int.catch, sa.dat$sa.catch), type="o", col=1, lty=1, pch=c(1, 19), xlab="Year", ylab="Catch (t)")
dev.off()

jpeg("series_catch.jpg", width=700, quality=100)
plot(sa.dat$year, sa.dat$sa.catch, type="o", pch=19, xlab="Year", ylab="Catch (t)")
dev.off()

jpeg("series_ndays.jpg", width=700, quality=100)
plot(sa.dat$year, sa.dat$days, type="o", pch=19, xlab="Year", ylab="Nb days")
dev.off()

jpeg("series_cpue.jpg", width=700, quality=100)
plot(sa.dat$year, sa.dat$cpue, type="o", pch=19, xlab="Year", ylab="CPUE (t/day)")
dev.off()

# subset without NA's

sa.dat <- subset(sa.dat, !is.na(sa.dat$cpue) & !is.na(sa.dat$sa.catch))

# number of years

n <- dim(sa.dat)[1]

####################################################################
####################################################################

# load library for new color palette

library(fields)

# load library for JAGS

library(rjags)

# load library to add text labels to kobe plots

library(calibrate)

# load library for function kde2d for bi-variate kernel density

library(MASS)

####################################################################
####################################################################

# load complementary functions

source(".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/code_surplus/functions_pellatom.r")

####################################################################
####################################################################

# grid search for max.lik.

# combination of r and k that fulfil the restrictions
# and to see the value of the log-likelihood
# p=1 (Schaeffer model) and P0<1 (initial population lower than carrying capacity)

#------------------------------------------------------------------

rvec <- seq(0, 5, by=0.05)
kvec <- seq(0, 50000, by=500)
P0vec <- seq(0.1, 1, by=0.1)

print(Sys.time())
mat0 <- array(NA, dim=c(length(rvec), length(kvec), length(P0vec)))
mat1 <- array(NA, dim=c(length(rvec), length(kvec), length(P0vec)))
out <- NULL
for (i in 1:length(rvec)){
  for (j in 1:length(kvec)){
    for (h in 1:length(P0vec)){
      mat0[i, j, h] <- as.numeric(calc.biomass(rvec[i],kvec[j],1,P0vec[h],var.logP=0,C=sa.dat$sa.catch)["ok"])
      mat1[i, j, h] <- ff(c(rvec[i], kvec[j]), P0=P0vec[h], Index=sa.dat$cpue, C=sa.dat$sa.catch)
      qsig <- calc.nuisance(rvec[i],kvec[j],1,P0vec[h],Index=sa.dat$cpue,C=sa.dat$sa.catch)
      out <- rbind(out, c(rvec[i],kvec[j],P0vec[h],mat0[i,j,h],mat1[i,j,h], qsig$q, qsig$var.logI))
    }
  }  
}
print(Sys.time())
out <- as.data.frame(out)
names(out) <- c("r","k","P0","ok","ll","q","var.logI")

save(out, file="res_loglik.RData")

# check magnitude of loglik before plotting

summary(out$ll[out$ok>0])

# figure of loglik surface 

jpeg("surface_loglik.jpg", width=700, quality=100)
par(mfrow=c(2,5), mar=c(3,3,2,0)+0.1)
for (h in 1:length(P0vec)){
#  image(rvec, kvec, mat0[,,h], col=grey(c(1,0)), xlab="r", ylab="k", main="Valid values") # valid in black
  image(rvec, kvec, mat1[,,h], zlim=c(-50,-20), breaks=seq(-50,-20,by=1),col=tim.colors(30), xlab="",ylab="", cex.axis=1.5, cex.main=1.5, main=paste("P0 =",P0vec[h]), xlim=range(rvec),ylim=range(kvec)) # darker colours indicate larger values
}
dev.off()

# first 20 maximum ll values found by grid search

out[order(out$ll, decreasing=T)[1:20],]

#------------------------------------------------------------------

# Maximum likelihood. Use optim for maximization (control parameter changed to maximise instead of minimise)

maxlik <- optim(c(0.5,30000,0.7), function(v,Index,C){ff(v[1:2], P0=v[3],Index=Index, C=C)}, Index=sa.dat$cpue, C=sa.dat$sa.catch, method="L-BFGS-B", control=list(fnscale=-1, parscale=c(1,10000,1)), lower=c(0,0,0), upper=c(5,50000,1))
maxlik

# it gives $par. On the upper bound for k!!
# 5.741415e-01 1.000000e+05 3.716242e-01

aux <- calc.biomass(maxlik$par[1],maxlik$par[2],1,maxlik$par[3],var.logP=0, C=sa.dat$sa.catch)
qsig <- calc.nuisance(maxlik$par[1],maxlik$par[2],1,maxlik$par[3],Index=sa.dat$cpue, C=sa.dat$sa.catch)

# comparison observed and fitted

jpeg("maxlik_fit.jpg", width=700, quality=100)
par(mfrow=c(1,1))
matplot(sa.dat$year, cbind(sa.dat$cpue, aux[grep("P", names(aux))]*maxlik$par[2]*qsig$q), type="l", ylab="CPUE")
legend("topleft", c("observed","fitted"), col=1:2, lty=1:2)
dev.off()

####################################################################
####################################################################

# Bayesian fit

#------------------------------------------------------------------

# run identifier:

run.name <- "run1"

# first set of prior distributions

mu.logr <- log(0.3)   # based on r=2*Fmsey=2*0.87*M (from Zhou et al 2012) with M=0.19 from Grandcourt 2005
psi.logr <- 1         # it corresponds to a CV (in normal scale) of 131% (cv=1.31=sqrt(exp(1/1)-1)  )
mu.logk <-  log(2000) 
psi.logk <- 1
mu.logq <- log(1e-7)
psi.logq <- 1
mu.logQ <- mu.logq + mu.logk  # reparameterization. the mean is the sum of means and variance is the sum of variances
psi.logQ <- (psi.logq + psi.logk) /(psi.logq*psi.logk)
a.psi.logI <- 1          # 1/sqrt(a) is the CV of the gamma prior distribution
b.psi.logI <- 0.01 #0.01
a.psi.logP <- 1
b.psi.logP <- 0.01 #0.01 # 0.0005 # 0.0025
a.P0 <- 0
b.P0 <- 1
mu.logp <- log(1)
psi.logp <- 1

calc.biomass(r=exp(mu.logr),k=exp(mu.logk),p=1,P0=0.4,var.logP=b.psi.logP/a.psi.logP,C=sa.dat$sa.catch)
calc.nuisance(r=exp(mu.logr),k=exp(mu.logk),p=1,P0=0.4,Index=sa.dat$cpue, C=sa.dat$sa.catch)

# vector to compute quantiles

qvec <- c(0.05, 0.1, 0.5, 0.9, 0.95)

tab.priors <- NULL
tab.priors <- rbind(tab.priors, c(mu.logr, psi.logr, qlnorm(qvec, mu.logr, 1/sqrt(psi.logr))))
tab.priors <- rbind(tab.priors, c(mu.logk, psi.logk, qlnorm(qvec, mu.logk, 1/sqrt(psi.logk))))
tab.priors <- rbind(tab.priors, c(mu.logq, psi.logq, qlnorm(qvec, mu.logq, 1/sqrt(psi.logq))))
tab.priors <- rbind(tab.priors, c(a.P0, b.P0, qunif(qvec, a.P0, b.P0)))
tab.priors <- rbind(tab.priors, c(mu.logp, psi.logp, qlnorm(qvec, mu.logp, 1/sqrt(psi.logp))))
tab.priors <- rbind(tab.priors, c(a.psi.logI, b.psi.logI, qgamma(qvec, a.psi.logI, b.psi.logI)))
tab.priors <- rbind(tab.priors, c(a.psi.logP, b.psi.logP, qgamma(qvec, a.psi.logP, b.psi.logP)))

tab.priors <- cbind(parameter=c("r","k","q","P0","p","psi.logI","psi.logP"),as.data.frame(tab.priors))

# save the prior table outside

write.table(tab.priors, file=paste("table_priors_",run.name,".txt",sep=""), row.names=F, col.names=F)

# The gamma distribution has mean=a/b, var=a/b^2, cv=1/sqrt(a) 

qgamma(qvec, a.psi.logI, b.psi.logI)   

# this leads to CV of Index: 

sqrt(exp(1/qgamma(qvec, a.psi.logI, b.psi.logI))-1)

#psi.logP

qgamma(qvec, a.psi.logP, b.psi.logP)   

# this leads to CV of process equations: 

sqrt(exp(1/qgamma(qvec, a.psi.logP, b.psi.logP))-1)

#------------------------------------------------------------------

# especificar longitudes de MCMC. se pueden dar más abajo al aplicar las funciones directamente

n.chains <- 1
n.adapt <- 10000
n.burnin <- 100000
n.thin <- 100
n.iter <- 1000000

#------------------------------------------------------------------

# write data to a txt file

dat <- list(n=n, I=sa.dat$cpue, C=sa.dat$sa.catch,
            mu.logr=mu.logr, psi.logr=psi.logr,
            mu.logk=mu.logk, psi.logk=psi.logk,
            mu.logQ=mu.logQ, psi.logQ=psi.logQ,
            #mu.logq=mu.logq, psi.logq=psi.logq,
            a.P0=a.P0, b.P0=b.P0,                           
            a.psi.logI=a.psi.logI, b.psi.logI=b.psi.logI,
            a.psi.logP=a.psi.logP, b.psi.logP=b.psi.logP)

# default starting values in JAGS

#inits <- NULL

# generate initial values that fulfill the catch restrictions

inits <- pellatom.inits(C=sa.dat$sa.catch,
                        mu.logr=mu.logr, psi.logr=psi.logr,
                        mu.logk=mu.logk, psi.logk=psi.logk,
                        mu.logq=mu.logq, psi.logq=psi.logq,
                        mu.logp=mu.logp, psi.logp=psi.logp,
                        a.P0=a.P0, b.P0=b.P0,                     
                        a.psi.logI=a.psi.logI, b.psi.logI=b.psi.logI,
                        a.psi.logP=a.psi.logP, b.psi.logP=b.psi.logP,
                        fx.logp=T, fx.P0=F,
                        n.chains = n.chains, maxit = 100)

# modify to include logQ

for (i in 1:n.chains){
  inits[[i]]$logQ <- inits[[i]]$logq + inits[[i]]$logk  
}

# list of parameters to be saved

param <- c("logr", "logk", "logQ", "P0", "psi.logI", "psi.logP", "r", "k", "Q", "q", "P", "Fmsy", "Bmsy", "MSY", "dif")

# path for the JAGS model

model.path <- ".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/code_surplus/s1.jags"

ptm.ini <- Sys.time()
print(Sys.time())
m <- jags.model(model.path, dat=dat, inits=inits, n.chains=n.chains, n.adapt=n.adapt)
print(Sys.time())

# burn-in

print(Sys.time())
update(m, n.iter=n.burnin)
print(Sys.time())

# jags object with MCMC

run.jags <- jags.samples(m, variable.names=param, n.iter=n.iter, thin=n.thin)
ptm.end <- Sys.time()
ptm <- ptm.end - ptm.ini
print(ptm)

# Time difference of 2.630691 mins for 1 chain

# save results externally

save(run.jags, file=paste("res_bay_",run.name,".RData",sep=""))

#------------------------------------------------------------------

# check that the max function did not have any effect. look at min and max values. they should be 0

apply(run.jags$dif, 1, quantile, c(0,0.5,1))

#------------------------------------------------------------------

# change from the jags.samples format to coda format

run.mcmc <- jags2coda(run.jags)

# additional selection

#run.mcmc <- window(run.mcmc, start=500)

# full list of names of parameters to be estimated

param.names <- c("logr", "logk", "logQ", "P0",
                 "psi.logI", "psi.logP",
                 "r","k","Q","q",   
                 paste("P[",1:n,"]",sep=""),
                 "Fmsy", "Bmsy","MSY")

param.names <- param.names[param.names %in% varnames(run.mcmc)] # to filter in case some of the parameters are missing (e.g. fixed)

# re-order the mcmc object according to the param.names vector

run.mcmc <- run.mcmc[,param.names]

#summary
summary(run.mcmc)

# check results from 2 chains

# cbind(summary(run.mcmc[[2]])$quantiles[, c(1,3,5)], summary(run.mcmc[[1]])$quantiles[,c(1,3,5)])

# check standard plots made with coda library

pdf(paste("plot_coda_",run.name,".pdf",sep=""))
plot(run.mcmc)
dev.off()

# select only one chain

run.mcmc <- run.mcmc[[1]]

# vector of quantiles to use to get the credible intervals and the estimates

qvec <- c(0.05, 0.5, 0.95)

# summary statistics using coda

summary(run.mcmc)
write.table(summary(run.mcmc, quantiles=qvec)$quantiles, file=paste("table_summary_parameters_",run.name,".txt",sep=""), col.names=F)

# check plots

plot.chain(run.mcmc, filename=paste("plot_",run.name,".pdf",sep=""))

#------------------------------------------------------------------

# pdf with figures for comparison prior and posterior

pdf(paste("prior_posterior_",run.name,".pdf",sep=""))

# log(r)
xmin <- min(qnorm(c(0.025, 0.975), mu.logr, sqrt(1/psi.logr)), run.jags$logr)
xmax <- max(qnorm(c(0.025, 0.975), mu.logr, sqrt(1/psi.logr)), run.jags$logr)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$logr, bw=(xmax-xmin)/100)
yprior <- dnorm(x, mu.logr, sqrt(1/psi.logr))
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main="log(r)", yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lty=2, lwd=2)

# log(k)
xmin <- min(qnorm(c(0.025, 0.975), mu.logk, sqrt(1/psi.logk)), run.jags$logk)
xmax <- max(qnorm(c(0.025, 0.975), mu.logk, sqrt(1/psi.logk)), run.jags$logk)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$logk, bw=(xmax-xmin)/100)
yprior <- dnorm(x, mu.logk, sqrt(1/psi.logk))
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main="log(k)", yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lty=2, lwd=2)

# log(q)
xmin <- min(qnorm(c(0.025, 0.975), mu.logq, sqrt(1/psi.logq)), log(run.jags$q))
xmax <- max(qnorm(c(0.025, 0.975), mu.logq, sqrt(1/psi.logq)), log(run.jags$q))
x <- seq(xmin, xmax, length=10000)
y <- density(log(run.jags$q), bw=(xmax-xmin)/100)
yprior <- dnorm(x, mu.logq, sqrt(1/psi.logq))
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main="log(q)", yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lty=2, lwd=2)

# psi.logI

xmin <- min(qgamma(c(0.025, 0.975), a.psi.logI, b.psi.logI), run.jags$psi.logI)
xmax <- max(qgamma(c(0.025, 0.975), a.psi.logI, b.psi.logI), run.jags$psi.logI)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$psi.logI, bw=(xmax-xmin)/100)
yprior <- dgamma(x, a.psi.logI, b.psi.logI)
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
ymin <- 0
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=expression(psi[log(I)]), yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lwd=2, lty=2)

# psi.logP

xmin <- min(qgamma(c(0.025, 0.975), a.psi.logP, b.psi.logP), run.jags$psi.logP)
xmax <- max(qgamma(c(0.025, 0.975), a.psi.logP, b.psi.logP), run.jags$psi.logP)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$psi.logP, bw=(xmax-xmin)/100)
yprior <- dgamma(x, a.psi.logP, b.psi.logP)
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
ymin <- 0
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=expression(psi[log(P)]), yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lwd=2, lty=2)

# P0 prior and posterior

xmin <- min(qunif(c(0.025, 0.975), a.P0, b.P0), run.jags$P0)
xmax <- max(qunif(c(0.025, 0.975), a.P0, b.P0), run.jags$P0)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$P0, bw=(xmax-xmin)/50)
yprior <- dunif(x, a.P0, b.P0)
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
ymin <- 0
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=expression(P[0]), yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lwd=2, lty=2)

dev.off()

#------------------------------------------------------------------

# calculate B, B/Bmsy, F, F/Fmsy, Indexfit time series

biomass <- array(NA, dim(run.jags$P))
b.bmsy <- array(NA, dim(run.jags$P))
indexfit <- array(NA, dim(run.jags$P))
f <- array(NA, dim(run.jags$P))
f.fmsy <- array(NA, dim(run.jags$P))
for (i in 1:dim(run.jags$P)[1]){
  biomass[i,,] <- run.jags$P[i,,] * run.jags$k[1,,]  
  b.bmsy[i,,] <- biomass[i,,] / run.jags$Bmsy[1,,]  
  indexfit[i,,] <- run.jags$P[i,,] * run.jags$k[1,,] * run.jags$q[1,,]  
  f[i,,] <- sa.dat$sa.catch[i] / biomass[i,,]  
  f.fmsy[i,,] <- f[i,,] / run.jags$Fmsy[1,,]  
}

tab.series <- cbind(sa.dat$year,
             t(apply(biomass[,,1], 1, quantile, qvec)),
             t(apply(b.bmsy[,,1], 1, quantile, qvec)),
             t(apply(f[,,1], 1, quantile, qvec)),
             t(apply(f.fmsy[,,1], 1, quantile, qvec)),
             t(apply(indexfit[,,1], 1, quantile, qvec)))
tab.series <- as.data.frame(tab.series)
names(tab.series) <- c("Year", paste(rep(c("B","B.Bmsy","F","F.Fmsy","Indexfit"), each=length(qvec)), rep(qvec, 5), sep="_"))       

# save the table outside

write.table(tab.series, file=paste("table_series_",run.name,".txt",sep=""), row.names=F)

# years for the plots

years <- sa.dat$year

# plot B, B.Bmsy, F and F.Fmsy

jpeg(paste("series_b_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

jpeg(paste("series_b.bmsy_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B.Bmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B.Bmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B.Bmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]/Bmsy), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
abline(h=1, lty=1, col="grey")
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

jpeg(paste("series_f_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

jpeg(paste("series_f.fmsy_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F.Fmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F.Fmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F.Fmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]/Fmsy), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
abline(h=1, lty=1, col="grey")
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()


jpeg(paste("series_b_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

jpeg(paste("series_b.bmsy_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B.Bmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B.Bmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B.Bmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]/Bmsy), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
abline(h=1, lty=2, col=1)
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

jpeg(paste("series_f_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

jpeg(paste("series_f.fmsy_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F.Fmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F.Fmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F.Fmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]/Fmsy), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
abline(h=1, lty=2, col=1)
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

#------------------------------------------------------------------

# compare observed and fitted

jpeg(paste("series_indexfit_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("Indexfit",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("Indexfit",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("Indexfit",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2,
        ylim=c(0,max( tmp.up, sa.dat$cpue)*1.1))
mtext(expression(Index[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, lwd=2)  
#lines(years, tmp.low, lty=2, lwd=2)  
#lines(years, tmp.up, lty=2, lwd=2)  
points(years, sa.dat$cpue, pch=19, cex=1.5, lwd=1.5)  
dev.off()

#------------------------------------------------------------------

# calculate residuals

res <- array(NA, dim(run.jags$P))
for (i in 1:dim(run.jags$P)[1]){
  res[i,,] <- (log(dat$I[i]) - log(run.jags$q[1,,]) - log(biomass[i,,])) * sqrt(run.jags$psi.logI[1,,])
}

tab.res <- cbind(years, 
              t(apply(res, 1, quantile, qvec)))
tab.res <- as.data.frame(tab.res)
names(tab.res) <- c("Year", paste(rep(c("res.I"), each=length(qvec)), rep(qvec, 1), sep="_"))       

# save the table outside

write.table(tab.res, file=paste("table_residuals_",run.name,".txt",sep=""), row.names=F)

# plot residuals

jpeg(paste("series_res_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 3, 4, 1)+0.1)
tmp.low <- tab.res[, paste("res.I",qvec[1], sep="_")]
tmp.med <- tab.res[, paste("res.I",qvec[2], sep="_")]
tmp.up <- tab.res[, paste("res.I",qvec[3], sep="_")]
matplot(years, data.frame(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
title("Residuals")
abline(h=c(-2,0,2), lty=c(2,1,2), col="grey")
segments(years, tmp.low, years, tmp.up, lwd=1.5)
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

#------------------------------------------------------------------

# Bmsy, Fmsy and MSY densities

jpeg(paste("pdf_msy_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,3))
plot(density(run.mcmc[,"Bmsy"], from=quantile(run.mcmc[,"Bmsy"],0.01), to=quantile(run.mcmc[,"Bmsy"],0.99)), main="Bmsy", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
plot(density(run.mcmc[,"Fmsy"], from=quantile(run.mcmc[,"Fmsy"],0.01), to=quantile(run.mcmc[,"Fmsy"],0.99)), main="Fmsy", xlab="", ylab="", cex.lab=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
plot(density(run.mcmc[,"MSY"], from=quantile(run.mcmc[,"MSY"],0.01), to=quantile(run.mcmc[,"MSY"],0.99)), main="MSY", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

# CV's of process error (on normal scale from the psi.logP)

aux1 <- sqrt(exp(1/run.mcmc[,"psi.logP"])-1)
summary(aux1)
plot(density(aux1, from=0, to=1), xlab="", ylab="", main="CV process error")

# CV's of observation error (on normal scale from the psi.logI)

aux2 <- sqrt(exp(1/run.mcmc[,"psi.logI"])-1)
summary(aux2)
plot(density(aux2, from=0, to=1), xlab="", ylab="", main="CV observation error")

jpeg(paste("pdf_cv_errors_",run.name,".jpg",sep=""), width=700, quality=100)
plot(density(aux1, from=0, to=1), xlab="", ylab="", main="CV errors")
lines(density(aux2, from=0, to=1), col=2)
legend("topright", c("process","observation"), lty=1, col=1:2)
dev.off()

#------------------------------------------------------------------

# kobe plot: stock trajectory

ymax <- max(3, tab.series$F.Fmsy_0.5)+0.1 # maximum for F/Fmsy
xmax <- max(3, tab.series$B.Bmsy_0.5)+0.1 # maximum for B/Bmsy

# define colors with transparency
# see also: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

col.yel <- as.vector(col2rgb(7))/255
col.yel <- rgb(col.yel[1], col.yel[2], col.yel[3], alpha=0.5)
col.red <- as.vector(col2rgb(2))/255
col.red <- rgb(col.red[1], col.red[2], col.red[3], alpha=0.5)
col.gre <- as.vector(col2rgb(3))/255
col.gre <- rgb(col.gre[1], col.gre[2], col.gre[3], alpha=0.5)


jpeg(paste("kobe_trajectory_",run.name,".jpg"), width=700, height=700, quality=100)
plot(c(0,xmax), c(0,ymax), type="n", xlab="B/Bmsy", ylab="F/Fmsy")
rect(0,0,1,1, col=col.yel)
rect(1,0,xmax,1, col=col.gre)
rect(0,1,1,ymax, col=col.red)
rect(1,1,xmax,ymax, col=col.yel)
lines(tab.series$B.Bmsy_0.5, tab.series$F.Fmsy_0.5, type="o", cex=2, pch=1)
points(tab.series$B.Bmsy_0.5[c(1,n)], tab.series$F.Fmsy_0.5[c(1,n)], cex=2, pch=19)
textxy(tab.series$B.Bmsy_0.5[c(1,n)], tab.series$F.Fmsy_0.5[c(1,n)], tab.series$Year[c(1,n)], cex=1)
dev.off()

# kobe plot contour plot last year

jpeg(paste("kobe_lastyear_",run.name,".jpg",sep=""), width=700, height=700, quality=100)
plot(c(0,xmax), c(0,ymax), type="n", xlab="B/Bmsy", ylab="F/Fmsy")
rect(0,0,1,1, col=col.yel)
rect(1,0,xmax,1, col=col.gre)
rect(0,1,1,ymax, col=col.red)
rect(1,1,xmax,ymax, col=col.yel)
f1 <- kde2d(as.vector(b.bmsy[n,,1]), as.vector(f.fmsy[n,,1] ), n=50)
contour(f1, nlevels=5, add=T)
points(tab.series$B.Bmsy_0.5[n], tab.series$F.Fmsy_0.5[n], cex=2, pch=19)
dev.off()

save.image("main.RData")

####################################################################
####################################################################

# run identifier:

run.name <- "run2"

# first set of prior distributions

mu.logr <- log(0.3)   # based on r=2*Fmsey=2*0.87*M (from Zhou et al 2012) with M=0.19 from Grandcourt 2005
psi.logr <- 5         
mu.logk <-  log(2000) 
psi.logk <- 5
mu.logq <- log(1e-7)
psi.logq <- 5
mu.logQ <- mu.logq + mu.logk  # reparameterization. the mean is the sum of means and variance is the sum of variances
psi.logQ <- (psi.logq + psi.logk) /(psi.logq*psi.logk)
a.psi.logI <- 4          # 1/sqrt(a) is the CV of the gamma prior distribution
b.psi.logI <- 0.16
a.psi.logP <- 4
b.psi.logP <- 0.01 
a.P0 <- 0
b.P0 <- 1
mu.logp <- log(1)
psi.logp <- 1

calc.biomass(r=exp(mu.logr),k=exp(mu.logk),p=1,P0=0.4,var.logP=b.psi.logP/a.psi.logP,C=sa.dat$sa.catch)
calc.nuisance(r=exp(mu.logr),k=exp(mu.logk),p=1,P0=0.4,Index=sa.dat$cpue, C=sa.dat$sa.catch)

# vector to compute quantiles

qvec <- c(0.05, 0.1, 0.5, 0.9, 0.95)

tab.priors <- NULL
tab.priors <- rbind(tab.priors, c(mu.logr, psi.logr, qlnorm(qvec, mu.logr, 1/sqrt(psi.logr))))
tab.priors <- rbind(tab.priors, c(mu.logk, psi.logk, qlnorm(qvec, mu.logk, 1/sqrt(psi.logk))))
tab.priors <- rbind(tab.priors, c(mu.logq, psi.logq, qlnorm(qvec, mu.logq, 1/sqrt(psi.logq))))
tab.priors <- rbind(tab.priors, c(a.P0, b.P0, qunif(qvec, a.P0, b.P0)))
tab.priors <- rbind(tab.priors, c(mu.logp, psi.logp, qlnorm(qvec, mu.logp, 1/sqrt(psi.logp))))
tab.priors <- rbind(tab.priors, c(a.psi.logI, b.psi.logI, qgamma(qvec, a.psi.logI, b.psi.logI)))
tab.priors <- rbind(tab.priors, c(a.psi.logP, b.psi.logP, qgamma(qvec, a.psi.logP, b.psi.logP)))

tab.priors <- cbind(parameter=c("r","k","q","P0","p","psi.logI","psi.logP"),as.data.frame(tab.priors))

# save the prior table outside

write.table(tab.priors, file=paste("table_priors_",run.name,".txt",sep=""), row.names=F, col.names=F)

# The gamma distribution has mean=a/b, var=a/b^2, cv=1/sqrt(a) 

qgamma(qvec, a.psi.logI, b.psi.logI)   

# this leads to CV of Index: 

sqrt(exp(1/qgamma(qvec, a.psi.logI, b.psi.logI))-1)

#psi.logP

qgamma(qvec, a.psi.logP, b.psi.logP)   

# this leads to CV of process equations: 

sqrt(exp(1/qgamma(qvec, a.psi.logP, b.psi.logP))-1)

#------------------------------------------------------------------

# especificar longitudes de MCMC. se pueden dar más abajo al aplicar las funciones directamente

n.chains <- 1
n.adapt <- 10000
n.burnin <- 100000
n.thin <- 100
n.iter <- 1000000

#------------------------------------------------------------------

# write data to a txt file

dat <- list(n=n, I=sa.dat$cpue, C=sa.dat$sa.catch,
            mu.logr=mu.logr, psi.logr=psi.logr,
            mu.logk=mu.logk, psi.logk=psi.logk,
            mu.logQ=mu.logQ, psi.logQ=psi.logQ,
            #mu.logq=mu.logq, psi.logq=psi.logq,
            a.P0=a.P0, b.P0=b.P0,                           
            a.psi.logI=a.psi.logI, b.psi.logI=b.psi.logI,
            a.psi.logP=a.psi.logP, b.psi.logP=b.psi.logP)

# default starting values in JAGS

#inits <- NULL

# generate initial values that fulfill the catch restrictions

inits <- pellatom.inits(C=sa.dat$sa.catch,
                        mu.logr=mu.logr, psi.logr=psi.logr,
                        mu.logk=mu.logk, psi.logk=psi.logk,
                        mu.logq=mu.logq, psi.logq=psi.logq,
                        mu.logp=mu.logp, psi.logp=psi.logp,
                        a.P0=a.P0, b.P0=b.P0,                     
                        a.psi.logI=a.psi.logI, b.psi.logI=b.psi.logI,
                        a.psi.logP=a.psi.logP, b.psi.logP=b.psi.logP,
                        fx.logp=T, fx.P0=F,
                        n.chains = n.chains, maxit = 100)

# modify to include logQ

for (i in 1:n.chains){
  inits[[i]]$logQ <- inits[[i]]$logq + inits[[i]]$logk  
}

# list of parameters to be saved

param <- c("logr", "logk", "logQ", "P0", "psi.logI", "psi.logP", "r", "k", "Q", "q", "P", "Fmsy", "Bmsy", "MSY", "dif")

# path for the JAGS model

model.path <- ".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/code_surplus/s1.jags"

ptm.ini <- Sys.time()
print(Sys.time())
m <- jags.model(model.path, dat=dat, inits=inits, n.chains=n.chains, n.adapt=n.adapt)
print(Sys.time())

# burn-in

print(Sys.time())
update(m, n.iter=n.burnin)
print(Sys.time())

# jags object with MCMC

run.jags <- jags.samples(m, variable.names=param, n.iter=n.iter, thin=n.thin)
ptm.end <- Sys.time()
ptm <- ptm.end - ptm.ini
print(ptm)

# Time difference of 2.504008 mins for 1 chain

# save results externally

save(run.jags, file=paste("res_bay_",run.name,".RData",sep=""))

#------------------------------------------------------------------

# check that the max function did not have any effect. look at min and max values. they should be 0

apply(run.jags$dif, 1, quantile, c(0,0.5,1))

#------------------------------------------------------------------

# change from the jags.samples format to coda format

run.mcmc <- jags2coda(run.jags)

# additional selection

#run.mcmc <- window(run.mcmc, start=500)

# full list of names of parameters to be estimated

param.names <- c("logr", "logk", "logQ", "P0",
                 "psi.logI", "psi.logP",
                 "r","k","Q","q",   
                 paste("P[",1:n,"]",sep=""),
                 "Fmsy", "Bmsy","MSY")

param.names <- param.names[param.names %in% varnames(run.mcmc)] # to filter in case some of the parameters are missing (e.g. fixed)

# re-order the mcmc object according to the param.names vector

run.mcmc <- run.mcmc[,param.names]

#summary
summary(run.mcmc)

# check results from 2 chains

# cbind(summary(run.mcmc[[2]])$quantiles[, c(1,3,5)], summary(run.mcmc[[1]])$quantiles[,c(1,3,5)])

# check standard plots made with coda library

pdf(paste("plot_coda_",run.name,".pdf",sep=""))
plot(run.mcmc)
dev.off()

# select only one chain

run.mcmc <- run.mcmc[[1]]

# vector of quantiles to use to get the credible intervals and the estimates

qvec <- c(0.05, 0.5, 0.95)

# summary statistics using coda

summary(run.mcmc)
write.table(summary(run.mcmc, quantiles=qvec)$quantiles, file=paste("table_summary_parameters_",run.name,".txt",sep=""), col.names=F)

# check plots

plot.chain(run.mcmc, filename=paste("plot_",run.name,".pdf",sep=""))

#------------------------------------------------------------------

# pdf with figures for comparison prior and posterior

pdf(paste("prior_posterior_",run.name,".pdf",sep=""))

# log(r)
xmin <- min(qnorm(c(0.025, 0.975), mu.logr, sqrt(1/psi.logr)), run.jags$logr)
xmax <- max(qnorm(c(0.025, 0.975), mu.logr, sqrt(1/psi.logr)), run.jags$logr)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$logr, bw=(xmax-xmin)/100)
yprior <- dnorm(x, mu.logr, sqrt(1/psi.logr))
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main="log(r)", yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lty=2, lwd=2)

# log(k)
xmin <- min(qnorm(c(0.025, 0.975), mu.logk, sqrt(1/psi.logk)), run.jags$logk)
xmax <- max(qnorm(c(0.025, 0.975), mu.logk, sqrt(1/psi.logk)), run.jags$logk)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$logk, bw=(xmax-xmin)/100)
yprior <- dnorm(x, mu.logk, sqrt(1/psi.logk))
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main="log(k)", yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lty=2, lwd=2)

# log(q)
xmin <- min(qnorm(c(0.025, 0.975), mu.logq, sqrt(1/psi.logq)), log(run.jags$q))
xmax <- max(qnorm(c(0.025, 0.975), mu.logq, sqrt(1/psi.logq)), log(run.jags$q))
x <- seq(xmin, xmax, length=10000)
y <- density(log(run.jags$q), bw=(xmax-xmin)/100)
yprior <- dnorm(x, mu.logq, sqrt(1/psi.logq))
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main="log(q)", yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lty=2, lwd=2)

# psi.logI

xmin <- min(qgamma(c(0.025, 0.975), a.psi.logI, b.psi.logI), run.jags$psi.logI)
xmax <- max(qgamma(c(0.025, 0.975), a.psi.logI, b.psi.logI), run.jags$psi.logI)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$psi.logI, bw=(xmax-xmin)/100)
yprior <- dgamma(x, a.psi.logI, b.psi.logI)
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
ymin <- 0
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=expression(psi[log(I)]), yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lwd=2, lty=2)

# psi.logP

xmin <- min(qgamma(c(0.025, 0.975), a.psi.logP, b.psi.logP), run.jags$psi.logP)
xmax <- max(qgamma(c(0.025, 0.975), a.psi.logP, b.psi.logP), run.jags$psi.logP)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$psi.logP, bw=(xmax-xmin)/100)
yprior <- dgamma(x, a.psi.logP, b.psi.logP)
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
ymin <- 0
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=expression(psi[log(P)]), yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lwd=2, lty=2)

# P0 prior and posterior

xmin <- min(qunif(c(0.025, 0.975), a.P0, b.P0), run.jags$P0)
xmax <- max(qunif(c(0.025, 0.975), a.P0, b.P0), run.jags$P0)
x <- seq(xmin, xmax, length=10000)
y <- density(run.jags$P0, bw=(xmax-xmin)/50)
yprior <- dunif(x, a.P0, b.P0)
ymin  <- min(y$y, yprior)
ymax  <- max(y$y, yprior)
ymin <- 0
plot(y, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="", main=expression(P[0]), yaxt="n", cex.main=2.5, cex.lab=1.5, cex.axis=1.5)
lines(x, yprior, lwd=2, lty=2)

dev.off()

#------------------------------------------------------------------

# calculate B, B/Bmsy, F, F/Fmsy, Indexfit time series

biomass <- array(NA, dim(run.jags$P))
b.bmsy <- array(NA, dim(run.jags$P))
indexfit <- array(NA, dim(run.jags$P))
f <- array(NA, dim(run.jags$P))
f.fmsy <- array(NA, dim(run.jags$P))
for (i in 1:dim(run.jags$P)[1]){
  biomass[i,,] <- run.jags$P[i,,] * run.jags$k[1,,]  
  b.bmsy[i,,] <- biomass[i,,] / run.jags$Bmsy[1,,]  
  indexfit[i,,] <- run.jags$P[i,,] * run.jags$k[1,,] * run.jags$q[1,,]  
  f[i,,] <- sa.dat$sa.catch[i] / biomass[i,,]  
  f.fmsy[i,,] <- f[i,,] / run.jags$Fmsy[1,,]  
}

tab.series <- cbind(sa.dat$year,
                    t(apply(biomass[,,1], 1, quantile, qvec)),
                    t(apply(b.bmsy[,,1], 1, quantile, qvec)),
                    t(apply(f[,,1], 1, quantile, qvec)),
                    t(apply(f.fmsy[,,1], 1, quantile, qvec)),
                    t(apply(indexfit[,,1], 1, quantile, qvec)))
tab.series <- as.data.frame(tab.series)
names(tab.series) <- c("Year", paste(rep(c("B","B.Bmsy","F","F.Fmsy","Indexfit"), each=length(qvec)), rep(qvec, 5), sep="_"))       

# save the table outside

write.table(tab.series, file=paste("table_series_",run.name,".txt",sep=""), row.names=F)

# years for the plots

years <- sa.dat$year

# plot B, B.Bmsy, F and F.Fmsy

jpeg(paste("series_b_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

jpeg(paste("series_b.bmsy_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B.Bmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B.Bmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B.Bmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]/Bmsy), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
abline(h=1, lty=1, col="grey")
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

jpeg(paste("series_f_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

jpeg(paste("series_f.fmsy_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F.Fmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F.Fmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F.Fmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]/Fmsy), side=2, line=3, cex=1.5)
segments(years, tmp.low, years, tmp.up, lwd=1.5)
abline(h=1, lty=1, col="grey")
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()


jpeg(paste("series_b_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

jpeg(paste("series_b.bmsy_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B.Bmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B.Bmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B.Bmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(B[y]/Bmsy), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
abline(h=1, lty=2, col=1)
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

jpeg(paste("series_f_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

jpeg(paste("series_f.fmsy_lines_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F.Fmsy",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F.Fmsy",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F.Fmsy",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
mtext(expression(F[y]/Fmsy), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
abline(h=1, lty=2, col=1)
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
dev.off()

#------------------------------------------------------------------

# compare observed and fitted

jpeg(paste("series_indexfit_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("Indexfit",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("Indexfit",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("Indexfit",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2,
        ylim=c(0,max( tmp.up, sa.dat$cpue)*1.1))
mtext(expression(Index[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, lwd=2)  
#lines(years, tmp.low, lty=2, lwd=2)  
#lines(years, tmp.up, lty=2, lwd=2)  
points(years, sa.dat$cpue, pch=19, cex=1.5, lwd=1.5)  
dev.off()

#------------------------------------------------------------------

# calculate residuals

res <- array(NA, dim(run.jags$P))
for (i in 1:dim(run.jags$P)[1]){
  res[i,,] <- (log(dat$I[i]) - log(run.jags$q[1,,]) - log(biomass[i,,])) * sqrt(run.jags$psi.logI[1,,])
}

tab.res <- cbind(years, 
                 t(apply(res, 1, quantile, qvec)))
tab.res <- as.data.frame(tab.res)
names(tab.res) <- c("Year", paste(rep(c("res.I"), each=length(qvec)), rep(qvec, 1), sep="_"))       

# save the table outside

write.table(tab.res, file=paste("table_residuals_",run.name,".txt",sep=""), row.names=F)

# plot residuals

jpeg(paste("series_res_bars_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 3, 4, 1)+0.1)
tmp.low <- tab.res[, paste("res.I",qvec[1], sep="_")]
tmp.med <- tab.res[, paste("res.I",qvec[2], sep="_")]
tmp.up <- tab.res[, paste("res.I",qvec[3], sep="_")]
matplot(years, data.frame(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2)
title("Residuals")
abline(h=c(-2,0,2), lty=c(2,1,2), col="grey")
segments(years, tmp.low, years, tmp.up, lwd=1.5)
points(years, tmp.med, pch=19, cex=1.5, lwd=1.5)  
dev.off()

#------------------------------------------------------------------

# Bmsy, Fmsy and MSY densities

jpeg(paste("pdf_msy_",run.name,".jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,3))
plot(density(run.mcmc[,"Bmsy"], from=quantile(run.mcmc[,"Bmsy"],0.01), to=quantile(run.mcmc[,"Bmsy"],0.99)), main="Bmsy", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
plot(density(run.mcmc[,"Fmsy"], from=quantile(run.mcmc[,"Fmsy"],0.01), to=quantile(run.mcmc[,"Fmsy"],0.99)), main="Fmsy", xlab="", ylab="", cex.lab=1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
plot(density(run.mcmc[,"MSY"], from=quantile(run.mcmc[,"MSY"],0.01), to=quantile(run.mcmc[,"MSY"],0.99)), main="MSY", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

# CV's of process error (on normal scale from the psi.logP)

aux1 <- sqrt(exp(1/run.mcmc[,"psi.logP"])-1)
summary(aux1)
plot(density(aux1, from=0, to=1), xlab="", ylab="", main="CV process error")

# CV's of observation error (on normal scale from the psi.logI)

aux2 <- sqrt(exp(1/run.mcmc[,"psi.logI"])-1)
summary(aux2)
plot(density(aux2, from=0, to=1), xlab="", ylab="", main="CV observation error")

jpeg(paste("pdf_cv_errors_",run.name,".jpg",sep=""), width=700, quality=100)
plot(density(aux1, from=0, to=1), xlab="", ylab="", main="CV errors")
lines(density(aux2, from=0, to=1), col=2)
legend("topright", c("process","observation"), lty=1, col=1:2)
dev.off()

#------------------------------------------------------------------

# kobe plot: stock trajectory

ymax <- max(3, tab.series$F.Fmsy_0.5)+0.1 # maximum for F/Fmsy
xmax <- max(3, tab.series$B.Bmsy_0.5)+0.1 # maximum for B/Bmsy

# define colors with transparency
# see also: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

col.yel <- as.vector(col2rgb(7))/255
col.yel <- rgb(col.yel[1], col.yel[2], col.yel[3], alpha=0.5)
col.red <- as.vector(col2rgb(2))/255
col.red <- rgb(col.red[1], col.red[2], col.red[3], alpha=0.5)
col.gre <- as.vector(col2rgb(3))/255
col.gre <- rgb(col.gre[1], col.gre[2], col.gre[3], alpha=0.5)


jpeg(paste("kobe_trajectory_",run.name,".jpg"), width=700, height=700, quality=100)
plot(c(0,xmax), c(0,ymax), type="n", xlab="B/Bmsy", ylab="F/Fmsy")
rect(0,0,1,1, col=col.yel)
rect(1,0,xmax,1, col=col.gre)
rect(0,1,1,ymax, col=col.red)
rect(1,1,xmax,ymax, col=col.yel)
lines(tab.series$B.Bmsy_0.5, tab.series$F.Fmsy_0.5, type="o", cex=2, pch=1)
points(tab.series$B.Bmsy_0.5[c(1,n)], tab.series$F.Fmsy_0.5[c(1,n)], cex=2, pch=19)
textxy(tab.series$B.Bmsy_0.5[c(1,n)], tab.series$F.Fmsy_0.5[c(1,n)], tab.series$Year[c(1,n)], cex=1)
dev.off()

# kobe plot contour plot last year

jpeg(paste("kobe_lastyear_",run.name,".jpg",sep=""), width=700, height=700, quality=100)
plot(c(0,xmax), c(0,ymax), type="n", xlab="B/Bmsy", ylab="F/Fmsy")
rect(0,0,1,1, col=col.yel)
rect(1,0,xmax,1, col=col.gre)
rect(0,1,1,ymax, col=col.red)
rect(1,1,xmax,ymax, col=col.yel)
f1 <- kde2d(as.vector(b.bmsy[n,,1]), as.vector(f.fmsy[n,,1] ), n=50)
contour(f1, nlevels=5, add=T)
points(tab.series$B.Bmsy_0.5[n], tab.series$F.Fmsy_0.5[n], cex=2, pch=19)
dev.off()

save.image("main.RData")

####################################################################
####################################################################
