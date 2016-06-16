# Short term forecast from surplus production model
# BayesianSchaeffer/examples/code_stf_sa_eni.R 

# Copyright 2015 AZTI Team. Distributed under the GPL 2 or later 
# Maintainer: Leire Ibaibarriaga (libaibarriaga@azti.es) 
# $Id: code_stf_sa_eni.R 24/09/2015 libaibarriaga $ 

####################################################################
####################################################################
#
# Orange spotted grouper (Epinephelus coioides)      - ENI - AZTI
# Catches from Saudi Arabia
#
####################################################################
####################################################################

# species

sp <- "ENI"

####################################################################
####################################################################

# set working directory

wd <- paste(".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/",tolower(sp),"/stf_sa/",sep="")

setwd(wd)

####################################################################
####################################################################

# load library for JAGS

library(rjags)

####################################################################
####################################################################

# load complementary functions

source(".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/code_surplus/functions_pellatom.r")

####################################################################
####################################################################

# load data file

sa.dat <- read.table(paste(".../IM12KFUPM/WP5.2_PopDyn&StockAss/data_analysis/data_",sp,".txt",sep=""), header=T)

# subset without NA's

sa.dat <- subset(sa.dat, !is.na(sa.dat$cpue) & !is.na(sa.dat$sa.catch))

# number of years

n <- dim(sa.dat)[1]

####################################################################
####################################################################

# load mcmc run results in jags format

run.name <- 'run2'
load(paste(".../IM12KFUPM/WP5.2_PopDyn&StockAss/assessments_by_species/",tolower(sp),"/surplus/","res_bay_",run.name,".RData",sep=""))

# change from the jags.samples format to coda format

run.mcmc <- jags2coda(run.jags)
run.mcmc <- run.mcmc[[1]]

# full list of names of parameters to be estimated

param.names <- c("logr", "logk", "logQ", "P0",
                 "psi.logI", "psi.logP",
                 "r","k","Q","q",   
                 paste("P[",1:n,"]",sep=""),
                 "Fmsy", "Bmsy","MSY")

param.names <- param.names[param.names %in% varnames(run.mcmc)] # to filter in case some of the parameters are missing (e.g. fixed)

# re-order the mcmc object according to the param.names vector

run.mcmc <- run.mcmc[,param.names]

# vector of quantiles to use to get the credible intervals and the estimates

qvec <- c(0.05, 0.5, 0.95)

# calculate B, B/Bmsy, F, F/Fmsy, Indexfit time series

biomass <- array(NA, dim(run.jags$P))
bmsy <- array(NA, dim(run.jags$P))
b.bmsy <- array(NA, dim(run.jags$P))
indexfit <- array(NA, dim(run.jags$P))
f <- array(NA, dim(run.jags$P))
fmsy <- array(NA, dim(run.jags$P))
f.fmsy <- array(NA, dim(run.jags$P))
for (i in 1:dim(run.jags$P)[1]){
  biomass[i,,] <- run.jags$P[i,,] * run.jags$k[1,,]  
  bmsy[i,,] <- run.jags$Bmsy[1,,]  
  b.bmsy[i,,] <- biomass[i,,] / run.jags$Bmsy[1,,]  
  indexfit[i,,] <- run.jags$P[i,,] * run.jags$k[1,,] * run.jags$q[1,,]  
  f[i,,] <- sa.dat$sa.catch[i] / biomass[i,,] 
  fmsy[i,,] <- run.jags$Fmsy[1,,]
  f.fmsy[i,,] <- f[i,,] / run.jags$Fmsy[1,,]  
}

tab.series <- cbind(sa.dat$year,
                    sa.dat$sa.catch,
                    t(apply(biomass[,,1], 1, quantile, qvec)),
                    t(apply(bmsy[,,1], 1, quantile, qvec)),
                    t(apply(b.bmsy[,,1], 1, quantile, qvec)),
                    t(apply(f[,,1], 1, quantile, qvec)), #! ERROR
                    t(apply(fmsy[,,1], 1, quantile, qvec)),
                    t(apply(f.fmsy[,,1], 1, quantile, qvec)))
tab.series <- as.data.frame(tab.series)
names(tab.series) <- c("Year", "Catch", paste(rep(c("B","Bmsy","B.Bmsy","F","Fmsy","F.Fmsy"), each=length(qvec)), rep(qvec, 4), sep="_"))       

####################################################################
####################################################################

# first we need to move the population to the end of the last year, so that it can be used as starting point

Cval <- sa.dat$sa.catch[n]
nproj <- 1

aux <- run.mcmc[,c("r","k",paste("P[",n,"]",sep=""),"psi.logP")]
fun <- function(vec, Cval){
  proj.biomass(r=vec[1], k=vec[2], p=1, Pini=vec[3], var.logP=1/vec[4], C=Cval, type="C", nproj=nproj)
}

proj.mat <- matrix(unlist(apply(aux, 1, fun, Cval=Cval)), byrow=T, ncol=3)
proj.mat <- as.data.frame(proj.mat)
names(proj.mat) <- paste(rep(c("B","C","F"),each=nproj), 1:nproj, sep="")

Bini <- proj.mat$B1 # biomass at the beginning of the projection period (or at the end of the assessment year)
Pini <- Bini/run.jags$k[1,,]
  
####################################################################
####################################################################

# NOTE THAT THE FUNCTION proj.biomass() gives B at the beginning of year y+1 and C (or F) during year y
# projection for a given F value
# current F

Fval <- Fcurr <- median(f[n,,])
nproj  <- 2

aux <- cbind(run.mcmc[,c("r","k","psi.logP")], Pini)
fun <- function(vec, Fval){
  proj.biomass(r=vec[1], k=vec[2], p=1, Pini=vec[4], var.logP=1/vec[3], F=Fval, type="F", nproj=nproj)
}

proj.mat <- matrix(unlist(apply(aux, 1, fun, Fval=Fval)), byrow=T, ncol=3*2)
proj.mat <- cbind(Bini, proj.mat)
proj.mat <- as.data.frame(proj.mat)
names(proj.mat) <- c("B0", paste(rep(c("B","C","F"),each=nproj), 1:nproj, sep=""))

proj.mat <- cbind(proj.mat,
                  proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]/run.jags$Bmsy[1,,],
                  proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]/run.jags$Fmsy[1,,],
                  (proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]-Bini)/Bini,
                  (proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]-f[n,,])/f[n,,])
names(proj.mat) <- c("B0", paste(rep(c("B","C","F","B.Bmsy","F.Fmsy","B.change","F.change"),each=2), 1:nproj, sep=""))

# quantiles

proj.q <- apply(proj.mat, 2, quantile, qvec)

# years for the plots

years <- sa.dat$year

# plot

jpeg(paste("series_b_lines_",run.name,"_Fcurr.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj+1))
mtext(expression(B[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
tmp <- t(proj.q[,match(paste("B",0:nproj,sep=""), colnames(proj.q))])
tmp <- rbind(matrix(c(tmp.low[n], tmp.med[n], tmp.up[n]),nrow=1), tmp)
matplot(years[n]+(0:(nproj+1)), tmp, add=T, type="o", pch=c(NA,1,NA), lty=c(2,1,2), lwd=2, col=1, cex=1.5)
msy.low <- tab.series[, paste("Bmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Bmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Bmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

jpeg(paste("series_f_lines_",run.name,"_Fcurr.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj))
mtext(expression(F[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
tmp <- proj.q[2,match(paste("F",1:nproj,sep=""), colnames(proj.q))]
tmp <- c(tmp.med[n], tmp)
matplot(years[n]+(0:nproj), tmp, add=T, type="o", pch=1, lty=1, lwd=2, col=1, cex=1.5)
msy.low <- tab.series[, paste("Fmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Fmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Fmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

# we could do also plots for Bmsy and Fmsy

# Table of results:

write.table(matrix(proj.q[qvec==0.5,], nrow=1), file="catch-table.txt", append=F, row.names=F, col.names=colnames(proj.q))

####################################################################
####################################################################

# projection for a given C value
# current C

Cval <- Ccurr <- sa.dat$sa.catch[n]
nproj  <- 2

aux <- cbind(run.mcmc[,c("r","k","psi.logP")], Pini)
fun <- function(vec, Cval){
  proj.biomass(r=vec[1], k=vec[2], p=1, Pini=vec[4], var.logP=1/vec[3], C=Cval, type="C", nproj=nproj)
}

proj.mat <- matrix(unlist(apply(aux, 1, fun, Cval=Cval)), byrow=T, ncol=3*2)
proj.mat <- cbind(Bini, proj.mat)
proj.mat <- as.data.frame(proj.mat)
names(proj.mat) <- c("B0", paste(rep(c("B","C","F"),each=nproj), 1:nproj, sep=""))

proj.mat <- cbind(proj.mat,
                  proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]/run.jags$Bmsy[1,,],
                  proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]/run.jags$Fmsy[1,,],
                  (proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]-Bini)/Bini,
                  (proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]-f[n,,])/f[n,,])
names(proj.mat) <- c("B0", paste(rep(c("B","C","F","B.Bmsy","F.Fmsy","B.change","F.change"),each=2), 1:nproj, sep=""))

# quantiles

proj.q <- apply(proj.mat, 2, quantile, qvec)

# years for the plots

years <- sa.dat$year

# plot

jpeg(paste("series_b_lines_",run.name,"_Ccurr.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj+1))
mtext(expression(B[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
tmp <- t(proj.q[,match(paste("B",0:nproj,sep=""), colnames(proj.q))])
tmp <- rbind(matrix(c(tmp.low[n], tmp.med[n], tmp.up[n]),nrow=1), tmp)
matplot(years[n]+(0:(nproj+1)), tmp, add=T, type="o", pch=c(NA,1,NA), lty=c(2,1,2), lwd=2, col=1, cex=1.5)
msy.low <- tab.series[, paste("Bmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Bmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Bmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

jpeg(paste("series_f_lines_",run.name,"_Ccurr.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj))
mtext(expression(F[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
tmp <- proj.q[2,match(paste("F",1:nproj,sep=""), colnames(proj.q))]
tmp <- c(tmp.med[n], tmp)
matplot(years[n]+(0:nproj), tmp, add=T, type="o", pch=1, lty=1, lwd=2, col=1, cex=1.5)
msy.low <- tab.series[, paste("Fmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Fmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Fmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

# we could do also plots for Bmsy and Fmsy

# Table of results:

write.table(matrix(proj.q[qvec==0.5,], nrow=1), file="catch-table.txt", append=T, row.names=F, col.names=F)

####################################################################
####################################################################

# projection for a range of F values for the catch options

Fvector <- c(0, 
             Fcurr, 
             median(fmsy),
             Fcurr*0.75,
             Fcurr*1.25)

nproj  <- 2

aux <- cbind(run.mcmc[,c("r","k","psi.logP")], Pini)
fun <- function(vec, Fval){
  proj.biomass(r=vec[1], k=vec[2], p=1, Pini=vec[4], var.logP=1/vec[3], F=Fval, type="F", nproj=nproj)
}

proj.all <- NULL
for (i in 1:length(Fvector)){
  Fval <- Fvector[i]
  proj.mat <- matrix(unlist(apply(aux, 1, fun, Fval=Fval)), byrow=T, ncol=3*2)
  proj.mat <- cbind(Bini, proj.mat)
  proj.mat <- as.data.frame(proj.mat)
  names(proj.mat) <- c("B0", paste(rep(c("B","C","F"),each=nproj), 1:nproj, sep=""))  
  proj.mat <- cbind(proj.mat,
                    proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]/run.jags$Bmsy[1,,],
                    proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]/run.jags$Fmsy[1,,],
                    (proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]-Bini)/Bini,
                    (proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]-f[n,,])/f[n,,])
  names(proj.mat) <- c("B0", paste(rep(c("B","C","F","B.Bmsy","F.Fmsy","B.change","F.change"),each=2), 1:nproj, sep=""))
  proj.q <- apply(proj.mat, 2, quantile, qvec)
  proj.all <- rbind(proj.all, proj.q[2,])
}

jpeg(paste("series_b_lines_",run.name,"_Frange.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj+1))
mtext(expression(B[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
for (i in 1:length(Fvector)){
  tmp <- proj.all[i, match(paste("B",0:nproj,sep=""), colnames(proj.all))]
  tmp <- c(tmp.med[n], tmp)
  matplot(years[n]+(0:(nproj+1)), tmp, add=T, type="l", lty=1, lwd=2, col=1, cex=1.5)
}
msy.low <- tab.series[, paste("Bmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Bmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Bmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

jpeg(paste("series_f_lines_",run.name,"_Frange.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj), ylim=c(min(tmp.low,Fvector),max(tmp.up,Fvector)))
mtext(expression(F[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)
for (i in 1:length(Fvector)){
  tmp <- proj.all[i, match(paste("F",1:nproj,sep=""), colnames(proj.all))]
  tmp <- c(tmp.med[n], tmp)
  matplot(years[n]+(0:nproj), tmp, add=T, type="l", lty=1, lwd=2, col=1, cex=1.5)
}
msy.low <- tab.series[, paste("Fmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Fmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Fmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

# Table of results:

write.table(proj.all, file="catch-table.txt", append=T, row.names=F, col.names=F)


####################################################################
####################################################################

# projection for a range of F values: from 0 to 1 by 0.1 

Fvector <- seq(0, 1, by=0.1)

nproj  <- 2

aux <- cbind(run.mcmc[,c("r","k","psi.logP")], Pini)
fun <- function(vec, Fval){
  proj.biomass(r=vec[1], k=vec[2], p=1, Pini=vec[4], var.logP=1/vec[3], F=Fval, type="F", nproj=nproj)
}

proj.all <- NULL
for (i in 1:length(Fvector)){
  Fval <- Fvector[i]
  proj.mat <- matrix(unlist(apply(aux, 1, fun, Fval=Fval)), byrow=T, ncol=3*2)
  proj.mat <- cbind(Bini, proj.mat)
  proj.mat <- as.data.frame(proj.mat)
  names(proj.mat) <- c("B0", paste(rep(c("B","C","F"),each=nproj), 1:nproj, sep=""))  
  proj.mat <- cbind(proj.mat,
                    proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]/run.jags$Bmsy[1,,],
                    proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]/run.jags$Fmsy[1,,],
                    (proj.mat[, match(paste("B",1:nproj,sep=""), names(proj.mat))]-Bini)/Bini,
                    (proj.mat[, match(paste("F",1:nproj,sep=""), names(proj.mat))]-f[n,,])/f[n,,])
  names(proj.mat) <- c("B0", paste(rep(c("B","C","F","B.Bmsy","F.Fmsy","B.change","F.change"),each=2), 1:nproj, sep=""))
  proj.q <- apply(proj.mat, 2, quantile, qvec)
  proj.all <- rbind(proj.all, proj.q[2,])
}

jpeg(paste("series_b_lines_",run.name,"_Frange_by_0.1.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("B",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("B",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("B",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj+1))
mtext(expression(B[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)  
for (i in 1:length(Fvector)){
  tmp <- proj.all[i, match(paste("B",0:nproj,sep=""), colnames(proj.all))]
  tmp <- c(tmp.med[n], tmp)
  matplot(years[n]+(0:(nproj+1)), tmp, add=T, type="l", lty=1, lwd=2, col=1, cex=1.5)
}
msy.low <- tab.series[, paste("Bmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Bmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Bmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

jpeg(paste("series_f_lines_",run.name,"_Frange_by_0.1.jpg",sep=""), width=700, quality=100)
par(mfrow=c(1,1), mar=c(2, 5, 2, 1)+0.1)
tmp.low <- tab.series[, paste("F",qvec[1], sep="_")]
tmp.med <- tab.series[, paste("F",qvec[2], sep="_")]
tmp.up <- tab.series[, paste("F",qvec[3], sep="_")]
matplot(years, cbind(tmp.low, tmp.up), type="n", xlab="", ylab="", cex.lab=1.5, cex.axis=1.5, main="", cex.main=2, xlim=c(min(years), max(years)+nproj), ylim=c(min(tmp.low,Fvector),max(tmp.up,Fvector)))
mtext(expression(F[y]), side=2, line=3, cex=1.5)
polygon(c(years,rev(years)),c(tmp.low, rev(tmp.up)), border=NA, col="grey")
lines(years, tmp.med, type="o",pch=19, cex=1.5, lwd=2)
for (i in 1:length(Fvector)){
  tmp <- proj.all[i, match(paste("F",1:nproj,sep=""), colnames(proj.all))]
  tmp <- c(tmp.med[n], tmp)
  matplot(years[n]+(0:nproj), tmp, add=T, type="l", lty=1, lwd=2, col=1, cex=1.5)
}
msy.low <- tab.series[, paste("Fmsy",qvec[1], sep="_")]
msy.med <- tab.series[, paste("Fmsy",qvec[2], sep="_")]
msy.up  <- tab.series[, paste("Fmsy",qvec[3], sep="_")]
abline(h=msy.low, col='red', lty=2)
abline(h=msy.med, col='red')
abline(h=msy.up, col='red', lty=2)
dev.off()

# Table of results:

write.table(proj.all, file="catch-table.txt", append=T, row.names=F, col.names=F)

# save

save.image("main.RData")


####################################################################
####################################################################
