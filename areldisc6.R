## Analysis of DHS Couple HIV Data
###################################################################### 
## Steve Bellan, September 2012
## steve.bellan@gmail.com
###################################################################### 
## The goal of this analysis is to determine probability of males &
## females being infected from within and outside their relationships.
###################################################################### 
## Set directory where output is to be written
if(getwd()=="/home/ubuntu/files")       # if on cloud
  {
    setwd("~/files/")
  }else{                                # if on sb mac
    setwd("~/Dropbox/disc mod backup/Couples Model Revision 121013/R scripts & Data Files/")
  }
library(plotrix)
library(ade4)
library(mvtnorm)
library(mnormt)
library(multicore)
library(rjags)
library(coda)
library(abind)
library(Hmisc)
rm(list=ls())
######################################################################
## data grouping indices are
######################################################################
##      [,1]  [,2]       [,3]    [,4]      [,5]     [,6]     [,7]        [,8]
## [1,] "DRC" "Ethiopia" "Kenya" "Lesotho" "Malawi" "Rwanda" "Swaziland" "WA"
##      [,9]     [,10]      [,11]             
## [1,] "Zambia" "Zimbabwe" "simulated"
######################################################################
## sample sizes
##              DRC         Ethiopia            Kenya          Lesotho 
##             1199             5347             1620             1099 
##           Malawi           Rwanda        Swaziland               WA 
##             3045             1749              431             7915 
##           Zambia         Zimbabwe          simulated 
##             1598             2825             K.sim*surv 
## west africa is pooled and analyzed in a separate script modified to
## deal with multiple countrie's prevalences simultaneously.
######################################################################
## copy below line for commands

## nohup R CMD BATCH '--args group.ind=4 short.test=F all.cores=T num.cores=8 d.nburn=300 d.nthin=3 d.niter=1000 nburn=500 nthin=1 niter=7500 survive=T tell=100 adapt=F term.on.finish=T simul=F K.sim=1000 heter=F low.coverage.arv=F partner.arv=F fsd.sens=F' areldisc6.R &
##
## First read in the arguments listed at the command line in the following way:
simpars <- c(bmb = .01, bfb = .03,      # for plotting/simulating purposes
             bme = .012, bfe=.01,
             bmp=.03, lrho=0)
args=(commandArgs(TRUE))
## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)>0)
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
## if arguments not supplied, set defaults
##  which country
aaaa <- ""               #need at least one object so ls() isn't empty
if(sum(ls() %in% "group.ind")==0)  group.ind <- 4 # default kenya
## short test on 40 couples
if(sum(ls() %in% "short.test")==0)  short.test <- F
## use multiple cores
if(sum(ls() %in% "all.cores")==0)  all.cores <- F
## how many cores to use
if(sum(ls() %in% "num.cores")==0)  num.cores <- 1
if(sum(ls() %in% "seed.bump")==0)  seed.bump <- 0 # change seed
if(sum(ls() %in% "adapt")==0)  adapt <- T# do adaptive iterations
if(sum(ls() %in% "nburn")==0)  d.nburn <-  500 # adaptive iterations
if(sum(ls() %in% "nthin")==0)  d.nthin <- 3
if(sum(ls() %in% "niter")==0)  d.niter <- 3500
if(sum(ls() %in% "nburn")==0)  nburn <-  500 # normal iterations
if(sum(ls() %in% "nthin")==0)  nthin <- 1
if(sum(ls() %in% "niter")==0)  niter <- 1500
if(sum(ls() %in% "survive")==0)  survive <- T # account for survival to sampling?
if(sum(ls() %in% "tell")==0) tell <- 200      # every tell-th iteration say how far it's gone
## terminate ec2 instance after finishing (after copying to s3 bucket)
if(sum(ls() %in% "term.on.finish")==0)  term.on.finish <- T
if(sum(ls() %in% "simul")==0)  simul <- F
if(sum(ls() %in% "K.sim")==0)  K.sim <- 1000
if(sum(ls() %in% "lrho.sd")==0) lrho.sd <- 1/2
if(sum(ls() %in% "heter")==0)  heter <- F #have individual heterogeneity in simulation
if(sum(ls() %in% "gofsim")==0)  gofsim <- 126 # how many gof reps to do (*8 for each core)
if(sum(ls() %in% "low.coverage.arv")==0)  low.coverage.arv <- F # assume 100% (F) or 50% (T) of ARV coverage reduces transm
if(sum(ls() %in% "partner.arv")==0)  partner.arv <- F # ART coverage affects within partnership transmission?
if(sum(ls() %in% "fsd.sens")==0)  fsd.sens <- F # conduct sens analys for female sexual debut?
trans.ratio <- 1
seed <- 8                               # set automatically to 1:8 if all.cores=T

if(simul)                               # if simulating
  {
    ## choose parameters to simulate if simulating
    ## load cleaned DHS data sets for all countries plus simulate one data set
    ## this also loads epicm, epicf, csurv, and sources pcalc.R
    source("sim dhs dat.R")
    load("alldhs plus sim.Rdata")
    ds <- levels(dat$ds)[grepl("simulated",levels(dat$ds))]
    group <- levels(dat$group)[group.ind]
    dat <- dat[dat$ds %in% ds,]           # only looking at that data set
    print(paste("analyzing",paste(ds, collapse = " & "),"couples data using", dat$epic.nm[1],
                "epidemic curve"))
    save(dat, file = paste(group,".Rdata")) # just to indicate which ds working on when looking in file dir
  }else{                                # otherwise load things manually
    load("alldhs.Rdata")
    ## load cleaned epidemic curves (w/ art-normalized prevalence for
    ## 15-49 by M & F
    if(low.coverage.arv)
      { # if we're assuming 50% of ARV coverage results in no transmission
        load("allepicm.5.Rdata")
        load("allepicf.5.Rdata")
      }else{ # if we're assuming 100% of ARV coverage results in no transmission
        load("allepicm.Rdata")
        load("allepicf.Rdata")
      }
    load("art.prev.Rdata")
    ## load survival curve cs, monthly survival
    load("csurv.Rdata")
    ## source probability of parameter calculator
    source("pcalc5.R")
    ## If not simulating
    wa <- levels(dat$group)[group.ind]=="WA"
    if(!wa)
      {
        group <- levels(dat$group)[group.ind]
        print(paste(nlevels(dat$ds),"country data sets"))
        ds <- levels(dat$ds)[grepl(group, levels(dat$ds))]   # data set we're working on
        dat <- dat[dat$ds %in% ds,]           # only looking at that data set
        print(paste("analyzing",paste(ds, collapse = " & "),"couples data using", dat$epic.nm[1],
                    "epidemic curve"))
        if(fsd.sens) # if doing sensitivity analysis to female sexual debut
          {
            ## Of females saying their sexual debut occured at
            ## marriage, assume 30% were lying and that actually
            ## sexual debut occurred one year earlier.
            print(paste("lowering female sexual debut by one year for 30% of females that stated they first started having sex at marriage (sensitivity analysis)"))
            fsd.mar.ind <- which(dat$tfs==dat$tmar) #which sexual debuts occurred at marriage
            alter.ind <- sample(fsd.mar.ind, size = round(length(fsd.mar.ind)*.3)) # ones to alter
            dat$tfs[alter.ind] <- dat$tfs[alter.ind] - 12            
          }
        save(dat, file = paste(group,".Rdata")) # just to indicate which ds working on when looking in file dir
      }else{
        group <- levels(dat$group)[group.ind]
        print(paste(nlevels(dat$ds),"country data sets"))
        dat <- dat[dat$group == group,]                 # west africa
        print("analyzing pooled west african couples data using their respective epidemic curves")
        save(dat, file = "wa.Rdata") # just to indicate which ds working on when looking in file dir
      }
  }

## Sets chain length etc, have options to do a short run or one core.
if(short.test)
  {
    nburn <- 200
    niter <- 400
    nthin <- 1
    nc <- ifelse(all.cores, num.cores, 1)        # number of MCMC chains
  }else{
    nburn <- nburn
    niter <- niter
    nthin <- nthin
    nc <- ifelse(all.cores, num.cores, 1)        # number of MCMC chains
  }

## Get before couple duration bd, where bd = max(mbd,fbd)
dat$bd <- apply(cbind(dat$tmar-dat$tms,dat$tmar-dat$tfs), 1, max)
## Get couple duration cd
dat$cd <- dat$tint - dat$tmar
K <- nrow(dat)
testpars <- simpars                     # for simulation runs

sd.props <- c(bmb.sd = .001, bfb.sd = .004, # sd of proposal distr before adaptive phase
              bme.sd = .0015, bfe.sd = .0015,
              bmp.sd = .003, lrho.sd = .15)
parnames <- c("bmb","bfb","bme","bfe","bmp","lrho")
sd.name <- paste(parnames, ".sd", sep = "") #proposal distr sd names

start.time1 <- Sys.time()
sigma.found <- "sigmas.Rdata" %in% list.files()
if(sigma.found)
  {
    print("loading covar matrix for multinormal sampling from file")
    load("sigmas.Rdata")            
    sigma <- sigmas[,,group.ind]        #  choose covar matrix for this data set (from earlier adaptations)
  }else{
    if(adapt)   print("no covar matrix found in file folder, starting with block univariate sampling")
    sigma <- NA
  }

if(all.cores)
  {
    wrp <- function(seed=1, multiv=F, covar=NULL,
                    niter, survive, browse,
                    nthin, nburn)
      {
        inits.temp <- init.fxn(seed = seed)
        sampler(sd.props = sd.props, inits = inits.temp, browse = browse,
                multiv = multiv, covar = covar, lrho.sd = lrho.sd,
                verbose = T, tell = tell, seed = seed,
                niter = niter, survive = survive,
                nthin = nthin,
                nburn = nburn)
      }
    if(adapt)        ## Adaptive (d.) phase to get multivariate normal sampler
      {
        print("beginning adaptive phase")
        d.out <- mclapply((seed.bump + 1:nc), wrp, 
                          multiv = sigma.found, covar = sigma, browse=F,
                          survive = survive, niter = d.niter, nthin = d.nthin, nburn = d.nburn)
        save.image(file="workspace.Rdata")
        ## reformat into mcmc object
        mcmc.d.out <- list(NA)
        d.aratio <- 0
        for(ii in 1:nc)
          {
            mcmc.d.out[[ii]] <- as.mcmc(t(d.out[[ii]][[1]]))
            d.aratio <- d.aratio + d.out[[ii]]$aratio
            if(ii==1)
              {
                init.adapt <- d.out[[ii]]$inits
              }else{
                init.adapt <- rbind(init.adapt, d.out[[ii]]$inits)
              }
          }
        mcmc.d.out <- mcmc.list(mcmc.d.out)
        d.aratio <- d.aratio/nc
        print(paste("adaptive phase aratio is", round(d.aratio,2), ".Rdata"))
        bmb.vec <- unlist(mcmc.d.out[,"bmb"])
        bfb.vec <- unlist(mcmc.d.out[,"bfb"])
        bmp.vec <- unlist(mcmc.d.out[,"bmp"])
        lrho.vec <- unlist(mcmc.d.out[,"lrho"])
        bme.vec <- unlist(mcmc.d.out[,"bme"])
        bfe.vec <- unlist(mcmc.d.out[,"bfe"])
        posts <- data.frame(bmb = bmb.vec, bfb = bfb.vec, bme = bme.vec, bfe = bfe.vec, bmp = bmp.vec, lrho = lrho.vec)
        sbpairs(posts, truepars = testpars, show.lines = simul,
                file.nm = "posterior pairs after adaptive phase",
                width = 12, height = 12,
                cex = 1, col = "black", nrpoints = 200, do.jpeg = T)
        mu <- mean(posts)
        sigma <- cov.wt(posts)$cov #estimate covariance matrix, then plot what proposal distr is gonna look like
        sbpairs(rmnorm(3000, mean = mu, varcov = sigma), truepars = testpars,
                 show.lines = simul,
                file.nm = "adapted proposal distr after adaptive phase",
                width = 12, height = 12,
                cex = 1, col = "black", nrpoints = 200, do.jpeg = T)
      } # end adaptive phase, sigma has been updated if this is run, otherwise it's been loaded
    print("beginning sampling")
    ## use multicore to run a chain on each core with seeds 1:nc (where nc is number chains/core)
    out <- mclapply((seed.bump + 1:nc), wrp,  multiv = T, covar = sigma, browse = F,
                    survive = survive, niter = niter, nthin = nthin, nburn = nburn) 
    save.image(file="workspace.Rdata")
    ## reformat into mcmc object
    mcmc.out <- list(NA)
    aratio <- 0
    for(ii in 1:nc)
      {
        mcmc.out[[ii]] <- as.mcmc(t(out[[ii]][[1]]))
        aratio <- aratio + out[[ii]]$aratio
        if(ii==1)
          {
            init.samp <- out[[ii]]$inits
          }else{
            init.samp <- rbind(init.samp, out[[ii]]$inits)
          }
      }
    mcmc.out <- mcmc.list(mcmc.out)
    aratio <- aratio/nc
    bmb.vec <- unlist(mcmc.out[,"bmb"])
    bfb.vec <- unlist(mcmc.out[,"bfb"])
    bmp.vec <- unlist(mcmc.out[,"bmp"])
    lrho.vec <- unlist(mcmc.out[,"lrho"])
    bme.vec <- unlist(mcmc.out[,"bme"])
    bfe.vec <- unlist(mcmc.out[,"bfe"])
    posts <- data.frame(bmb = bmb.vec, bfb = bfb.vec, bme = bme.vec, bfe = bfe.vec, bmp = bmp.vec, lrho = lrho.vec)
    sigma.e <- cov.wt(posts)$cov #estimate covariance matrix, then plot what proposal distr is gonna look like
    save(sigma.e, file = "sigma.Rdata")
    sbpairs(posts, truepars = testpars, show.lines = simul,
            file.nm = "pairs after sample phase",
            width = 12, height = 12,
            cex = 1, col = "black", nrpoints = 200, do.jpeg = T)
  }else{                                # if just doing one core
    inits <- init.fxn(seed = seed)
    out <- sampler(sd.props = sd.props, inits = inits,
            verbose = T, tell = 20, seed = seed, lrho.sd = lrho.sd,
            niter = niter, survive = survive,
            nthin = nthin,
            nburn = nburn)
    mcmc.out <- as.mcmc(t(out[[1]]))
    aratio <- out$aratio
  }  
save.image(file="workspace.Rdata")
save(aratio, file = paste("aratio is", round(aratio,2), ".Rdata"))
print(paste("aratio is", round(aratio,2), ".Rdata"))

###################################################################### 
## Extract parameters' posterior and save to file (use parallel processing)
sum.wrap <- function(col.ind, xx)
  {
    summary(xx[,col.ind], quantiles = c(.025, .5, .975))$quant
  }
pars <- abind(mclapply(1:ncol(mcmc.out[[1]]), sum.wrap,  xx = mcmc.out), along = 0, new.names=colnames(mcmc.out[[1]]))
## never save over old files
file.name <- paste("pars", format(Sys.time(), "%Y%m%d"), sep = "-")
ii <- 1
while(file.exists(paste(file.name, ".Rdata",sep="")))
  {
    file.name <- paste(file.name, "-",ii)
    ii <- ii+1
  }
pars.name <- paste(file.name, ".Rdata",sep="")
save(pars, file = pars.name)
chain.name <- paste(file.name, "chains.Rdata",sep="")
save(mcmc.out, file = chain.name)
out.csv <- signif(pars,3)
out.csv <- data.frame(median = out.csv[,2], CI95 = paste("(",out.csv[,1],", ",out.csv[,3],")", sep=""))
write.csv(out.csv, file="out.csv")
######################################################################
save.image(file="workspace.Rdata")

######################################################################
###################################################################### 
## Figure 0 - MCMC diagnostics (for beta's only)
###################################################################### 
pdf("Fig 0 - mcmc diagnostics.pdf")
show.cols <- colnames(mcmc.out[[1]]) %in% c(parnames,"bfp")
plot(mcmc.out[,show.cols])
dev.off()
######################################################################

## Gelman-Rubin diagnotics (save to file)
######################################################################
library(coda)
betpars <- c("bmb", "bfb", "bme", "bfe", "bmp", "bfp")
beta.out <- mcmc.out[, betpars, drop=FALSE]
gelout <- gelman.diag(beta.out)
gelout

###################################################################### 
## Figure 1 - plot aa and bb by mysa and fysa and gg and dd by mardur
###################################################################### 
pdf("Fig 1 - abgd by yr.pdf", width = 7, height = 8)
par(mfrow=c(3,2))
nn <- nrow(dat)
mysa.vec <- 1:40 
fysa.vec <- 1:40 
mardur.vec <- 1:40 
###################################################################### 
xlim <- c(0, 46)
xlab <- "years sexually active before relationship"
plot(0,0, xlim = xlim, ylim = c(0,1),
     type = "n", bty = "n",
     ylab = expression(1-exp(-beta[M]*(x[M,i]-r[i]))), xlab = xlab,
     main = "male risk of infection prior to relationship")
polygon(c(mysa.vec, rev(mysa.vec)),
        c(1-exp(-mysa.vec*pars["bmb","2.5%"]), rev(1-exp(-mysa.vec*pars["bmb","97.5%"]))),
        col = "gray", border = NA)
lines(mysa.vec, 1-exp(-mysa.vec*pars["bmb","50%"]), lwd = 2)
plot(0,0, xlim = xlim, ylim = c(0,1),
     type = "n", bty = "n",
     ylab = expression(1-exp(-beta[F]*(x[F,i]-r[i]))), xlab = xlab,
     main = "female risk of infection prior to mariage")
polygon(c(mysa.vec, rev(mysa.vec)),
        c(1-exp(-mysa.vec*pars["bfb","2.5%"]), rev(1-exp(-mysa.vec*pars["bfb","97.5%"]))),
        col = "gray", border = NA)
lines(mysa.vec, 1-exp(-mysa.vec*pars["bfb","50%"]), lwd = 2)
###################################################################### 
xlab <- "relationship duration"
plot(0,0, xlim = xlim, ylim = c(0,1),
     type = "n", bty = "n",
     ylab = expression(1-exp(-beta[MR]*r[i])), xlab = xlab,
     main = "male risk of infection from extra-couple sex")
polygon(c(mysa.vec, rev(mysa.vec)),
        c(1-exp(-mysa.vec*pars["bme","2.5%"]), rev(1-exp(-mysa.vec*pars["bme","97.5%"]))),
        col = "gray", border = NA)
lines(mysa.vec, 1-exp(-mysa.vec*pars["bme","50%"]), lwd = 2)
plot(0,0, xlim = xlim, ylim = c(0,1),
     type = "n", bty = "n",
     ylab = expression(1-exp(-beta[FR]*r[i])), xlab = xlab,
     main = "female risk of infection from extra-couple sex")     
polygon(c(mysa.vec, rev(mysa.vec)),
        c(1-exp(-mysa.vec*pars["bfe","2.5%"]), rev(1-exp(-mysa.vec*pars["bfe","97.5%"]))),
        col = "gray", border = NA)
lines(mysa.vec, 1-exp(-mysa.vec*pars["bfe","50%"]), lwd = 2)
######################################################################      
xlim <- c(0, 45)
plot(0,0, xlim = xlim, ylim = c(0,1),
     type = "n", bty = "n",
     ylab = expression(1-exp(beta[MF]*r[i])), xlab = xlab,
     main = "prob of m->f transmission in couple")
polygon(c(mardur.vec, rev(mardur.vec)),
        c(1-exp(-mardur.vec*pars["bfp","2.5%"]), rev(1-exp(-mardur.vec*pars["bfp","97.5%"]))),
        col = "gray", border = NA)
lines(mardur.vec, 1-exp(-mardur.vec*pars["bfp","50%"]), lwd = 2)
plot(0,0, xlim = xlim, ylim = c(0,1),
     type = "n", bty = "n",
     ylab = expression(1-exp(beta[MF]*r[i])), xlab = xlab,
     main = "prob of f->m transmission in couple")
polygon(c(mardur.vec, rev(mardur.vec)),
        c(1-exp(-mardur.vec*pars["bmp","2.5%"]), rev(1-exp(-mardur.vec*pars["bmp","97.5%"]))),
        col = "gray", border = NA)
lines(mardur.vec, 1-exp(-mardur.vec*pars["bmp","50%"]), lwd = 2)
dev.off()
######################################################################



###################################################################### 
## Figure 2 - prob an infection is from extracouple sex by serostatus
## by ysa and reldur
###################################################################### 
cex <- .8
medpars <- pars[parnames,2]
pis <- pcalc(medpars, dat = dat, trace = T, give.pis=T, survive=survive, lrho.sd = lrho.sd)$pis
######################################################################
breaks <- seq(0,1, by = .1)
xlim <- c(0, 35)
ylim <- c(0, 35)
xlab <- "years sexually active before relationship"
ylab <- "relationship duration"
mains <- c("M+F+","M+F-","M-F+","M-F-")

######################################################################
cex <- .65
rmp <- colorRamp(c("yellow","red"))     #create color ramp
pal <- colorRampPalette(c("yellow","red"))     #create color palette (for legend)
## pdf("Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.pdf",
##     width = 6.5, height = 4
tiff <- F
if(tiff)
  {
    tiff("Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.tiff",
         width = 5.5, height = 4, units = "in", res = 300)
  }else{
    pdf("Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.pdf",
         width = 5.5, height = 4)
  }
par(mar = c(.5,1,1,.5), oma = c(5,5,0,9))
                                        # do it for each data set since the interview times were different for
                                        # same countries and so the plot shold show how long the couples were
                                        # together too
for(dd in unique(dat$epic.nm)) 
  {
    cc <- dat$epic.ind[dat$epic.nm==dd][1]   #find epidemic curve for that data set
    epic.col <- "blue"
    layout(t(matrix(1:4,2,2)))
    ylim <- c(0, 30) #max(dat$m.bef.pm, dat$f.bef.pm))/12
    xlim <- c(1975,2012)
    ylab <-  "YSA before couple formation"
    xlab <- "date of couple formation"
    mains <- c("M+F+","M+F-","M-F+","M-F-")
    mains <- c("B","A","","C")
    for(ii in 2:1)
      {
        if(ii!=4) cols.show <- rgb(rmp(pis$piCe.A[dat$ser==ii & dat$epic.ind==cc]), max = 255)
        if(ii==4) cols.show <- "black"
        plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii & dat$epic.ind==cc],
             1/12*(dat$tmar-dat$tms)[dat$ser==ii & dat$epic.ind==cc],
             col = cols.show, xlab = "",
             pch = 19, cex = cex, axes = F,
             xlim = xlim, ylim = ylim, bty = "n",
             main = mains[ii])
        xs <- epicm[,1]/12 + 1900
        if(ii==2) axis(2, at = seq(0, 30, by = 10), las = 2)
        if(ii==1) axis(2, at = seq(0, 30, by = 10), labels = NA, las = 2)
        axis(1, at = seq(1980, 2010, by = 10), labels = NA)
        lines(xs[xs>1975], epicf[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
        if(ii==1) axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
      }
    mains <- c("F+M+","","F+M-","F-M-")
    mains <- c("D","","C","F")    
    for(ii in c(3,1))
      {
        if(ii!=4) cols.show <- rgb(rmp(pis$piC.eA[dat$ser==ii & dat$epic.ind==cc]), max = 255)
        if(ii==4) cols.show <- "black"
        plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii & dat$epic.ind==cc],
             1/12*(dat$tmar-dat$tfs)[dat$ser==ii & dat$epic.ind==cc],
             ##          1/12*(dat$f.bef.pm)[dat$ser==ii],
             col = cols.show, pch = 19, las = 2, axes = F,
             xlim = xlim, ylim = ylim, bty = "n", xlab = "", cex = cex,
             main = mains[ii])
        if(ii==3) axis(2, at = seq(0, 30, by = 10), las = 2)
        if(ii==1) axis(2, at = seq(0, 30, by = 10), labels = NA, las = 2)
        axis(1, at = seq(1980, 2010, by = 10), las = 2)
        lines(xs[xs>1975], epicm[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
        if(ii==1) axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
      }
    cex.ax <- .7
##     mtext("male YSA before couple formation", side = 2, outer = T, line = 2, cex = cex.ax, adj = .92)
##     mtext("female YSA before couple formation", side = 2, outer = T, line = 2, cex = cex.ax, adj = .08)
##     mtext(paste(xlab), side = 1, outer = T, line = 2, cex = cex.ax)
##     mtext("F HIV pop. prevalence", side = 4, outer = T , line = 2, adj = .9, cex = cex.ax)
##     mtext("M HIV pop. prevalence", side = 4, outer = T , line = 2, adj = .15, cex = cex.ax)
##     mtext(dd, side = 3, outer = T, line = .5) # data set title
    ##       mtext(paste("male",xlab), side = 1, outer = T, line = -21, cex = 1.2)
  } # end loop through pooled data sets
dev.off()
######################################################################

library(plotrix)
###################################################################### 
tiff("fig 2 legend.tiff", width = .9, height = 3, units = "in", res = 300)
par(mar=rep(0,4))
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9))
cols <- pal(100)
color.legend(0,.1,.05,.8, seq(0,1, length=11), rect.col = cols, gradient = "y", cex = .8)
dev.off()
######################################################################


###################################################################### 
pdf("fig 2 legend.pdf", width = .9, height = 3)
par(mar=rep(0,4))
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9))
cols <- pal(100)
color.legend(0,.1,.05,.8, seq(0,1, length=11), rect.col = cols, gradient = "y", cex = .8)
dev.off()
######################################################################


###################################################################### 
## Figure 3 - Hazard posteriors
###################################################################### 
pdf("Fig 3 - hazard posteriors.pdf", width = 6, height = 3.5)
par(mar=c(5,10,0.5,.5))
show <- match(c("bmb","bfb",
                "bme","bfe",
                "bmp","bfp"), rownames(pars))
labs <- c("M before relationship", "F before relationship",
          "M extra-couple sex", "F extra-couple sex",
          "M from partner", "F from partner")
plot(12*pars[show,2], 6:1, pch = 15, cex = 1,
     ylim = c(.8,6.2), xlim=c(0,12*max(pars[show,])), bty = "n", yaxt = "n",
     ylab="", xlab = expression(beta[yearly]))
arrows(12*pars[show,1],6:1,  12*pars[show,3],  6:1, angle = 90, length=.1, code = 3, lwd = 2)
axis(2, at = 6:1, label = labs, las = 2, cex=2)
## axis(1, at = seq(0,.3, by = .05))
dev.off()
###################################################################### 


###################################################################### 
## Figure 6 - bfp/bmp posterior & prior
###################################################################### 
pdf("Figure 6 - lmf to lfm post.pdf", width = 6, height = 5)
par(mar = c(4.5,4,2.5,0))
xx <- seq(0, 10, length.out=1000)
yy <- dlnorm(xx, mean = log(trans.ratio), sd = 1/2)
hist(exp(unlist(mcmc.out[,"lrho"])), breaks = 100, col="gray", border = NA,
     xlab = expression(beta[mf]/beta[fm]), ylab = "probability density", main = "",
     xlim = c(0,10), freq = F, ylim = c(0, 1))
lines(xx, yy, lwd = 2)
legend("topright", leg = c("prior","posterior"), col = c("black","gray"), lwd = 2, bty = "n")
dev.off()

###################################################################### 
## Figure 7 - bfp/bmp prior
###################################################################### 
pdf("Figure 7 - prior beta out.pdf", width =4, height = 3.5)
par(mar = c(4.5,4,2.5,0))
xx <- seq(0, 6, length.out=1000)
yy <- dlnorm(xx, mean = log(trans.ratio), sd = 1/2)
plot(xx,yy, lwd = 2, type = "l",
     xlab = expression(beta[Fpartner]/beta[Mpartner]), ylab = "probability density", main = "",
     xlim = c(0,6), ylim = c(0, 1), bty = "n")
dev.off()
######################################################################


###################################################################### 
## Figure 8 - years sexally active before relationship distribution
###################################################################### 
pdf("Figure 8 - ysa for disc couples.pdf", width = 6, height = 8)
par(mfrow=c(2,1))
hist(1/12*(dat$tmar-dat$tms)[dat$ser==2], breaks = seq(0,55,by=1/2), xlab = "MYSA before relationship",
     main="M+ F- couples", col = "black")
hist(1/12*(dat$tmar-dat$tfs)[dat$ser==3], breaks = seq(0,55,by=1/2), xlab = "FYSA before relationship",
     main="M- F+ couples", col = "black")
dev.off()
###################################################################### 

###################################################################### 
## Figure 9 - beta ratios
###################################################################### 
pdf("Figure 9 - beta ratios.pdf", width = 5, height = 3.5)
par(mar=c(4,7,.5,.5))
show <- match(c("rr.mf.bef", "rr.mf.exc",
                "rr.m.out", "rr.f.out"), rownames(pars)) # SHOWING THE INVERSES (SEE LABELS)
labs <- c(expression(beta[Fbefore] / beta[Mbefore]),
          expression(beta[Fduring] / beta[Mduring]),
          expression(beta[Mbefore] / beta[Mduring]),
          expression(beta[Fbefore] / beta[Fduring]))
plot(1/pars[show,2], 4:1, pch = 15, cex = 2, log="x",
     ylim = c(.8,4.2), xlim=c(.05,20), bty = "n", axes=F, ylab="", xlab = "rate ratio")
segments(1, 4.2, 1, .8, lty = 2, lwd = 2)
arrows(1/pars[show,1],4:1,  1/pars[show,3],  4:1, angle = 90, length=.1, code = 3, lwd = 2)
axis(2, at = 4:1, label = labs, las = 2, cex=2)
axis(1, at = c(.05,.1,.2,.5,1,2,5,10,20), las = 2,
     label = c("1/20","1/10","1/5","1/2","1","2","5","10","20"))
dev.off()
######################################################################

###################################################################### 
## Figure 9b - beta ratios
###################################################################### 
pdf("Figure 9b - contact ratios.pdf", width = 5, height = 3.5)
par(mar=c(4,7,.5,.5))
show <- match(c("rr.mf.bef.cont", "rr.mf.exc.cont",
                "rr.m.out", "rr.f.out"), rownames(pars)) # SHOWING THE INVERSES (SEE LABELS)
labs <- c(expression(beta[Fbefore] / beta[Mbefore]),
          expression(beta[Fduring] / beta[Mduring]),
          expression(beta[Mbefore] / beta[Mduring]),
          expression(beta[Fbefore] / beta[Fduring]))
plot(1/pars[show,2], 4:1, pch = 15, cex = 2, log="x",
     ylim = c(.8,4.2), xlim=c(.05,20), bty = "n", axes=F, ylab="", xlab = "rate ratio")
segments(1, 4.2, 1, .8, lty = 2, lwd = 2)
arrows(1/pars[show,1],4:1,  1/pars[show,3],  4:1, angle = 90, length=.1, code = 3, lwd = 2)
axis(2, at = 4:1, label = labs, las = 2, cex=2)
axis(1, at = c(.05,.1,.2,.5,1,2,5,10,20), las = 2,
     label = c("1/20","1/10","1/5","1/2","1","2","5","10","20"))
dev.off()
######################################################################

###################################################################### 
pdf("Figure 11 - # next year of transmission.pdf", width = 4, height = 5)
  ## Calculate probability each uninfected person is infected in next
  ## year using estimated beta's as well as current prevalence.
  par(mar=c(10,4,0.5,.5))
show <- match(c("n.m.part.tot", "n.m.exc.tot", "n.f.part.tot", "n.f.exc.tot"),
              rownames(pars))                
nms <- c("male: partner", "male: extracouple", "female: partner", "female: extracouple")
bp <- barplot(pars[show,2], names.arg = nms, col = c("blue","red","blue","red"),
ylab = "# incident infections in next 12 months",
las = 2, ylim = c(0, max(pars[show,])))
arrows(bp,pars[show,1],bp,pars[show,3], angle = 90, length=.1, code = 3, lwd = 2)
mtext(paste("N =",c(sum(dat$ser %in% c(3:4)))), side = 1, adj = .22, line = 0)
mtext(paste("N =",c(sum(dat$ser %in% c(2,4)))), side = 1, adj = .82, line = 0)
dev.off()
######################################################################


###################################################################### 
pdf("Figure 12 - prob any new inf is extracouple.pdf", width = 3, height = 5)
## Calculate probability each uninfected person is infected in next
## year using estimated beta's as well as current prevalence.
par(mar=c(3,4,0.5,.5))
nms <- c("male","female")
show <- match(c("prop.exc.m", "prop.exc.f"), rownames(pars))
bp <- barplot(pars[show,2], names.arg = nms,
        ylab = "probability new infection is from extracouple intercourse",
        las = 1, ylim = c(0, 1))
arrows(bp,pars[show,1],bp,pars[show,3], angle = 90, length=.1, code = 3, lwd = 2)
dev.off()
######################################################################


###################################################################### 
pdf("Fig 16 - survival times.pdf", width = 4.5, height = 3)
par(mar=c(4,4,1,1))
## create discrete cumulative probability of mortality
## it's age dependent, fit by eyeballing to CASCADE study
xseq <- 1:(12*60)
age.seq <- seq(20*12,60*12, by = 10*12)
cols <- c("orange","green","pink","light blue", "dark gray")
for(ii in 1:length(age.seq))
  {
    aa <- age.seq[ii]
    shp <- 2.3
   scl <- 2000/shp/(aa/12)^.53
    cmort <- pweibull(xseq, shape = shp, scale = scl)
    csurv.temp <- 1-cmort
    if(ii==1) plot(xseq/12, csurv.temp, type = "l", xlim = c(0,25), col=cols[ii], lwd = 3, yaxt = "n",
         xlab = "years since seroconversion", ylab="probability of survival", bty = "n")
    if(ii!=1) lines(xseq/12, csurv.temp, type = "l", col=cols[ii], lwd = 3)
  }
legend("topright",paste(age.seq/12,"yrs old"), col = cols, pch = 15, bty = "n", cex = 1, title = "age at seroconversion")
axis(2, seq(0,1,l=5), las = 2)
dev.off()
######################################################################


###################################################################### 
ctraj(medpars, dat, browse =F,
      plot.cpls = sample(1:nrow(dat),10,replace=F),
      surv = survive,                   # show surv curves o plot
      nsurv = F,                        # don't show marginal curves 
      lty.surv = 1,                     # line type for surv curves
      dead =T, col.dead = "black", lty.dead = 1,
      survive = survive,                # plot point at survival curve
      pdf.name = "Figure 14 - Probability trajectories for some couples.pdf")
###################################################################### 

######################################################################
## Model fit - simulate serostatus outcomes based on random draws from
## posterior parameter distributions, calculate p(simulated data|pars)
## for lots of simulations, compare this distribution to p(observed
## data|pars)

pdf("Figure 15 - model fit.pdf")
calgof <- function(seed = 1, nsim)
  {
    randpars <- posts[sample(1:nrow(posts), nsim),]
    llreal <- rep(NA, nsim)
    llsim <-  rep(NA, nsim)
    for(ii in 1:nsim)
      {
        if(ii %% 5 == 0) print(paste("working on sim", ii, "of", nsim))
        if(survive)
          {
            fser <- pcalc(randpars[ii,], dat = dat, trace = T, give.pis=T, survive = survive, lrho.sd = lrho.sd)$pser.a
          }else{
            fser <- pcalc(randpars[ii,], dat = dat, trace = T, give.pis=T, survive = survive, lrho.sd = lrho.sd)$pser
          }    
        fser <- fser / as.matrix(rowSums(fser))[,rep(1,4)] # normalize so sums to 1
        colnames(fser) <- 1:4
        temp.sers <- as.numeric(as.vector(rMultinom(fser, 1)))     #
        ## Likelihood of sampled parameters | sim data
        llsim[ii] <- sum(log(fser[cbind(1:nrow(fser), temp.sers)]))
        ## Likelihood of sampled parameters | real data
        llreal[ii] <- sum(log(fser[cbind(1:nrow(fser), dat$ser)]))
      }
    out <- data.frame(llsim, llreal)
    return(out)
  }
calgof.out <- mclapply(1:8, calgof, nsim = gofsim)
for(ii in 1:8)
  {
    if(ii==1)
      {
        cgout <- calgof.out[[ii]]
      }else{
        cgout <- rbind(cgout, calgof.out[[ii]])
      }
  }
attach(cgout)
pdata <- mean(llreal)
nsim <- length(llsim)
hist(llsim, col = "black", xlim = c(min(llsim, pdata), max(llsim,pdata)),
     xlab = "log probability of simulated data given median parameters",
     ylab = "frequency under multinomial sampling")
segments(pdata,0,pdata, nsim/6, col = "red")
text(pdata,nsim/6, "observed data", pos=4, col = "red")
dev.off()
ecdf.lp <- ecdf(llsim)
p.val <- ecdf.lp(pdata)
save(p.val, file = paste(signif(p.val,4),"prob of data fitting this poorly to model.Rdata"))

## Show bias in estimates for simulated data
if(simul)
  {
    ## divide simulated data by routes of transmission for males
    pibUA.lg <- dat$cat.nm %in% c("mb.A", "hbeA", "hbpA", "hb1b2A", "hb2b1A")
    pieUA.lg <- dat$cat.nm %in% c("me.A", "hepA", "hebA", "he1e2A", "he2e1A")
    pipUA.lg <- dat$cat.nm %in% c("hpeA", "hpbA")
    ## 
    inf.males <- sum(pibUA.lg) + sum(pieUA.lg) + sum(pipUA.lg)
    pibUA.r <- sum(pibUA.lg) / inf.males
    pieUA.r <- sum(pieUA.lg) / inf.males
    pipUA.r <- sum(pipUA.lg) / inf.males
    ## divide simulated data by routes of transmission for females
    piUbA.lg <- dat$cat.nm %in% c("f.bA", "hebA", "hpbA", "hb1b2A", "hb2b1A")
    piUeA.lg <- dat$cat.nm %in% c("f.eA", "hpeA", "hbeA", "he1e2A", "he2e1A")
    piUpA.lg <- dat$cat.nm %in% c("hepA", "hbpA")
    ## 
    inf.females <- sum(piUbA.lg) + sum(piUeA.lg) + sum(piUpA.lg)
    piUbA.r <- sum(piUbA.lg) / inf.females
    piUeA.r <- sum(piUeA.lg) / inf.females
    piUpA.r <- sum(piUpA.lg) / inf.females
    ## Create a data frame that compares true routes of transmission
    ## for observed individuals to estimates
    bdis <- c(pibUA.r, pieUA.r, pipUA.r, piUbA.r, piUeA.r, piUpA.r)
    names(bdis) <- c("pibUA","pieUA","pipUA", "piUbA","piUeA","piUpA")
    truepars <- c(simpars[1:5], simpars[5]*exp(simpars[6]))
    names(truepars)[6] <- "bfp"
    show <- c(rownames(pars)[c(1:5,7)], "pibUA","pieUA","pipUA", "piUbA","piUeA","piUpA")
    propbias <- signif(data.frame(pars[show,], c(truepars,bdis)),3)
    colnames(propbias) <- c("2.5%",    "50%",     "97.5%",   "trueval")
    if(heter)
      {
        scalars <- c(rep(sqrt(exp(1)),4), rep(1,8))
      }else{
        scalars <- rep(1,12)
      }
    biass <- (propbias[,2,] - scalars*propbias[,4,]) / propbias[,2,]
    pdf("trans distr bias.pdf", w = 8, h = 5)
    layout(matrix(c(1,2,7,10,3,4,8,11,5,6,9,12),4,3))
    par(mar = c(3,.5,2,.5))
    for(jjj in 1:nrow(propbias))
      {
        if(jjj %in% 1:2)            xlim <- c(0, .1)
        if(jjj %in% 3:4)            xlim <- c(0, .02)
        if(jjj %in% 5:6)            xlim <- c(0, .05)
        if(jjj < 7)     xlim <- c(0,.05)
        if(jjj > 6)     xlim <- c(0,1)
        plot(0,0, type = "n", xlim = xlim, ylim = c(.5, 2), bty = "n", ylab = "", yaxt = "n",
             main = paste(rownames(propbias)[jjj], ": bias =", signif(biass[jjj],2)))
        points(scalars[jjj]*propbias[jjj,4], 1:1, pch = 19, col = "red", cex = 1.5)
        arrows(propbias[jjj,1], 1:1, propbias[jjj,3], 1:1, code = 3, angle = 90, length = .03)
        points(propbias[jjj,2], 1:1, pch = 19, cex = .8)
        if(jjj < 6) points(init.samp[,jjj], rep(.8, nrow(init.samp)), col = "blue", cex = .5, pch = 19)
        if(jjj==6) points(exp(init.samp[,6])*init.samp[,5], rep(.8, nrow(init.samp)), col = "blue", cex = .5, pch = 19)
      }
    legend("topright", c("true", "estimates (95% CIs)"), col = c("red","black"), pch = 19, bty = "n")
    dev.off()
  }

if(simul)       print(paste("lprob at testpars = ", pcalc(testpars, dat = dat, trace = T, sim = F, survive = survive, give.pis=F, lrho.sd = lrho.sd)$lprob))
print(paste("lprob at medpars = ", pcalc(medpars, dat = dat, trace = T, sim = F, survive = survive, give.pis=F, lrho.sd = lrho.sd)$lprob))

######################################################################
## Print processing time to file
hours <- round(as.numeric(difftime(Sys.time(), start.time1, unit = "hour")),3)
save(hours, file = paste("took",hours,"hrs.Rdata"))
save.image(file="workspace.Rdata")      #resave workspace

print(getwd())
    
if(getwd()=="/home/ubuntu/files")       # if on cloud
  {
    ## Upload all output to S3 Bucket

    group <- sub(" ", "-", group)
    if(survive) dirnm <- paste(group, "-Deaths-", format(Sys.time(), "%Y%m%d-%H:%M:%S"), sep = "")
    if(!survive) dirnm <- paste(group, "-NoDeaths-", format(Sys.time(), "%Y%m%d-%H:%M:%S"), sep = "")
    if(simul) dirnm <- paste(dirnm, "-sim", sep = "")
    if(heter) dirnm <- paste(dirnm, "-het", sep = "")
    if(partner.arv) dirnm <- paste(dirnm, "-parv", sep = "")
    if(low.coverage.arv) dirnm <- paste(dirnm, "-lcarv", sep = "")
    if(fsd.sens) dirnm <- paste(dirnm, "-fsd.sens", sep = "")    
    dirnm <- paste(dirnm, "-rho",lrho.sd, sep = "")
    print(paste("sending to S3:", dirnm))
    comnd <- paste("s3cmd put /home/ubuntu/files s3://disc-output/",dirnm,"/ --recursive", sep ="")
    system(comnd)
    if(term.on.finish)
      {
        ## Shutdown instance in 15 minutes, give time for copying and for email alert
        system("sudo shutdown -h 15")
      }
  }
  
