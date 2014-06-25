## load cleaned DHS data sets for all countries
load("alldhs.Rdata")                    
## load cleaned epidemic curves (w/ art-normalized prevalence for
## 15-49 by M & F
load("allepicm.Rdata")
load("allepicf.Rdata")
## load survival curve cs, monthly survival
load("csurv.Rdata")
source("pcalc3.R")

if(!"group.ind" %in% ls())        group.ind <- 4
if(!"survive" %in% ls())        survive <- T
## Use this ds' couple & sex duration distributions in simulation.
group <- levels(dat$group)[group.ind]
ds <- levels(dat$ds)[grepl(group, levels(dat$ds))]   # data set we're working on
dat <- dat[dat$ds %in% ds,]           # only looking at that data set
print(paste("real",group,"serostatus"))
print(table(dat$ser))

if(length(K.sim)==0) K.sim <- 1000
odat <- dat[sample(1:nrow(dat), K.sim, replace = T),]
odat$ds <- paste(group, "simulated")
odat$group <- odat$ds
odat$ser <- NA
## Get before couple duration bd, where bd = max(mbd,fbd)+1
odat$bd <- apply(cbind(odat$tmar-odat$tms,odat$tmar-odat$tfs), 1, max) + 1
## Get couple duration cd
odat$cd <- odat$tint - odat$tmar
parnames <- c("bmb","bfb","bme","bfe","bmp","lrho")
simpars <- c(0.03250257, 0.06735142, 0.01688723, 0.01882287, 0.03143116, 0.18118647)
names(simpars) <- parnames
simpars


its <- 100                              # iterations
out <- array(NA, c(5, 4, its))          # array to store results
rownames(out) <- c("true","v1","v2","v3","v4") 
colnames(out) <- c("piGb1.", "piGe1.", "piG.b1", "piG.e1")
for(ii in 1:its)
  {
    if(ii %% 5 == 0) print(paste(ii,"of",its))
    simout <- pcalc(simpars, dat = odat, trace = F, give.pis=T,
                    survive = survive, sim = T, browse=F)
    simdat <- simout$dat
    simdat$ser <- NA
    simdat$ser[simdat$cat.nm == "s..A"] <- 4
    simdat$ser[simdat$cat.nm %in% c("mb.D", "me.D", "mb.A", "me.A")] <- 2
    simdat$ser[simdat$cat.nm %in% c("f.bD", "f.eD", "f.bA", "f.eA")] <- 3
    simdat$ser[simdat$cat.nm %in% c("hb1b2D", "hb2b1D", "hbeD","hebD",   "hbpD",   "hpbD",   "hepD",   "hpeD",   "he1e2D", "he2e1D",
                                    "hb1b2A", "hb2b1A", "hbeA","hebA",   "hbpA",   "hpbA",   "hepA",   "hpeA",   "he1e2A", "he2e1A")] <- 1
    allstates <- simout$allstates
    exp.index <- sum(1-allstates$s..)       # expected total # couples w/ 1 or more infection
    ## exp.index
    ## act.index <- sum(simdat$ser!=4)
    ## act.index
    exp.mb1. <- sum(allstates[,c("mb.", "hb1b2", "hbe", "hbp")])
    exp.me1. <- sum(allstates[,c("me.", "he1e2", "hep")])
    exp.f.b1 <- sum(allstates[,c("f.b", "hb2b1", "heb", "hpb")])
    exp.f.e1 <- sum(allstates[,c("f.e", "he2e1", "hpe")])
    ## sum(exp.mb1. + exp.me1. + exp.f.b1 + exp.f.e1) # total index cases
    exp.piGmb1. <-  sum(exp.mb1.) / exp.index
    exp.piGme1. <-  sum(exp.me1.) / exp.index
    exp.piGf.b1 <-  sum(exp.f.b1) / exp.index
    exp.piGf.e1 <-  sum(exp.f.e1) / exp.index
    out[1,,ii] <- c(exp.piGmb1.,exp.piGme1.,exp.piGf.b1,exp.piGf.e1)
    ## analyze live people
    adat <- simdat[grepl("A",simdat$cat.nm),]
    asA <- simout$allstates[grepl("A",simdat$cat.nm),]
    hh.log <- adat$ser==1
    mm.log <- adat$ser==2
    ff.log <- adat$ser==3
    ss.log <- adat$ser==4
    ## version 1
    mb1. <- rowSums(asA[mm.log|hh.log, c("mb.","hb1b2","hbe","hbp")])
    me1. <- rowSums(asA[mm.log|hh.log, c("me.","he1e2","hep")])
    m1.A <- rowSums(asA[mm.log|hh.log, c("mb.A","hb1b2A","hbeA","hbpA","me.A","he1e2A","hepA")])
    ## female
    f.b1 <- rowSums(asA[ff.log|hh.log, c("f.b","hb2b1","heb","hpb")])
    f.e1 <- rowSums(asA[ff.log|hh.log, c("f.e","he2e1","hpe")])
    f.1A <- rowSums(asA[ff.log|hh.log, c("f.bA","hb2b1A","hebA","hpbA","f.eA","he2e1A","hpeA")])
    ## number of inflated male before-couple index infections (in any couple)
    mb1.Infl <- sum( mb1. / m1.A)
    ## number of inflated male extra-couple index infections (in any couple)
    me1.Infl <- sum( me1. / m1.A)
    ## nufber of inflated male before-couple index infections (in any couple)
    f.b1Infl <- sum( f.b1 / f.1A)
    ## number of inflated male extra-couple index infections (in any couple)
    f.e1Infl <- sum( f.e1 / f.1A)
    IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
    piGb1. <- mb1.Infl / IndInfl
    piGe1. <- me1.Infl / IndInfl
    piG.b1 <- f.b1Infl / IndInfl
    piG.e1 <- f.e1Infl / IndInfl
    out[2,,ii] <- c(piGb1.,piGe1.,piG.b1,piG.e1)
    ## ## version 2
    mb1. <- rowSums(asA[mm.log|hh.log, c("mb.","hb1b2","hbe","hbp")])
    me1. <- rowSums(asA[mm.log|hh.log, c("me.","he1e2","hep")])
    m.A <- rowSums(asA[mm.log|hh.log, names(asA)[(grepl("m", names(asA)) | grepl("h", names(asA))) & grepl("A", names(asA))]])
    ## female
    f.b1 <- rowSums(asA[ff.log|hh.log, c("f.b","hb2b1","heb","hpb")])
    f.e1 <- rowSums(asA[ff.log|hh.log, c("f.e","he2e1","hpe")])
    f.A <- rowSums(asA[ff.log|hh.log, names(asA)[(grepl("f", names(asA)) | grepl("h", names(asA))) & grepl("A", names(asA))]])
    ## number of inflated male before-couple index infections (in any couple)
    mb1.Infl <- sum( mb1. / m.A)
    ## number of inflated male extra-couple index infections (in any couple)
    me1.Infl <- sum( me1. / m.A)
    ## nufber of inflated male before-couple index infections (in any couple)
    f.b1Infl <- sum( f.b1 / f.A)
    ## number of inflated male extra-couple index infections (in any couple)
    f.e1Infl <- sum( f.e1 / f.A)
    ## number of inflated index infections (1 in each couple with an infection)
    IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
######################################################################
    ## Proportion of index infections pooling across gender
    piGb1. <- mb1.Infl / IndInfl
    piGe1. <- me1.Infl / IndInfl
    piG.b1 <- f.b1Infl / IndInfl
    piG.e1 <- f.e1Infl / IndInfl
    out[3,,ii] <- c(piGb1.,piGe1.,piG.b1,piG.e1)
    ## ## version 3
    mb1. <- rowSums(asA[!ss.log, c("mb.","hb1b2","hbe","hbp")])
    me1. <- rowSums(asA[!ss.log, c("me.","he1e2","hep")])
    ## female
    f.b1 <- rowSums(asA[!ss.log, c("f.b","hb2b1","heb","hpb")])
    f.e1 <- rowSums(asA[!ss.log, c("f.e","he2e1","hpe")])
    all.infA <- rowSums(asA[!ss.log,names(asA)[(grepl("m", names(asA)) | grepl("f", names(asA)) | grepl("h", names(asA))) & grepl("A", names(asA))]])
    ## number of inflated male before-couple index infections (in any couple)
    mb1.Infl <- sum( mb1. / all.infA)
    ## number of inflated male extra-couple index infections (in any couple)
    me1.Infl <- sum( me1. / all.infA)
    ## nufber of inflated male before-couple index infections (in any couple)
    f.b1Infl <- sum( f.b1 / all.infA)
    ## number of inflated male extra-couple index infections (in any couple)
    f.e1Infl <- sum( f.e1 / all.infA)
    ## number of inflated index infections (1 in each couple with an infection)
    IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
######################################################################
    ## Proportion of index infections pooling across gender
    piGb1. <- mb1.Infl / IndInfl
    piGe1. <- me1.Infl / IndInfl
    piG.b1 <- f.b1Infl / IndInfl
    piG.e1 <- f.e1Infl / IndInfl
    out[4,,ii] <- c(piGb1.,piGe1.,piG.b1,piG.e1)

    ## ## version 4
    mb1. <- rowSums(asA[, c("mb.","hb1b2","hbe","hbp")])
    me1. <- rowSums(asA[, c("me.","he1e2","hep")])
    ## female
    f.b1 <- rowSums(asA[, c("f.b","hb2b1","heb","hpb")])
    f.e1 <- rowSums(asA[, c("f.e","he2e1","hpe")])
    all.infA <- rowSums(asA[,c("s..", names(asA)[(grepl("m", names(asA)) | grepl("f", names(asA)) | grepl("h", names(asA))) & grepl("A", names(asA))])])
    ## number of inflated male before-couple index infections (in any couple)
    mb1.Infl <- sum( mb1. / all.infA)
    ## number of inflated male extra-couple index infections (in any couple)
    me1.Infl <- sum( me1. / all.infA)
    ## nufber of inflated male before-couple index infections (in any couple)
    f.b1Infl <- sum( f.b1 / all.infA)
    ## number of inflated male extra-couple index infections (in any couple)
    f.e1Infl <- sum( f.e1 / all.infA)
    ## number of inflated index infections (1 in each couple with an infection)
    IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
######################################################################
    ## Proportion of index infections pooling across gender
    piGb1. <- mb1.Infl / IndInfl
    piGe1. <- me1.Infl / IndInfl
    piG.b1 <- f.b1Infl / IndInfl
    piG.e1 <- f.e1Infl / IndInfl
    out[5,,ii] <- c(piGb1.,piGe1.,piG.b1,piG.e1)
  }

bias <- apply(out[2:5,,] - out[rep(1,4),,], 1:2, mean)
variance <- apply(out[2:5,,] - out[rep(1,4),,], 1:2, var)

par(mfrow=c(1,4))
for(ii in 2:5)
  {
    hist(out[ii,1,] - out[1,1,1], col="black")
    abline(v=0, col ="red")
  }

tab <- signif(abind(bias, sqrt(variance), along = 3),3)
tab <- apply(tab, 1:2, function(x) paste(x[1], " (",x[2],")",sep=""))
write.csv(tab, file = "~/Documents/R files/discordant couples/reldurmod/pi estimator sim.csv")

write.csv(signif(bias,3), file = "~/Documents/R files/discordant couples/reldurmod/pi est bias.csv")
write.csv(signif(variance,3), file = "~/Documents/R files/discordant couples/reldurmod/pi est variance.csv")


