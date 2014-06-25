## Simulate DHS like data to be fit to make sure everything is working
## fine.
###################################################################### 
## Steve Bellan, February 2012
## sbellan@berkeley.edu / steve.bellan@gmail.com
###################################################################### 
## Set directory where output is to be written
if(getwd()=="/home/ubuntu/files")       # if on cloud
  {
    setwd("~/files/")
  }else{                                # if on sb mac
    setwd("~/Dropbox/disc mod backup/Couples Model Revision 120811/R scripts & Data Files/")
  }
## load cleaned DHS data sets for all countries
load("alldhs.Rdata")                    
## load cleaned epidemic curves (w/ art-normalized prevalence for
## 15-49 by M & F
load("allepicm.Rdata")
load("allepicf.Rdata")
## load survival curve cs, monthly survival
load("csurv.Rdata")
######################################################################
## non western-africa indices are
######################################################################
##      [,1]  [,2]       [,3]    [,4]      [,5]     [,6]     [,7]        [,8]
## [1,] "DRC" "Ethiopia" "Kenya" "Lesotho" "Malawi" "Rwanda" "Swaziland" "WA"
##      [,9]     [,10]     
## [1,] "Zambia" "Zimbabwe"
######################################################################
## sample sizes
##       DRC  Ethiopia     Kenya   Lesotho    Malawi    Rwanda Swaziland        WA 
##      1199      1650      1620      1099      3045      1749       431      7915 
##    Zambia  Zimbabwe 
##      1598      1271 

source("pcalc5.R")
if(!"group.ind" %in% ls())        group.ind <- 9
if(!"survive" %in% ls())        survive <- T
if(!"K.sim" %in% ls())  K.sim <- 1000
## Use this ds' couple & sex duration distributions in simulation.
group <- levels(dat$group)[group.ind]
ds <- levels(dat$ds)[grepl(group, levels(dat$ds))]   # data set we're working on
dat <- dat[dat$ds %in% ds,]           # only looking at that data set
print(paste("real",group,"serostatus"))
print(table(dat$ser))


odat <- dat[sample(1:nrow(dat), K.sim, replace = T),]
odat$ds <- paste(group, "simulated")
odat$group <- odat$ds
odat$ser <- NA
## Get before couple duration bd, where bd = max(mbd,fbd)+1
odat$bd <- apply(cbind(odat$tmar-odat$tms,odat$tmar-odat$tfs), 1, max) + 1
## Get couple duration cd
odat$cd <- odat$tint - odat$tmar
## Sexual partners per year of activity (s p acquisition rate: mspar, fspar)
odat$mspar <- odat$mlsp / (odat$tint - odat$tms + 1) * 12
odat$fspar <- odat$flsp / (odat$tint - odat$tfs + 1) * 12


if(!"survive" %in% ls())  survive <- T
## if(!"simpars" %in% ls()) # if sim pars is not already loaded (because normally sourcing from areldisc
##   {
        simpars <- c(bmb = .01, bfb = .02,
                 bme = .005, bfe=.005,
                 bmp=.015, lrho=0)
##   }
load("pars.arr.Rdata")                  # load simpars from last fitted analysis for that country

if(!"heter" %in% ls())  heter <- F      # default for individual heterogeneity in simulation
if(heter)
  {
    simpars[1:6] <- pars.arr[1:6,2, group.ind]
    simpars[5] <- .9* simpars[5]
    simpars[1:2] <- .5*simpars[1:2]
    simpars[3:4] <- .7*simpars[3:4]
    simpars[6] <- 0
  }

K <- nrow(odat)
if(survive)
  {
    simout <- pcalc(simpars, dat = odat, trace = F, give.pis=T, heter = heter, lrho.sd = lrho.sd,
                    survive = survive, sim = T, browse=F)
  }else{
    simout <- pcalc(simpars, dat = odat, trace = F, give.pis=F, heter = heter, lrho.sd = lrho.sd,
                    survive = survive, sim = T, browse=F)
  }
simdat <- simout$dat
simdat$ser <- NA
simdat$ser[simdat$cat.nm == "s..A"] <- 4
simdat$ser[simdat$cat.nm %in% c("mb.D", "me.D", "mb.A", "me.A")] <- 2
simdat$ser[simdat$cat.nm %in% c("f.bD", "f.eD", "f.bA", "f.eA")] <- 3
simdat$ser[simdat$cat.nm %in% c("hb1b2D", "hb2b1D", "hbeD","hebD",   "hbpD",   "hpbD",   "hepD",   "hpeD",   "he1e2D", "he2e1D",
                                "hb1b2A", "hb2b1A", "hbeA","hebA",   "hbpA",   "hpbA",   "hepA",   "hpeA",   "he1e2A", "he2e1A")] <- 1
simdat$alive <- T
simdat$alive[simdat$cat.nm %in% levels(simdat$cat.nm)[grepl("D",levels(simdat$cat.nm))]] <- F
allstates <- simout$allstates
## print(paste("simulated",group,"serostatus"))
## print(table(simdat$ser))
## print("dead couples")
## print(table(simdat$ser[!simdat$alive]))
## print("live couples")
## print(table(simdat$ser[simdat$alive]))
K <- nrow(simdat)
print("real dat serostatus breakdown")
print(round(xtabs(~ser,dat) / nrow(dat), 3))
print("simulated dat serostatus breakdown")
print(round(xtabs(~ser, simdat, subset = alive) / sum(simdat$alive), 3))
simdat <- data.frame(simdat, m.het = simout$m.het, f.het = simout$f.het)
save(simdat, file = "simulated dat.Rata")
save.image("workspace.Rdata")
simdat <- simdat[simdat$alive,]

cols <- c("yellow","green","purple","dark gray")
par(mar= rep(5,4))
plot(log(simdat$m.het), log(simdat$f.het), pch = 21, col =  cols[simdat$ser], cex = 1, bty = "n",
     xlab = "male log RF", ylab = "female log RF", lwd = 1.5)

## reload dhs data and add the columns necessry to bind it with the simulation
load("alldhs.Rdata")
dat$bd <- apply(cbind(dat$tmar-dat$tms,dat$tmar-dat$tfs), 1, max) + 1
dat$cd <- dat$tint - dat$tmar
dat$cat <- NA
dat$cat.nm <- NA
dat$alive <- NA
dat$mspar <- NA
dat$fspar <- NA
dat$m.het <- NA
dat$f.het <- NA
dat <- rbind(dat, simdat)
save(dat, file="alldhs plus sim.Rdata")





