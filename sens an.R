## Conduct a sensitivity analysis of results to couples that were
## missing serostatus. Predict their serostatus using the fitted
## parameters and add them to the results.
rm(list=ls())


## load raw data
load("~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/R scripts & Data Files/alldhs raw.Rdata")
show <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","mage","fage",
          "epic.ind","epic.nm","group", "remser", "noser","rem")
## Need to have ysa be >= mardur to avoid a nonnegative time b/w start
## of sex & start of marriage create serostatus of couple so for
## individuals where time of first sex was after marriage, set it to
## time of marriage.
nna <- !is.na(allraw$tms) & !is.na(allraw$tfs) & !is.na(allraw$tmar)
allraw$tms[nna & allraw$tms>allraw$tmar] <- allraw$tmar[nna & allraw$tms>allraw$tmar]
allraw$tfs[nna & allraw$tfs>allraw$tmar] <- allraw$tmar[nna & allraw$tfs>allraw$tmar]
allraw$fysa[nna & !allraw$fysa>=allraw$mardur] <- allraw$mardur[nna & !allraw$fysa>=allraw$mardur]
allraw$mysa[nna & !allraw$mysa>=allraw$mardur] <- allraw$mardur[nna & !allraw$mysa>=allraw$mardur]
araw <- allraw[,show]

groups <- unique(allraw$group)

result.dir <- "~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/sims/runs/"
result.out.dir <- "~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/"
setwd(result.dir)
dirnms <- list.files()
dirnms <- dirnms[!grepl("sim",dirnms)]
ndir <- length(dirnms)

stab.arr <- array(NA, c(16, 8, 10))
for(c.ind in 1:ndir)
  {
    ## load data set
    tempdir <- paste(result.dir, dirnms[c.ind], "/files/", sep ="")
    setwd(tempdir)
    load("workspace.Rdata")
    source("~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 120811/R scripts & Data Files/pcalc5.R")    
    tgr <- dat$group[1]                 # which group
    if(c.ind==1)
      {
        groups <- tgr
      }else{
        groups <- c(groups,tgr)
      }
    tdat <- araw[araw$group==tgr & (araw$remser | !araw$rem),] # select augmented data
    tdat$bd <- apply(cbind(tdat$tmar-tdat$tms,tdat$tmar-tdat$tfs), 1, max) + 1
    ## Get couple duration cd
    tdat$cd <- tdat$tint - tdat$tmar
    tsim <- pcalc(pars[1:6,2], tdat, sim = T, survive = T, browse = F, # augmented
                  give.pis=T, trace=T)
    dsim <- pcalc(pars[1:6,2], dat, sim = T, survive = T, browse = F, # normal
                  give.pis=T, trace=T)
    pis <- names(tsim$pop.avs)
    rownames(stab.arr) <- pis
    stab.arr[,,c.ind] <- cbind(pars[pis,], andat = t(dsim$pop.avs), augdat = t(tsim$pop.avs),
                                    absbias = (t(tsim$pop.avs) - t(dsim$pop.avs)),
                                    bias = (t(tsim$pop.avs) - t(dsim$pop.avs)) / t(dsim$pop.avs),
                                    in.ci = t(tsim$pop.avs) < pars[pis,3] & t(tsim$pop.avs) > pars[pis,1])
    pdf(paste("~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/sens",tgr,".pdf"))
    par(mfrow=c(3,2))
    breaks <- seq(-36,800, by = 12*5)
    plot(density(tdat$mage, bw = 24), col = "red")
    lines(density(dat$mage, bw = 24), col = "black") 
    plot(density(tdat$fage, bw = 24), col = "red")
    lines(density(dat$fage, bw = 24), col = "black")
    plot(density(tdat$tmar - tdat$tms, bw = 24), col = "red")
    lines(density(dat$tmar - dat$tms, bw = 24), col = "black")
    plot(density(tdat$tmar - tdat$tfs, bw = 24), col = "red")
    lines(density(dat$tmar - dat$tfs, bw = 24), col = "black")
    plot(density(tdat$mardur.mon, bw = 24), col = "red")
    lines(density(dat$mardur.mon, bw = 24), col = "black")
    dev.off()
##         mean(tdat$mage)
##         mean(dat$mage)
##         mean(tdat$mardur.mon)
##         mean(dat$mardur.mon)
##         mean(tdat$tmar - tdat$tms)
##         mean(dat$tmar - dat$tms)        
  }
colnames(stab.arr) <- c(colnames(pars), "andat", "augdat", "absbias", "bias", "in.ci")

stab.arr <- signif(stab.arr,3)

sbias <- stab.arr[,"bias",]
colnames(sbias) <- levels(dat$group)
sbias


save(stab.arr, file="~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/stab.arr.Rdata")
write.csv(t(sbias), file="~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/sens bias remser.csv")


######################################################################
## Sensitivity analyses to various things
rm(list=ls())
dirnms <- c("main","parv","lcarv","parv-lcarv","fsd")
flds <- paste("~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/",dirnms,"/",sep="")
nms <- dirnms
nms[4] <- "parv.lcarv"

ah.m <- data.frame(matrix(NA,30,6))
ah.f <- data.frame(matrix(NA,30,6))
fb.m <- data.frame(matrix(NA,10,6))
fb.f <- data.frame(matrix(NA,10,6))
for(ii in 1:length(flds))
  {
    temp <- read.csv(paste(flds[ii],"all historical breakdown.csv", sep = ""))
    #assign(paste("ah.",nms[ii], sep=""), temp)
    if(ii==1)
      {
        ah.m[,c(1,2)] <- temp[,c(1,3)]
        ah.f[,c(1,2)] <- temp[,c(1,6)]        
      }else{
        whicher <- ah.m[,1] %in% temp[,1]
        ah.m[whicher,ii+1] <- as.character(temp[,3])
        ah.f[whicher,ii+1] <- as.character(temp[,6])
      }
    temp <- read.csv(paste(flds[ii],"proportion incidence future breakdown.csv", sep = ""))
    if(ii==1)
      {
        fb.m[,c(1,2)] <- temp[,c(1,2)]
        fb.f[,c(1,2)] <- temp[,c(1,3)]        
      }else{
        whicher <- fb.m[,1] %in% temp[,1]
        fb.m[whicher,ii+1] <- as.character(temp[,2])
        fb.f[whicher,ii+1] <- as.character(temp[,3])        
      }
    #assign(paste("fb.",nms[ii], sep=""), temp)
  }
colnames(ah.m) <- c("country",nms)
colnames(ah.f) <- c("country",nms)
ah <- cbind(ah.m,ah.f[,-1])
colnames(fb.m) <- c("country",nms)
colnames(fb.f) <- c("country",nms)
fb <- cbind(fb.m,fb.f[,-1])

write.csv(ah, "~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/all historical sens an.csv")
write.csv(fb, "~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/Output/future breakdown sens an.csv")
