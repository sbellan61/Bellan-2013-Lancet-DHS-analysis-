rm(list=ls())
library(coda)
library(abind)
library(multicore)


## setwd("~/Downloads/cloudjun/WA-Deaths-20120609-05:18/files/")
## load("workspace.Rdata")
## mcmc.out1 <- mcmc.out
## setwd("~/Downloads/cloudjun/WA-Deaths-20120609-05:22/files/")
## load("workspace.Rdata")
## mcmc.out2 <- mcmc.out
## setwd("~/Downloads/cloudjun/WA-Deaths-20120609-05:25/files/")
## load("workspace.Rdata")
## mcmc.out3 <- mcmc.out
## mcmc.out <- as.mcmc.list(c(mcmc.out1,mcmc.out2, mcmc.out3))


setwd("~/Downloads/cloudjun/Ethiopia-Deaths-20120608-14:36/files/")
load("workspace.Rdata")
mcmc.out1 <- mcmc.out
setwd("~/Downloads/cloudjun/Ethiopia-Deaths-20120608-23:39/files/")
load("workspace.Rdata")
mcmc.out2 <- mcmc.out
mcmc.out <- as.mcmc.list(c(mcmc.out1,mcmc.out2))

sum.wrap <- function(col.ind, xx)
  {
    summary(xx[,col.ind], quantiles = c(.025, .5, .975))$quant
  }

pars <- abind(mclapply(1:ncol(mcmc.out[[1]]), sum.wrap,  xx = mcmc.out, mc.cores=1), along = 0, new.names=colnames(mcmc.out[[1]]))

## never save over old files
file.name <- paste("pars", format(Sys.time(), "%Y%m%d"), sep = "-")

pars.name <- paste(file.name, ".Rdata",sep="")
chain.name <- paste(file.name, "chains.Rdata",sep="")

out.csv <- signif(pars,3)
out.csv <- data.frame(median = out.csv[,2], CI95 = paste("(",out.csv[,1],", ",out.csv[,3],")", sep=""))
write.csv(out.csv, file="out.csv")

save(pars, file = pars.name)
save(mcmc.out, file = chain.name)
