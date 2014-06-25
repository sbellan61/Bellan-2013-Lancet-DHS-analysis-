load("~/Documents/R files/discordant couples/cloud final/deaths/Zambia-Deaths-20120608-12:05/files/workspace.Rdata")

###################################################################### 
## Figure 2 - prob an infection is from extracouple sex by serostatus
## by ysa and reldur
###################################################################### 

medpars <- pars[parnames,2]
pis <- pcalc(medpars, dat = dat, trace = T, give.pis=T, survive=survive)$pis
######################################################################

######################################################################
cex <- .8
cex.ax <- 1.2
rmp <- colorRamp(c("yellow","red"))     #create color ramp
pal <- colorRampPalette(c("yellow","red"))     #create color palette (for legend)
## pdf("Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.pdf",
##     width = 6.5, height = 4
pdf("~/Documents/R files/discordant couples/cloud final/zambia 2b.pdf", width = 6, h = 5)
par(mar = c(4,4,4,4), oma = rep(0,4))
                                        # do it for each data set since the interview times were different for
                                        # same countries and so the plot shold show how long the couples were
                                        # together too
for(dd in unique(dat$epic.nm)) 
  {
    cc <- dat$epic.ind[dat$epic.nm==dd][1]   #find epidemic curve for that data set
    epic.col <- "blue"
    ylim <- c(0, 30) #max(dat$m.bef.pm, dat$f.bef.pm))/12
    xlim <- c(1975,2012)
    ylab <-  "YSA before couple formation"
    xlab <- "date of couple formation"
    ii <- 2
    cols.show <- rgb(rmp(pis$piCe.A[dat$ser==ii & dat$epic.ind==cc]), max = 255)

    plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii & dat$epic.ind==cc],
         1/12*(dat$tmar-dat$tms)[dat$ser==ii & dat$epic.ind==cc],
         col = cols.show, ylab = "male YSA before couple formation",
         xlab = xlab,
         pch = 19, cex = cex, axes = F,
         xlim = xlim, ylim = ylim, bty = "n", 
         main = "Zambia: male HIV+ serodiscordant couples")
    xs <- epicm[,1]/12 + 1900
    axis(2, at = seq(0, 30, by = 10), las = 2)
    axis(1, at = seq(1980, 2010, by = 10), las = 2)
    lines(xs[xs>1975], epicf[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
    axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
    mtext("Female infectious HIV pop. prevalence (%)", side = 4, outer = F , line = 2.5, adj = .5, cex = 1)
} # end loop through pooled data sets
dev.off()
######################################################################
