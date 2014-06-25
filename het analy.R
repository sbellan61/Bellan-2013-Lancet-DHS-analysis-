rm(list=ls())

simstr <- "runs"
dirnm <- paste("~/Documents/R files/discordant couples/sims/", simstr, sep="")

flnms <- c("Zambia Hom","Zambia Het","Kenya Hom","Kenya Het")
setwd(dirnm)
dirs <- list.files()
to.do <- data.frame(!grepl("het",dirs) & grepl("Zambia",dirs)& grepl("sim",dirs),
           grepl("het",dirs) & grepl("Zambia",dirs) & grepl("sim",dirs),
           !grepl("het",dirs) & grepl("Kenya",dirs) & grepl("sim",dirs),
           grepl("het",dirs) & grepl("Kenya",dirs) & grepl("sim",dirs))


                     
for(ff in 1:4)
  {
    tdirs <- dirs[to.do[,ff]]

    nnn <- length(tdirs)
    parsout <- array(NA, c(12,4,nnn))

    for(iii in 1:nnn)
      {
        to.load <- paste(tdirs[iii],"/files/workspace.Rdata", sep = "")
        load(to.load)
        ## estimated breakdown in simulated data (alive)
        show <- c(rownames(pars)[c(1:5,7)], "pibUA","pieUA","pipUA", "piUbA","piUeA","piUpA")
        pars[show,]
        parsout[,1:3,iii] <- pars[show,]
        if(iii==1)
          {
            colnames(parsout) <- c(colnames(pars), "trueval")
            rownames(parsout) <- show
          }
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
        parsout[,4,iii] <- c(simpars[1:5],exp(simpars[6])*simpars[5], pibUA.r, pieUA.r, pipUA.r, piUbA.r, piUeA.r, piUpA.r)
      }

    if(heter)
      {
        scalars <- c(rep(sqrt(exp(1)),4), rep(1,8))
      }else{
        scalars <- rep(1,12)
      }
    
    biass <- rowMeans((parsout[,2,] - scalars*parsout[,4,]) / parsout[,2,])

    coverage <- rowMeans(parsout[,1,] < scalars*parsout[,4,] & parsout[,3,] > scalars*parsout[,4,])

    if(ff==1)
      {
        out.tab <- c(paste(signif(biass,2), " (",signif(coverage,2)*100,"%)", sep = ""))
      }else{
        out.tab <- rbind(out.tab, c(paste(signif(biass,2), " (",signif(coverage,2)*100,"%)", sep = "")))
      }
    ## Create a data frame that compares true routes of transmission
    ## for observed individuals to estimates
    rownames(parsout)[7:12] <- paste("prop", rep(c("male","female"), each = 3), rep(c("b","e","p"), 2))

    pdf(paste("~/Documents/R files/discordant couples/sims/",flnms[ff], ".pdf", sep = ""), w = 8, h = 10)
    layout(matrix(c(1,2,7,10,3,4,8,11,5,6,9,12),4,3))
    par(mar = c(3,.5,2,.5))
    for(jjj in 1:nrow(parsout))
      {
        if(jjj < 7)
          {
            xlim <- range(parsout[jjj,,])
          }else{
            xlim <- c(0,1)
          }
        plot(0,0, type = "n", xlim = xlim, ylim = c(.5, nnn + .5), bty = "n", ylab = "", yaxt = "n",
             main = paste(rownames(parsout)[jjj], ": bias =", signif(biass[jjj],2)))
        arrows(parsout[jjj,1,], 1:nnn, parsout[jjj,3,], 1:nnn, code = 3, angle = 90, length = .03)
        points(parsout[jjj,2,], 1:nnn, pch = 19)
        if(jjj %in% 1:4 & heter)
          {
            scalar <- sqrt(exp(1))
          }else{
            scalar <- 1
          }
        points(scalar*parsout[jjj,4,], 1:nnn, pch = 19, col = "red")
      }
    dev.off()
  }

rownames(out.tab) <- flnms
colnames(out.tab) <- rownames(parsout)
write.csv(out.tab, file="~/Documents/R files/discordant couples/sims/bias cov.csv")
