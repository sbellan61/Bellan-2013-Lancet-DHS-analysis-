rm(list=ls())
setwd("~/Dropbox/disc mod backup/Couples Model Revision 120811/R scripts & Data Files/")
load("alldhs.Rdata")
load("allepicm.Rdata")
load("allepicf.Rdata")
load("pars.arr.Rdata")
load("csurv.Rdata")
load("ds.name.Rdata")
######################################################################
## data grouping indices are
######################################################################
##      [,1]  [,2]       [,3]    [,4]      [,5]     [,6]     [,7]        [,8]
## [1,] "DRC" "Ethiopia" "Kenya" "Lesotho" "Malawi" "Rwanda" "Swaziland" "WA"
##      [,9]     [,10]      [,11]             
## [1,] "Zambia" "Zimbabwe" "simulated"
## Get before couple duration bd, where bd = max(mbd,fbd)
dat$bd <- apply(cbind(dat$tmar-dat$tms,dat$tmar-dat$tfs), 1, max)
## Get couple duration cd
dat$cd <- dat$tint - dat$tmar

country <- 9
ds.nm[country]

hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
pars <- pars.arr[hazs,2,country]
pars


## load("~/Documents/R files/discordant couples/cloud final/deaths/Malawi-Deaths-20120516-15:31/files/pars-20120516.Rdata")
## load("~/Documents/R files/discordant couples/cloud final/deaths/Lesotho-Deaths-20120516-03:25/files/pars-20120516.Rdata")
## pars[hazs,2]

dat <- dat[dat$group==ds.nm[country],]
dat$mds <- dat$tmar - dat$tms
dat$fds <- dat$tmar - dat$tfs
head(dat)


## Age-at-seroconversion dependent Weibull survival times.
ageweib <- function(age)
  {
    shp <- 2.3
    scl <- 2000/shp/(age/12)^.53
    round(rweibull(length(age), shape = shp, scale = scl),0)
  }

## Simulation function
sim.fn <- function(pars, dat, browse = F,
                   vfreq = 50)      # how often to show progress
  {
    K <- nrow(dat)
    bmb <- as.numeric(pars["bmb"])
    bfb <- as.numeric(pars["bfb"])
    bme <- as.numeric(pars["bme"])
    bfe <- as.numeric(pars["bfe"])
    bmp <- as.numeric(pars["bmp"])
    bfp <- as.numeric(pars["bfp"])
    ## serostatus at end
    mser <- rep(0, nrow(dat))
    fser <- rep(0, nrow(dat))
    ## date of death (CMC)
    mdod <- rep(NA, nrow(dat))
    fdod <- rep(NA, nrow(dat))
    ## date of infection (CMC)
    mdoi <- rep(NA, nrow(dat))
    fdoi <- rep(NA, nrow(dat))
    ## route of infection
    mcoi <- factor(rep(NA, nrow(dat)), levels = c("b","e","p","-"))
    fcoi <- factor(rep(NA, nrow(dat)), levels = c("b","e","p","-"))
    if(browse) browser()
    for(ii in 1:nrow(dat))              # loop over couples
      {
        if(ii %% vfreq == 0) print(paste("On couple", ii, "of", nrow(dat)))
        epic.ind.temp <- dat$epic.ind[ii]
        ## keep simulating until date of death > date of couple
        ## formation, i.e. condition on couple forming.
        finished <- FALSE
        ## add IF STATEMENT FOR IF HE HAD SEX BEFORE MARRIAGE
        while(!finished) 
          {
            ## male before marriage: get bernoulli probability of
            ## infection in each month of sexual activity before
            ## marriage
            m.inf.bef <- rbinom(dat$tmar[ii] - dat$tms[ii],1,
                                prob = 1 - exp(-bmb*epicf[dat$tms[ii]:(dat$tmar[ii]-1), epic.ind.temp]))
            ## if he gets infected in 1 or more months
            if(sum(m.inf.bef)>0)
              {
                ## change serostatus to HIV+ & use earliest
                ## infection as date of infection. Then figure out
                ## date of death from age-at-seroconversion
                ## dependent Weibull survival times.
                mser[ii] <- 1
                mdoi[ii] <- dat$tms[ii] + min(which(m.inf.bef==1)) - 1
                mdod[ii] <- mdoi[ii] + ageweib(dat$mage[ii] - (dat$tint[ii] - mdoi[ii]))
                mcoi[ii] <- "b"
                if(mdod[ii] > dat$tmar[ii]) finished <- TRUE # if death is after couple formation, finished
              }else{                # if no infections, end while loop
                finished <- TRUE
              }
          }
        ## keep simulting until date of death > date of couple
        ## formation, i.e. condition on couple forming.
        finished <- FALSE
        while(!finished) 
          {
            ## female before marriage: get bernoulli probability of
            ## infection in each month of sexual activity before
            ## marriage
            f.inf.bef <- rbinom(dat$tmar[ii] - dat$tfs[ii],1,
                                prob = 1 - exp(-bfb*epicm[dat$tfs[ii]:(dat$tmar[ii]-1), epic.ind.temp]))
            ## if she gets infected in 1 or more months
            if(sum(f.inf.bef)>0)
              {
                ## change serostatus to HIV+ & use earliest
                ## infection as date of infection. Then figure out
                ## date of death from age-at-seroconversion
                ## dependent Weibull survival times.
                fser[ii] <- 1
                fdoi[ii] <- dat$tfs[ii] + min(which(f.inf.bef==1)) - 1
                fdod[ii] <- fdoi[ii] + ageweib(dat$fage[ii] - (dat$tint[ii] - fdoi[ii]))
                fcoi[ii] <- "b"
                if(fdod[ii] > dat$tmar[ii]) finished <- TRUE # if death is after couple formation, finished
              }else{                # if no infections, end while loop
                finished <- TRUE
              }
          }
        if(mser[ii] + fser[ii] < 2)     # if both aren't infected, simulate marriage stuff
          {
            end <- min(mdod[ii], fdod[ii], dat$tint[ii], na.rm = TRUE)
            tt <- dat$tmar[ii]
            ## While both are alive and both aren't infected, run infection loop.
            while((tt <= end) & (mser[ii] + fser[ii] < 2))
              {
                ## infections from partners, length 2 vector,
                ## bernoulli, c(male inf by part, female inf by part);
                from.part <- rbinom(2, 1, prob = 1 - exp(-c(bmp, bfp)))
                from.part <- from.part * c(fser[ii], mser[ii]) # condition on partner being infected
                ## if new m inf from partner
                if(mser[ii]==0 & from.part[1] == 1)
                  {
                    mser[ii] <- 1
                    mdoi[ii] <- tt
                    mdod[ii] <- tt + ageweib(dat$mage[ii] - (dat$tint[ii] - tt))
                    mcoi[ii] <- "p"
                  }
                ## if new f inf from partner
                if(fser[ii]==0 & from.part[2] == 1)
                  {
                    fser[ii] <- 1
                    fdoi[ii] <- tt
                    fdod[ii] <- tt + ageweib(dat$fage[ii] - (dat$tint[ii] - tt))
                    fcoi[ii] <- "p"
                  }
                ## extracouple infections c(m,f)
                exc <- rbinom(2, 1, prob = c(1 - exp(-bme*epicf[tt, epic.ind.temp]),
                                             1 - exp(-bfe*epicm[tt, epic.ind.temp])))
                ## if new m inf from extracouple
                if(mser[ii]==0 & exc[1] == 1)
                  {
                    mser[ii] <- 1
                    mdoi[ii] <- tt
                    mdod[ii] <- tt + ageweib(dat$mage[ii] - (dat$tint[ii] - tt))
                    mcoi[ii] <- "e"
                  }
                ## if new f inf from extracouple
                if(fser[ii]==0 & exc[2] == 1)
                  {
                    fser[ii] <- 1
                    fdoi[ii] <- tt
                    fdod[ii] <- tt + ageweib(dat$fage[ii] - (dat$tint[ii] - tt))
                    fcoi[ii] <- "e"
                  }
                end <- min(mdod[ii], fdod[ii], dat$tint[ii], na.rm = TRUE)
                tt <- tt+1
              }
          }                             # end marriage if statement
      }                                 # end couple loop
    dat <- data.frame(dat, mser, fser, mdoi, fdoi, mdod, fdod, mcoi, fcoi)
    dat$malive <- rep(T, nrow(dat))
    dat$malive[dat$mdod <= dat$tint] <- F
    dat$falive <- rep(T, nrow(dat))
    dat$falive[dat$fdod <= dat$tint] <- F
    dat$alive <- dat$malive & dat$falive
    dat$mcoi[is.na(dat$mcoi)] <- "-"
    dat$fcoi[is.na(dat$fcoi)] <- "-"
    dat$ser[dat$mser == 1 & dat$fser == 1] <- 1
    dat$ser[dat$mser == 1 & dat$fser == 0] <- 2
    dat$ser[dat$mser == 0 & dat$fser == 1] <- 3
    dat$ser[dat$mser == 0 & dat$fser == 0] <- 4    
    return(dat)
  }


## Create a function that takes the output from sim.fn and makes it
## look like the sample proportions matching the output of
## pcalc$allstates
transf <- function(stuff)
  {
    with(stuff,
         {
           K <- nrow(stuff)
           s.. <- sum(mser+fser==0) / K
           mb. <- sum(mcoi == "b" & fser==0) / K
           me. <- sum(mcoi == "e" & fser==0) / K
           f.b <- sum(fcoi == "b" & mser==0) / K
           f.e <- sum(fcoi == "e" & mser==0) / K
           hb1b2 <- sum(mcoi == "b" & fcoi == "b" & (mdoi < fdoi)) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi)) / K
           hb2b1 <- sum(mcoi == "b" & fcoi == "b" & (mdoi > fdoi)) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi)) / K
           hbe <- sum(mcoi == "b" & fcoi == "e") / K
           heb <- sum(mcoi == "e" & fcoi == "b") / K
           hbp <- sum(mcoi == "b" & fcoi == "p") / K
           hpb <- sum(mcoi == "p" & fcoi == "b") / K
           hep <- sum(mcoi == "e" & fcoi == "p") / K
           hpe <- sum(mcoi == "p" & fcoi == "e") / K
           he1e2 <- sum(mcoi == "e" & fcoi == "e" & (mdoi < fdoi)) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi)) / K
           he2e1 <- sum(mcoi == "e" & fcoi == "e" & (mdoi > fdoi)) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi)) / K
           ## And alive
           mb.A <- sum(mcoi == "b" & fser==0 & alive) / K
           me.A <- sum(mcoi == "e" & fser==0 & alive) / K
           f.bA <- sum(fcoi == "b" & mser==0 & alive) / K
           f.eA <- sum(fcoi == "e" & mser==0 & alive) / K
           hb1b2A <- sum(mcoi == "b" & fcoi == "b" & (mdoi < fdoi) & alive) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi) & alive) / K
           hb2b1A <- sum(mcoi == "b" & fcoi == "b" & (mdoi > fdoi) & alive) / K + .5 * sum(mcoi == "b" & fcoi == "b" & (mdoi == fdoi) & alive) / K
           hbeA <- sum(mcoi == "b" & fcoi == "e" & alive) / K
           hebA <- sum(mcoi == "e" & fcoi == "b" & alive) / K
           hbpA <- sum(mcoi == "b" & fcoi == "p" & alive) / K
           hpbA <- sum(mcoi == "p" & fcoi == "b" & alive) / K
           hepA <- sum(mcoi == "e" & fcoi == "p" & alive) / K
           hpeA <- sum(mcoi == "p" & fcoi == "e" & alive) / K
           he1e2A <- sum(mcoi == "e" & fcoi == "e" & (mdoi < fdoi) & alive) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi) & alive) / K
           he2e1A <- sum(mcoi == "e" & fcoi == "e" & (mdoi > fdoi) & alive) / K + .5 * sum(mcoi == "e" & fcoi == "e" & (mdoi == fdoi) & alive) / K
           out <- c(s.., mb., me., f.b, f.e, hb1b2, hb2b1, hbe, heb, hbp, hpb, hep, hpe, he1e2, he2e1, mb.A,
                    me.A,  f.bA,  f.eA,  hb1b2A, hb2b1A, hbeA,  hebA,  hbpA, hpbA,  hepA,  hpeA,  he1e2A, he2e1A)
           return(out)
         })
  }


## State probability formulation
pcalc <- function(pars, dat, browse = F, compars = NULL, # compare lprob for true pars with other pars (fitted)
                  give.pis = F,              # return individual pi values for couples (for outside mcmc)
                  sim = F,                   # just use this to simulate data given parameters? (outputs pser.a)
                  survive = T,               # account for survival in analysis
                  cond.sim = F,              # only simulate individuals that will live
                  trace = T) # only do certain calculations when tracing parameters (i.e. for non-thinned versions)
  {
    if(browse) browser()
    K <- nrow(dat)
    if(!sim)
      {
        hh.log <- dat$ser==1
        mm.log <- dat$ser==2
        ff.log <- dat$ser==3
        ss.log <- dat$ser==4
      }
    if(sum(pars[1:5]<0)>0)                   #if any parameters are <0 then the model must be rejected so we return logprob =-Inf
      {
        probs <- NA
        lprob <- -Inf
        pop.avs <- NA
        proj12 <- NA
        pser.a <- NA
        pser <- NA
        rrs <- NA
      }else{
        bmb <- pars[["bmb"]]
        bfb <- pars[["bfb"]]
        bme <- pars[["bme"]]
        bfe <- pars[["bfe"]]
        bmp <- pars[["bmp"]]
        if("lrho" %in% names(pars))
          {
            rho <- exp(pars[["lrho"]]) # feeding in log(rho)
            bfp <- bmp * rho
          }else{
            bfp <- pars[["bfp"]]
          }
        # L stands for *L*ast iteration
        s..L <- rep(1,K)                # s: concordant negative
        mb.L <- rep(0, K)               # m: male positive
        me.L <- rep(0, K)        
        f.bL <- rep(0, K)               # f: female positive
        f.eL <- rep(0, K)        
        hb1b2L <- rep(0, K)               # h: concordant positive
        hb2b1L <- rep(0, K)               # b1b2 is both inf before, but female 1st, & vice versa
        hbeL <- rep(0, K)
        hebL <- rep(0, K)               # 2nd & 3rd character give route of transmission for M & F, respectively
        hbpL <- rep(0, K) 
        hpbL <- rep(0, K)
        hepL <- rep(0, K)
        hpeL <- rep(0, K)
        he2e1L <- rep(0, K)
        he1e2L <- rep(0, K)
        ## initiate vectors to update based on *L*ast state
        s.. <- rep(1,K)                 
        mb. <- rep(0, K)               
        me. <- rep(0, K)        
        f.b <- rep(0, K)               
        f.e <- rep(0, K)        
        hb1b2 <- rep(0, K)
        hb2b1 <- rep(0, K)                       
        hbe <- rep(0, K)
        heb <- rep(0, K)               
        hbp <- rep(0, K) 
        hpb <- rep(0, K)
        hep <- rep(0, K)
        hpe <- rep(0, K)
        he1e2 <- rep(0, K)
        he2e1 <- rep(0, K)        
        # i.e., hbp is a ++ couple in which the male was inf *b*efore
        # couple formation & the female by her *p*artner
        if(survive)
          {
            # A stands for *A*live, i.e. joint probability of serostatus
            # and both partners being alive at sampling
            mb.AL <- rep(0, K)
            me.AL <- rep(0, K)        
            f.bAL <- rep(0, K)
            f.eAL <- rep(0, K)        
            hb1b2AL <- rep(0, K)
            hb2b1AL <- rep(0, K)
            hbeAL <- rep(0, K)
            hebAL <- rep(0, K)
            hbpAL <- rep(0, K)
            hpbAL <- rep(0, K)
            hepAL <- rep(0, K)
            hpeAL <- rep(0, K)
            he1e2AL <- rep(0, K)
            he2e1AL <- rep(0, K)
            ## initiate vectors to update based on *L*ast state
            mb.A <- rep(0, K)
            me.A <- rep(0, K)        
            f.bA <- rep(0, K)
            f.eA <- rep(0, K)        
            hb1b2A <- rep(0, K)
            hb2b1A <- rep(0, K)
            hbeA <- rep(0, K)
            hebA <- rep(0, K)
            hbpA <- rep(0, K)
            hpbA <- rep(0, K)
            hepA <- rep(0, K)
            hpeA <- rep(0, K)
            he1e2A <- rep(0, K)
            he2e1A <- rep(0, K)                    
          }
        for(tt in 1:max(dat$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- dat$tmar-dat$bd+tt-1 >= dat$tms & dat$tmar-dat$bd+tt-1 < dat$tmar
            f.sex <- dat$tmar-dat$bd+tt-1 >= dat$tfs & dat$tmar-dat$bd+tt-1 < dat$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- rep(0,K)
            p.f.bef <- rep(0,K)    
            p.m.bef[m.sex] <- (1 - exp(-bmb * epicf[cbind(dat$tmar[m.sex]-dat$bd[m.sex]+tt-1, dat$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * epicm[cbind(dat$tmar[f.sex]-dat$bd[f.sex]+tt-1, dat$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            if(survive)
              {
                p.m.bef.a <- rep(0,K)
                p.f.bef.a <- rep(0,K)
                ## csurv[time til interview, age in months in this month]
                p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(dat$mage[m.sex]-dat$cd[m.sex]-dat$bd[m.sex]+tt-1, dat$cd[m.sex]+dat$bd[m.sex]-tt+1)]
                p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(dat$fage[f.sex]-dat$cd[f.sex]-dat$bd[f.sex]+tt-1, dat$cd[f.sex]+dat$bd[f.sex]-tt+1)]
              }
            ## iterate probabilities based on previous values for only cases where it needs updating
            s..[e.sex] <- s..L[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mb.[e.sex] <- mb.L[e.sex]*(1 - p.f.bef[e.sex]) + s..L[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            f.b[e.sex] <- f.bL[e.sex]*(1 - p.m.bef[e.sex]) + s..L[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            ## for individuals infected in the same month, assign
            ## the order of infection based on competing risks
            ## formula, but if the denominator is 0, replace both
            ## with 0 to avoid errors.
            p.mfirst <- p.m.bef[e.sex] / (p.m.bef[e.sex]+p.f.bef[e.sex])
            p.ffirst <- 1-p.mfirst
            p.mfirst[is.na(p.mfirst)] <- 0
            p.ffirst[is.na(p.ffirst)] <- 0                
            hb1b2[e.sex] <- hb1b2L[e.sex] + p.mfirst * s..L[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
                                            mb.L[e.sex]*p.f.bef[e.sex]
            hb2b1[e.sex] <- hb2b1L[e.sex] + p.ffirst * s..L[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
                                            f.bL[e.sex]*p.m.bef[e.sex]
            ## iterate joint probabilities with survival
            if(survive)
              {
                mb.A[e.sex] <- mb.AL[e.sex]*(1 - p.f.bef[e.sex]) + s..L[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
                f.bA[e.sex] <- f.bAL[e.sex]*(1 - p.m.bef[e.sex]) + s..L[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
                ## for individuals infected in the same month, assign
                ## the order of infection based on competing risks
                ## formula, but if the denominator is 0, replace both
                ## with 0 to avoid errors.
                p.mfirst.a <- p.m.bef.a[e.sex] / (p.m.bef.a[e.sex]+p.f.bef.a[e.sex])
                p.ffirst.a <- 1-p.mfirst.a
                p.mfirst.a[is.na(p.mfirst.a)] <- 0
                p.ffirst.a[is.na(p.ffirst.a)] <- 0                
                hb1b2A[e.sex] <- hb1b2AL[e.sex] + p.mfirst.a * s..L[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
                                                  mb.AL[e.sex] * p.f.bef.a[e.sex]
                hb2b1A[e.sex] <- hb2b1AL[e.sex] + p.ffirst.a * s..L[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
                                                  f.bAL[e.sex] * p.m.bef.a[e.sex]
                ## Update *A*live *L*ast states
                mb.AL[e.sex] <- mb.A[e.sex]
                f.bAL[e.sex] <- f.bA[e.sex]
                hb1b2AL[e.sex] <- hb1b2A[e.sex]
                hb2b1AL[e.sex] <- hb2b1A[e.sex]                
              }
            ## Update other *L*ast states
            s..L[e.sex] <- s..[e.sex]
            mb.L[e.sex] <- mb.[e.sex]
            f.bL[e.sex] <- f.b[e.sex]
            hb1b2L[e.sex] <- hb1b2[e.sex]
            hb2b1L[e.sex] <- hb2b1[e.sex]                        
          }
        ## probability of being infected by partner (constant, used inside loop)
        p.m.part <- 1 - exp(-bmp)
        p.f.part <- 1 - exp(-bfp)
        ## Now loop through marriage
        for(tt in 1:max(dat$cd-1))
          {
            ## are partners formed in a couple?
            fmd <- dat$cd >= tt              
            ######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme*epicf[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe*epicm[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            if(survive)
              {
                ## Survival probabilities
                s.p.m <- csurv[cbind(dat$mage[fmd]-dat$cd[fmd]+tt-1, dat$tint[fmd] - (dat$tmar[fmd] + tt - 1))]
                s.p.f <- csurv[cbind(dat$fage[fmd]-dat$cd[fmd]+tt-1, dat$tint[fmd] - (dat$tmar[fmd] + tt - 1))]
                ## Transmission probabilities from partner (jointly with survival)
                p.m.part.a <- p.m.part * s.p.m
                p.f.part.a <- p.f.part * s.p.f
                p.m.exc.a <- p.m.exc * s.p.m
                p.f.exc.a <- p.f.exc * s.p.f
              }
            ######################################################################
            ## iterate probabilities
            s..[fmd] <- s..L[fmd]*(1-p.m.exc)*(1-p.f.exc)
            mb.[fmd] <- mb.L[fmd]*(1-p.f.exc)*(1-p.f.part)
            me.[fmd] <- me.L[fmd]*(1-p.f.exc)*(1-p.f.part) + s..L[fmd]*p.m.exc*(1-p.f.exc)
            f.b[fmd] <- f.bL[fmd]*(1-p.m.exc)*(1-p.m.part)
            f.e[fmd] <- f.eL[fmd]*(1-p.m.exc)*(1-p.m.part) + s..L[fmd]*p.f.exc*(1-p.m.exc)
##          hb1b2[fmd] <- hb1b2L[fmd] # Doesn't change during couple duration
##          hb2b1[fmd] <- hb2b1L[fmd] # Doesn't change during couple duration                
            hbe[fmd] <- hbeL[fmd] + mb.L[fmd]*(1-p.f.part)*p.f.exc
            heb[fmd] <- hebL[fmd] + f.bL[fmd]*(1-p.m.part)*p.m.exc
            hbp[fmd] <- hbpL[fmd] + mb.L[fmd]*p.f.part
            hpb[fmd] <- hpbL[fmd] + f.bL[fmd]*p.m.part
            hep[fmd] <- hepL[fmd] + me.L[fmd]*p.f.part
            hpe[fmd] <- hpeL[fmd] + f.eL[fmd]*p.m.part
            ## for individuals infected in the same month, assign
            ## the order of infection based on competing risks
            ## formula, but if the denominator is 0, replace both
            ## with 0 to avoid errors.
            p.mfirst <- p.m.exc / (p.m.exc+p.f.exc)
            p.ffirst <- 1-p.mfirst
            p.mfirst[is.na(p.mfirst)] <- 0
            p.ffirst[is.na(p.ffirst)] <- 0                
            he1e2[fmd] <- he1e2L[fmd] + p.mfirst * s..L[fmd]*p.m.exc*p.f.exc +
                                        me.L[fmd]*(1-p.f.part)*p.f.exc
            he2e1[fmd] <- he2e1L[fmd] + p.ffirst * s..L[fmd]*p.m.exc*p.f.exc +
                                        f.eL[fmd]*(1-p.m.part)*p.m.exc
            ######################################################################
            ## Iterate probabilities jointly with survival until survey.
            ## Note for probabilities of not being infected, we don't
            ## use the joint probability with being alive at sampling.
            if(survive)
              {
                mb.A[fmd] <- mb.AL[fmd]*(1-p.f.exc)*(1-p.f.part)
                me.A[fmd] <- me.AL[fmd]*(1-p.f.exc)*(1-p.f.part) + s..L[fmd]*p.m.exc.a*(1-p.f.exc)
                f.bA[fmd] <- f.bAL[fmd]*(1-p.m.exc)*(1-p.m.part)
                f.eA[fmd] <- f.eAL[fmd]*(1-p.m.exc)*(1-p.m.part) + s..L[fmd]*p.f.exc.a*(1-p.m.exc)
##              hb1b2A[fmd] <- hb1b2AL[fmd] # Doesn't change during couple duration
##              hb2b1A[fmd] <- hb2b1AL[fmd] # Doesn't change during couple duration                
                hbeA[fmd] <- hbeAL[fmd] + mb.AL[fmd]*(1-p.f.part)*p.f.exc.a
                hebA[fmd] <- hebAL[fmd] + f.bAL[fmd]*(1-p.m.part)*p.m.exc.a
                hbpA[fmd] <- hbpAL[fmd] + mb.AL[fmd]*p.f.part.a
                hpbA[fmd] <- hpbAL[fmd] + f.bAL[fmd]*p.m.part.a
                hepA[fmd] <- hepAL[fmd] + me.AL[fmd]*p.f.part.a
                hpeA[fmd] <- hpeAL[fmd] + f.eAL[fmd]*p.m.part.a
                ## for individuals infected in the same month, assign
                ## the order of infection based on competing risks
                ## formula, but if the denominator is 0, replace both
                ## with 0 to avoid errors.
                p.mfirst.a <- p.m.exc.a / (p.m.exc.a+p.f.exc.a)
                p.ffirst.a <- 1-p.mfirst.a
                p.mfirst.a[is.na(p.mfirst.a)] <- 0
                p.ffirst.a[is.na(p.ffirst.a)] <- 0                
                he1e2A[fmd] <- he1e2AL[fmd] + p.mfirst.a * s..L[fmd]*p.m.exc.a*p.f.exc.a +
                                          me.AL[fmd]*(1-p.f.part)*p.f.exc.a
                he2e1A[fmd] <- he2e1AL[fmd] + p.ffirst.a * s..L[fmd]*p.m.exc.a*p.f.exc.a +
                                          f.eAL[fmd]*(1-p.m.part)*p.m.exc.a
                ## update *L*ast month states for *A*live states
                mb.AL[fmd] <- mb.A[fmd]
                me.AL[fmd] <- me.A[fmd]
                f.bAL[fmd] <- f.bA[fmd]
                f.eAL[fmd] <- f.eA[fmd]
##              hb1b2AL[fmd] <- hb1b2A[fmd]
##              hb2b1AL[fmd] <- hb2b1A[fmd]                
                hbeAL[fmd] <- hbeA[fmd]
                hebAL[fmd] <- hebA[fmd]
                hbpAL[fmd] <- hbpA[fmd]
                hpbAL[fmd] <- hpbA[fmd]
                hepAL[fmd] <- hepA[fmd]
                hpeAL[fmd] <- hpeA[fmd]              
                he1e2AL[fmd] <- he1e2A[fmd]
                he2e1AL[fmd] <- he2e1A[fmd]                              
              }
            ## update other *L*ast month states
            s..L[fmd] <-  s..[fmd]
            mb.L[fmd] <-  mb.[fmd]
            me.L[fmd] <-  me.[fmd]
            f.bL[fmd] <-  f.b[fmd]
            f.eL[fmd] <-  f.e[fmd]
##          hb1b2L[fmd] <-  hb1b2[fmd]
##          hb2b1L[fmd] <-  hb2b1[fmd]            
            hbeL[fmd] <-  hbe[fmd]
            hebL[fmd] <-  heb[fmd]
            hbpL[fmd] <-  hbp[fmd]
            hpbL[fmd] <-  hpb[fmd]
            hepL[fmd] <-  hep[fmd]
            hpeL[fmd] <-  hpe[fmd]         
            he1e2L[fmd] <-  he1e2[fmd]
            he2e1L[fmd] <-  he2e1[fmd]                     
          }
        allstates <- data.frame(s..,mb.,me.,f.b,f.e,hb1b2,hb2b1,hbe,heb,hbp,hpb,hep,hpe,he1e2,he2e1,
                                mb.A,me.A,f.bA,f.eA,hb1b2A,hb2b1A,hbeA,hebA,hbpA,hpbA,hepA,hpeA,he1e2A,he2e1A)
        ss <- s..
        mm <- mb. + me.
        ff <- f.b + f.e
        hh <- hb1b2 + hb2b1 + hbe + heb + hbp + hpb + hep + hpe + he1e2 + he2e1
        ## Calculate probability of data given parameters * priors of paramters
        if(survive)
          {
            mmA <- mb.A + me.A
            ffA <- f.bA + f.eA
            hhA <- hb1b2A + hb2b1A + hbeA + hebA + hbpA + hpbA + hepA + hpeA + he1e2A + he2e1A
            pser.a <- cbind(hhA, mmA, ffA, ss)
          }
        pser <- cbind(hh, mm, ff, ss)
        if(trace & !sim)
          {
            ######################################################################
            ## Route of transmission breakdowns for observed couples
            ## (conditional on survival)
            ######################################################################
            ## male breakdown amongst observed M+F- couples, partner *N*egative
            pibNA <- sum(mb.A[mm.log] / mmA[mm.log]) / sum(mm.log)
            pieNA <- sum(me.A[mm.log] / mmA[mm.log]) / sum(mm.log)
            ## female breakdown amongst observed M-F+ couples, partner *N*egative
            piNbA <- sum(f.bA[ff.log] / ffA[ff.log]) / sum(ff.log)
            piNeA <- sum(f.eA[ff.log] / ffA[ff.log]) / sum(ff.log)
            ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
            pibPA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hbeA[hh.log] + hbpA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piePA <- sum((hebA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
            pipPA <- sum((hpbA[hh.log] + hpeA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
            piPbA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hebA[hh.log] + hpbA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPeA <- sum((hbeA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPpA <- sum((hbpA[hh.log] + hepA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
            pibUA <- (pibNA*sum(mm.log) + pibPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pieUA <- (pieNA*sum(mm.log) + piePA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pipUA <- (pipPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
            piUbA <- (piNbA*sum(mm.log) + piPbA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            piUeA <- (piNeA*sum(mm.log) + piPeA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            piUpA <- (piPpA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            ######################################################################
            ## Give pieNA, piNeA, piePA, piPeA, for each couple
            if(give.pis)
              {
                ## probability infection was extracouple given ser
                piCe.A <- rep(NA, K)
                piC.eA <- rep(NA, K)
                piCe.A[mm.log] <- me.A[mm.log] / mmA[mm.log]
                piCe.A[hh.log] <- (hebA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
                piC.eA[ff.log] <- f.eA[ff.log] / ffA[ff.log]
                piC.eA[hh.log] <- (hbeA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
                pis <- data.frame(piCe.A, piC.eA)
              }
            ######################################################################
            ## Route of transmission breakdowns for inferred
            ## pseudopopulation (unconditional on survival)
            ######################################################################
            ######################################################################
            ## Index infections, do with estimators summing over all
            ## infected couples and over all couples
            ## version 1 - all infected couples
            mb1. <- rowSums(allstates[!ss.log, c("mb.","hb1b2","hbe","hbp")])
            me1. <- rowSums(allstates[!ss.log, c("me.","he1e2","hep")])
            ## female
            f.b1 <- rowSums(allstates[!ss.log, c("f.b","hb2b1","heb","hpb")])
            f.e1 <- rowSums(allstates[!ss.log, c("f.e","he2e1","hpe")])
            all.infA <- rowSums(allstates[!ss.log,names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))]])
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
            piGb1.sumI <- mb1.Infl / IndInfl
            piGe1.sumI <- me1.Infl / IndInfl
            piG.b1sumI <- f.b1Infl / IndInfl
            piG.e1sumI <- f.e1Infl / IndInfl
            ## ## version 2 - sum over all couples
            mb1. <- rowSums(allstates[, c("mb.","hb1b2","hbe","hbp")])
            me1. <- rowSums(allstates[, c("me.","he1e2","hep")])
            ## female
            f.b1 <- rowSums(allstates[, c("f.b","hb2b1","heb","hpb")])
            f.e1 <- rowSums(allstates[, c("f.e","he2e1","hpe")])
            all.infA <- rowSums(allstates[,c("s..", names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))])])
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
            piGb1.sumIS <- mb1.Infl / IndInfl
            piGe1.sumIS <- me1.Infl / IndInfl
            piG.b1sumIS <- f.b1Infl / IndInfl
            piG.e1sumIS <- f.e1Infl / IndInfl
            ## put them all in a dataframe
            pop.avs <- data.frame(
                                  ## conditional on survival
                                  pibNA, pieNA, # b/e in +- given A
                                  piNbA, piNeA, # b/e in -+ given A
                                  pibPA, piePA, pipPA, # b/e/p in male in ++ given A
                                  piPbA, piPeA, piPpA, # b/e/p in female in ++ given A
                                  pibUA, pieUA, pipUA, # b/e/p in male in any given A
                                  piUbA, piUeA, piUpA, # b/e/p in female in any given A
                                  ## unconditional on survival version 1
                                  piGb1.sumI, piGe1.sumI, # b/e index in males amongst all infected 
                                  piG.b1sumI, piG.e1sumI, # b/e index in females amongst all infected
                                  piGb1.sumIS, piGe1.sumIS, # b/e index in males amongst all infected 
                                  piG.b1sumIS, piG.e1sumIS) # b/e index in females amongst all infected
            ######################################################################
            ## Project incidence forward 12 months for each of the 3
            ## couple types (ss, mh, fh) for each country in the data
            ## set (because they have different population prevalences
            num.country <- length(unique(dat$epic.ind))
            cc.inds <- unique(dat$epic.ind)
            ## concordant negative
            ss12.ssL <- rep(1, num.country)
            mm12.ssL <- rep(0, num.country)
            ff12.ssL <- rep(0, num.country)
            hh12.ssL <- rep(0, num.country)            
            ## male positive discordant
            mm12.mmL <- rep(1, num.country)
            hh12.mmL <- rep(0, num.country)            
            ## female positive discordant
            ff12.ffL <- rep(1, num.country)
            hh12.ffL <- rep(0, num.country)
            ## initialize pis
            pi.m.part12.ss <- 0
            pi.f.part12.ss <- 0
            pi.m.exc12.ss <- 0
            pi.f.exc12.ss <- 0
            pi.f.part12.mm <- 0
            pi.f.exc12.mm <- 0            
            pi.m.part12.ff <- 0
            pi.m.exc12.ff <- 0            
            for(tt in 1:12)
              {
                ######################################################################
                ## Transmission probabilities
                ## probability infected extracouply in various months of 2011
                p.m.exc <- 1 - exp(-bme*epicf[1332+tt-1, cc.inds])
                p.f.exc <- 1 - exp(-bfe*epicm[1332+tt-1, cc.inds])
                ## concordant negative couples
                ss12.ss <- ss12.ssL*(1-p.m.exc)*(1-p.f.exc)
                mm12.ss <- mm12.ssL*(1-p.f.exc)*(1-p.f.part) + ss12.ssL*p.m.exc*(1-p.f.exc)
                ff12.ss <- ff12.ssL*(1-p.m.exc)*(1-p.m.part) + ss12.ssL*p.f.exc*(1-p.m.exc)
                hh12.ss <- hh12.ssL + ss12.ssL* p.m.exc*p.f.exc +
                    mm12.ssL*(p.f.part + (1-p.f.part)*p.f.exc) +
                    ff12.ssL*(p.m.part + (1-p.m.part)*p.m.exc)
                pi.m.part12.ss <- pi.m.part12.ss + ff12.ssL*p.m.part
                pi.f.part12.ss <- pi.f.part12.ss + mm12.ssL*p.f.part        
                pi.m.exc12.ss <- pi.m.exc12.ss + (ss12.ssL + ff12.ssL*(1-p.m.part))*p.m.exc
                pi.f.exc12.ss <- pi.f.exc12.ss + (ss12.ssL + mm12.ssL*(1-p.f.part))*p.f.exc
                ## male positive couples & female seroconversion
                mm12.mm <- mm12.mmL*(1-p.f.exc)*(1-p.f.part)
                hh12.mm <- hh12.mmL + mm12.mmL*(p.f.part + (1-p.f.part)*p.f.exc) 
                pi.f.part12.mm <- pi.f.part12.mm + mm12.mmL*p.f.part        
                pi.f.exc12.mm <- pi.f.exc12.mm + mm12.mmL*(1-p.f.part)*p.f.exc
                ## female positive couples & male seroconversion                  
                ff12.ff <- ff12.ffL*(1-p.m.exc)*(1-p.m.part)
                hh12.ff <- hh12.ffL + ff12.ffL*(p.m.part + (1-p.m.part)*p.m.exc)
                pi.m.part12.ff <- pi.m.part12.ff + ff12.ffL*p.m.part
                pi.m.exc12.ff <- pi.m.exc12.ff + ff12.ffL*(1-p.m.part)*p.m.exc
                ss12.ssL <- ss12.ss
                mm12.ssL <- mm12.ss
                ff12.ssL <- ff12.ss
                hh12.ssL <- hh12.ss
                ## male positive discordant
                mm12.mmL <- mm12.mm
                hh12.mmL <- hh12.mm
                ## female positive discordant
                ff12.ffL <- ff12.ff
                hh12.ffL <- hh12.ff
              }
            n.m.part.dc <- 0
            n.m.part.cc <- 0
            n.f.part.dc <- 0
            n.f.part.cc <- 0
            n.m.exc.dc <- 0
            n.m.exc.cc <- 0
            n.f.exc.dc <- 0
            n.f.exc.cc <- 0
            ## add all the different countries incidence by scaling by serotype
            for(cc in 1:num.country)
              {
                n.m.part.dc <- n.m.part.dc + pi.m.part12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
                n.f.part.dc <- n.f.part.dc + pi.f.part12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
                n.m.part.cc <- n.m.part.cc + pi.m.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.part.cc <- n.f.part.cc + pi.f.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
                n.m.exc.dc <- n.m.exc.dc + pi.m.exc12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
                n.f.exc.dc <- n.f.exc.dc + pi.f.exc12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
                n.m.exc.cc <- n.m.exc.cc + pi.m.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.exc.cc <- n.f.exc.cc + pi.f.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
              }
            n.m.part.tot <- n.m.part.dc + n.m.part.cc
            n.f.part.tot <- n.f.part.dc + n.f.part.cc            
            n.m.exc.tot <- n.m.exc.dc + n.m.exc.cc
            n.f.exc.tot <- n.f.exc.dc + n.f.exc.cc
            proj12 <- data.frame(n.m.part.dc, n.f.part.dc, # incidence per 1000
                                 n.m.part.cc, n.f.part.cc,
                                 n.m.exc.dc, n.f.exc.dc,
                                 n.m.exc.cc, n.f.exc.cc,
                                 n.m.part.tot, n.f.part.tot,
                                 n.m.exc.tot, n.f.exc.tot) / sum(!hh.log) * 1000
            prop.exc.m <- n.m.exc.tot / (n.m.exc.tot + n.m.part.tot)
            prop.exc.f <- n.f.exc.tot / (n.f.exc.tot + n.f.part.tot)
            proj12 <- data.frame(proj12, prop.exc.m, prop.exc.f)            
            ## relative rate of transmission coefficient extracouply vs before relationship
            rr.m.out <- bme/bmb
            rr.f.out <- bfe/bfb
            rr.m.in <- bmp/bme
            rr.f.in <- bfp/bfe
            rr.m.pbef <- bmp/bmb        #partner to before
            rr.f.pbef <- bfp/bfb
            ## relative rate of transmission coefficient extracouply and before relationship between males and females
            rr.mf.bef <- bmb/bfb
            rr.mf.exc <- bme/bfe
            ## rho is the last one
            ## relative rate of contact/risk paramter (i.e. accounting
            ## for difference in per coital act probability as estimated
            ## from within partnership transmission.
            rr.mf.bef.cont <- rr.mf.bef * rho
            rr.mf.exc.cont <- rr.mf.exc * rho
            rrs <- data.frame(rr.m.out = rr.m.out, rr.f.out = rr.f.out,
                              rr.m.in = rr.m.in, rr.f.in = rr.f.in,
                              rr.m.pbef = rr.m.pbef, rr.f.pbef = rr.f.pbef,
                              rr.mf.bef = rr.mf.bef, rr.mf.exc = rr.mf.exc,
                              rr.mf.bef.cont = rr.mf.bef.cont, rr.mf.exc.cont = rr.mf.exc.cont)
          }
        if(sim) # if simulating data
          {
            probs <- NA
            lprob <- NA
            ## create couple state probability *A*live & *D*ead
            sim.probs <- data.frame(s..A = s..,
                                    mb.A,   mb.D   =  mb.- mb.A,
                                    me.A,   me.D   =  me. - me.A,
                                    f.bA,   f.bD   =  f.b - f.bA,
                                    f.eA,   f.eD   =  f.e - f.eA,
                                    hb1b2A, hb1b2D = hb1b2 - hb1b2A, # note some of the dead cases were infected by dead partners, so can only use the index cases in any h couple
                                    hb2b1A, hb2b1D = hb2b1 - hb2b1A,
                                    hbeA,   hbeD   =  hbe - hbeA,
                                    hebA,   hebD   =  heb - hebA,
                                    hepA,   hepD   =  hep - hepA,
                                    hpeA,   hpeD   =  hpe - hpeA,
                                    hbpA,   hbpD   =  hbp - hbpA,
                                    hpbA,   hpbD   =  hpb - hpbA,
                                    he1e2A, he1e2D = he1e2 - he1e2A,
                                    he2e1A, he2e1D = he2e1 - he2e1A)
            for(ii in 1:nrow(dat))
              {
                        dat$cat[ii] <- which(rmultinom(1, 1, sim.probs[ii,])==1)
              }
            dat$cat.nm <- names(sim.probs)[dat$cat]
            dat$cat.nm <- factor(dat$cat.nm, levels = names(sim.probs))
            K <- nrow(dat)
            if(!survive) pser.a <- NA
          }else{ ## if not simulating data calculate likelihood p(data|pars)
            if(survive)
              {                         # must NORMALIZE probabilities to 1 for likelihood!
                probs <- pser.a[cbind(1:K,dat$ser)] / rowSums(pser.a) # accounting for survival
              }else{
                probs <- pser[cbind(1:K,dat$ser)] # if not accounting for survival
                pser.a <- NA
              }
            if(sum(probs==0)==0) # if non of the serotatuses occur with 0 probability in the current model
              {
                lprob <- sum(log(probs)) + dnorm(log(rho), log(trans.ratio), 1/2, log = T)
                if(length(compars)>0) clprob <- sum(log(cprobs)) + dnorm(as.numeric(compars["lrho"]),
                                                                         log(trans.ratio), 1/2, log = T)
              }else{ # if some of the serostatuses are 0, then the current parameters have 0 probability
                lprob <- -Inf
              }
          }
      }
    if(length(compars)==0)
      {
        clprob <- NA
        cprobs <- NA
      }
    if(sim)
      {
        if(trace)
          {
            if(give.pis)
              {
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs,  proj12=proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser, dat = dat, sim.probs, allstates = allstates,
                        clprob = clprob, probs = probs, cprobs = cprobs))
          }
      }else{                            # if not simulating
        if(trace)
          {
            if(give.pis)
              {
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs,  proj12=proj12, pis = pis, allstates = allstates,
                            pser.a = pser.a, pser = pser, probs = probs))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12,
                            pser.a = pser.a, pser = pser, probs = probs))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser))
          }
      }
  }



## ## Simulate with state probability model
## sim <- pcalc(pars, dat, sim = T, survive = T, browse = F, trace = F)$allstates

## out <- sim.fn(pars, dat, browse=F, vfreq = 500)
## ## real serostatus data
## print("real dat serostatus breakdown")
## round(xtabs(~ser,dat) / nrow(dat), 3)
## print("simulated dat serostatus breakdown")
## round(xtabs(~ser, out, subset = alive) / sum(out$alive), 3)


## ## Compare event-driven and multinomial models
## stuff <- rbind(colMeans(sim), transf(out))
## rownames(stuff) <- c("multinomial state prob", "stochastic event-driven")
## stuff[,c(1,16:29)]*nrow(dat)


## stuff[2,c(1,16:29)] / stuff[1,c(1,16:29)]

## ## simulated serostatus data
## xtabs(~ser,out)
## xtabs(~ser + alive,out)



## xtabs(alive~ mser + fser, aggregate(alive~ mser + fser, out, mean))

## ## Look at causes of death
## tabm <- xtabs(~mcoi + alive, out)
## tabf <- xtabs(~fcoi + alive, out)
## tabm[,1]/rowSums(tabm)
## tabf[,1]/rowSums(tabf)
## ## shh <- c("tms", "tfs", "tmar", "tint","mds", "fds","mardur.mon","mcoi","fcoi","mdod","fdod","alive")
## ## tail(out[,shh])



## layout(matrix(1:6, 3,2))
## par(oma = c(0,0,2,0))
## breaks <- 0:50
## mains <- c("precouple","within-couple","extracouple")
## hist((dat$mdod-dat$mdoi)[dat$mcoi =="b"]/12, col = "black", breaks = breaks, main = mains[1], xlab = "survival time")
## hist((dat$mdod-dat$mdoi)[dat$mcoi =="p"]/12, col = "black", breaks = breaks, main = mains[2], xlab = "survival time")
## hist((dat$mdod-dat$mdoi)[dat$mcoi =="e"]/12, col = "black", breaks = breaks, main = mains[3], xlab = "survival time")
## hist((dat$fdod-dat$fdoi)[dat$fcoi =="b"]/12, col = "black", breaks = breaks, main = mains[1], xlab = "survival time")
## hist((dat$fdod-dat$fdoi)[dat$fcoi =="p"]/12, col = "black", breaks = breaks, main = mains[2], xlab = "survival time")
## hist((dat$fdod-dat$fdoi)[dat$fcoi =="e"]/12, col = "black", breaks = breaks, main = mains[3], xlab = "survival time")
## mtext("males", side = 3, outer = T, line = 0, adj = .3)
## mtext("females", side = 3, outer = T, line = 0, adj = .8)
## mtext("survival time distributions", side = 3, outer = T, line = 2, adj = .5)
