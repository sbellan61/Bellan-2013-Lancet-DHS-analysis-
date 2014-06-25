pcalc <- function(pars, dat, browse = F, compars = NULL, # compare lprob for true pars with other pars (fitted)
                  give.pis = F,              # return individual pi values for couples (for outside mcmc)
                  sim = F,                   # just use this to simulate data given parameters? (outputs pser.a)
                  spar.log = F,                  # use sexual partner acquisition rate
                  heter = F,                     # individual heterogeneity; simulation only
                  ## using global variables for these
#                  partner.arv = F,               # ART coverage affects within-partnership transmission?
#                  low.coverage.arv = F,          # if T, only 50% of those on ART are not infectious
                  survive = T,               # account for survival in analysis
                  cond.sim = F,              # only simulate individuals that will live
                  lrho.sd = 1/2,             # sd on lrho prior
                  trace = T) # only do certain calculations when tracing parameters (i.e. for non-thinned versions)
  {
    if(browse) browser()
    if(low.coverage.arv)    # if assuming only 50% of those on ART are not-infectious
      {
        cov.scalar <- .5
      }else{
        cov.scalar <- 1
      }
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
        bmb <- as.numeric(pars["bmb"])
        bfb <- as.numeric(pars["bfb"])
        bme <- as.numeric(pars["bme"])
        bfe <- as.numeric(pars["bfe"])
        bmp <- as.numeric(pars["bmp"])
        rho <- exp(as.numeric(pars["lrho"])) # feeding in log(rho)
        bfp <- bmp * rho
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
            ## spar
            if(spar.log)
              {
                m.spar <- rep(1, K)
                f.spar <- rep(1, K)
              }else{
                m.spar <- rep(1, K)
                f.spar <- rep(1, K)
              }
            if(heter)
              {
                m.het <- exp( rnorm(K, mean = 0, sd = 1) )
                f.het <- exp( rnorm(K, mean = 0, sd = 1) )
              }else{
                m.het <- rep(1, K)
                f.het <- rep(1, K)
              }
            p.m.bef[m.sex] <- (1 - exp(-bmb * m.spar[m.sex] * m.het[m.sex] *epicf[cbind(dat$tmar[m.sex]-dat$bd[m.sex]+tt-1, dat$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * f.spar[f.sex] * f.het[f.sex] *epicm[cbind(dat$tmar[f.sex]-dat$bd[f.sex]+tt-1, dat$epic.ind[f.sex])]))
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
            p.m.exc <- (1 - exp(-bme * m.spar[fmd] * m.het[fmd] *epicf[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe * f.spar[fmd] * f.het[fmd] *epicm[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            ## adjust probability of being infected by partner by probability partner is infectious (i.e. not on ART)
            if(partner.arv)
              {
                p.m.part <- 1 - exp(-bmp * (1 - cov.scalar * art.prev[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
                p.f.part <- 1 - exp(-bfp * (1 - cov.scalar * art.prev[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
              }
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
        if(trace & sim) # calculate expected route of transmission breakdowns for couples with unknown (or simulated serostatus)
          {
######################################################################
            ## Route of transmission breakdowns for observed couples
            ## (conditional on survival)
######################################################################
            ## male breakdown amongst observed M+F- couples, partner *N*egative
            pibNA <- mean(mb.A / mmA, na.rm = T) 
            pieNA <- mean(me.A / mmA, na.rm = T) 
            ## female breakdown amongst observed M-F+ couples, partner *N*egative
            piNbA <- mean(f.bA / ffA, na.rm = T)
            piNeA <- mean(f.eA / ffA, na.rm = T)
            ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
            pibPA <- mean((hb1b2A + hb2b1A + hbeA + hbpA) / hhA, na.rm = T) 
            piePA <- mean((hebA + hepA + he1e2A + he2e1A) / hhA, na.rm = T) 
            pipPA <- mean((hpbA + hpeA) / hhA, na.rm = T)                  
            ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
            piPbA <- mean((hb1b2A + hb2b1A + hebA + hpbA) / hhA, na.rm = T) 
            piPeA <- mean((hbeA + hpeA + he1e2A + he2e1A) / hhA, na.rm = T) 
            piPpA <- mean((hbpA + hepA) / hhA, na.rm = T)                  
            ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
            pibUA <- mean((mb.A + hb1b2A + hb2b1A + hbeA + hbpA) / (mmA + hhA), na.rm = T)
            pieUA <- mean((me.A + hebA + hepA + he1e2A + he2e1A) / (mmA + hhA), na.rm = T)
            pipUA <- mean((hpbA + hpeA) / (mmA + hhA), na.rm = T)
            ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
            piUbA <- mean((f.bA + hb1b2A + hb2b1A + hebA + hpbA) / (ffA + hhA), na.rm = T)
            piUeA <- mean((f.eA + hbeA + hpeA + he1e2A + he2e1A) / (ffA + hhA), na.rm = T)
            piUpA <- mean((hbpA + hepA) / (ffA + hhA), na.rm = T)
######################################################################
            pop.avs <- data.frame(
                                  ## conditional on survival
                                  pibNA, pieNA, # b/e in +- given A
                                  piNbA, piNeA, # b/e in -+ given A
                                  pibPA, piePA, pipPA, # b/e/p in male in ++ given A
                                  piPbA, piPeA, piPpA, # b/e/p in female in ++ given A
                                  pibUA, pieUA, pipUA, # b/e/p in male in any given A
                                  piUbA, piUeA, piUpA) # b/e/p in female in any given A
          }
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
            piUbA <- (piNbA*sum(ff.log) + piPbA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
            piUeA <- (piNeA*sum(ff.log) + piPeA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
            piUpA <- (piPpA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
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
                if(partner.arv)
                  {
                    p.m.part <- 1 - exp(-bmp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds]))
                    p.f.part <- 1 - exp(-bfp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds]))
                  }
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
            n.m.dc <- n.m.part.dc + n.m.exc.dc
            n.f.dc <- n.f.part.dc + n.f.exc.dc
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
            prop.exc.m.dc <- n.m.exc.dc / n.m.dc
            prop.exc.f.dc <- n.f.exc.dc / n.f.dc
            proj12 <- data.frame(proj12, prop.exc.m, prop.exc.f, prop.exc.m.dc, prop.exc.f.dc) 
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
                lprob <- sum(log(probs)) + dnorm(log(rho), log(trans.ratio), lrho.sd, log = T)
                if(length(compars)>0) clprob <- sum(log(cprobs)) + dnorm(as.numeric(compars["lrho"]),
                                                                         log(trans.ratio), lrho.sd, log = T)
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
                return(list(lprob = lprob,pop.avs = pop.avs, proj12 = proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser, dat = dat, sim.probs, allstates = allstates,
                        clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
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




init.fxn <- function(seed = 1)
  {
    set.seed(seed)
    ## sample uniform dispersed on log scale
    lpars <- c(bmb = runif(1, -6, -3),
               bfb = runif(1, -6, -3),
               bme = runif(1, -6, -3),
               bfe = runif(1, -6, -3),
               bmp = runif(1, -6, -3))
    pars <- c(exp(lpars), lrho = runif(1, -1, 1))
    return(pars)
  }

    
## MCMC SAMPLER
sampler <- function(sd.props = sd.props, inits,
                    multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
                    verbose = T, tell = 100, seed = 1, lrho.sd,
                    niter = 6*1000, survive,
                    nthin = 5,
                    nburn = 1000, browse=F)
  {
    if(browse)  browser()
    set.seed(seed)
    pars <- inits
    vv <- 2
    accept <- 0                  #track each parameters acceptance individually
    cur <- pcalc(pars, dat = dat, trace = T, survive = survive)     #calculate first log probability
    lprob.cur <- cur$lprob
    out <- t(data.frame(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"])),
                          cur$pop.avs, cur$rrs, cur$proj12)))
    last.it <- 0
    start <- Sys.time()
    while(vv < niter + 1)
      {
        if(verbose & vv%%tell+1==1) print(paste("on iteration",vv,"of",last.it + niter + 1))
        pars.prop <- pars              #initialize proposal parameterr vector
        ## propose new parameter vector
        if(multiv)
          {
            pars.prop <- pars.prop + rmnorm(1, mean = 0, varcov = covar)
            pars.prop <- as.vector(pars.prop) #otherwise is a matrix
            names(pars.prop) <- parnames
          }else{
            pars.prop <- pars.prop + rnorm(length(pars), mean = 0, sd = sd.props)
          }
        ## trace = T if in non-thinned iteration, or the previous one (in case of rejection)
        ## calculate proposal par log probability
        prop <- pcalc(pars.prop, dat = dat, survive = survive, trace = vv%%nthin + 1 %in% c(nthin,1))
        lprob.prop <- prop$lprob
        lmh <- lprob.prop - lprob.cur       # log Metropolis-Hastings ratio
        ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
        if(lmh >= 0 | runif(1,0,1) <= exp(lmh))
          {
            pars <- pars.prop
            if(vv>nburn) accept <- accept + 1 #only track acceptance after burn-in
            lprob.cur <- lprob.prop
            cur <- prop
          }
        if(vv%%nthin + 1 ==1)
          {
            out <- cbind(out,t(data.frame(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"])),
                                            cur$pop.avs, cur$rrs, cur$proj12))))
          }
        vv <- vv+1
      }
    if(verbose) print(paste("took", difftime(Sys.time(),start, units = "mins"),"mins"))
    aratio <- accept/((vv-nburn))    
    return(list(out = out[,1:ncol(out)>(nburn+1)/nthin], aratio = aratio, inits = inits))
  }



## plot trajectory for selected couples
ctraj <- function(pars, dat, browse = F, width = 11, height = 8.5, last.year = 2012,
                  plot.cpls = NULL, ylab = "probability in couple serostatus / pop prevalence",
                  surv = T,             # show curves joint with survival
                  nsurv = F,            # show curves not joint with survival
                  survive = T,          # show point on surv curve
                  dead = T,          # show curve for prob of couple dying
                  cprob = T,         # show conditional probabilities
                  cols = c("orange", "darkgreen", "purple","gray"),
                  col.dead = "black",
                  lty.surv = 1, lty.dead = 1,
                  lty.nsurv = 2, blob = T, blob.cex = 3, show.susc = T, cex.leg = 1,
                  show.age = F, do.leg = T, drp = 1, x.ticks = T, t.lwd = 1,
                  do.pdf = T, pdf.name = "surv traj.pdf")
  {
    if(browse) browser()
    K <- nrow(dat)
    if(do.pdf)
      {
      pdf(pdf.name, width = width, height = height)
      par(mar=c(0,4,.5,2))
    }
    for(cpl in plot.cpls)
      {
        temp <- dat[cpl,]
        xlim <- c(80*12, (last.year-1900)*12)
        if(show.susc)
          {
            ylim <- c(-.3,1)
          }else{
            ylim <- c(-.3,.5)
          }
        plot(0,0, type = "n", xlim = xlim, axes = F,
             ylim = ylim, xlab = "",  ylab = ylab, bty = "n")
        lines(xlim[1]:xlim[2], epicm[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 2, lwd = 3)
        lines(xlim[1]:xlim[2], epicf[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 3, lwd = 3)
        if(surv & !nsurv & dead & do.leg)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F- & alive", "dead before DHS sampling"),
                                                        col = c(cols,col.dead), lwd = 3,
                                                        lty = c(rep(lty.surv,3),1,lty.dead), bty = "n", cex = cex.leg)
        if(surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+","M+F+ & alive",
                                                    "M+F- & alive","M-F+ & alive","M-F-"),
                                       col = c(cols[1:3],cols), lwd = 3,
                                       lty = c(rep(lty.nsurv,3),rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(surv & !nsurv & do.leg & !dead)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F-"),
                                       col = cols, lwd = 3,
                                       lty = c(rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(!surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+", "M-F-"),
                                       col = cols, lwd = 3,
                                       lty = c(rep(lty.nsurv,3),1), bty = "n", cex = cex.leg)
        if(do.leg) legend("left", c("HIV population prev in M", "HIV population prev in F"),
               lty = 2:3, lwd = 3, col = "red", bty = "n", cex = cex.leg)
        if(x.ticks)
          {
            axis(1, seq(xlim[1],xlim[2], by = 60), seq(xlim[1],xlim[2], by = 60)/12 + 1900, pos = 0, las = 2)
          }else{
            axis(1, seq(xlim[1],xlim[2], by = 60), labels = NA, pos = 0, las = 2)
          }
        axis(2, seq(0,1,l=5), las = 2)
        bmb <- as.numeric(pars["bmb"])
        bfb <- as.numeric(pars["bfb"])
        bme <- as.numeric(pars["bme"])
        bfe <- as.numeric(pars["bfe"])
        bmp <- as.numeric(pars["bmp"])
        rho <- exp(as.numeric(pars["lrho"]))
        bfp <- bmp * rho
        out <- data.frame(cmc = temp$tmar - temp$bd, ss = 1, mm = 0, ff = 0, hh = 0, mm.a = 0, ff.a = 0, hh.a = 0)
        ##                           pi.m.bef = 0, pi.f.bef = 0, pi.m.part = 0, pi.f.part = 0, pi.m.exc = 0, pi.f.exc = 0,
        ##                           pi.m.bef.a = 0, pi.f.bef.a = 0, pi.m.part.a = 0, pi.f.part.a = 0, pi.m.exc.a = 0, pi.f.exc.a = 0) 
        ssL <- 1
        mmL <- 0
        ffL <- 0
        hhL <- 0
        mm.aL <- 0
        ff.aL <- 0
        hh.aL <- 0
        ss <- 1
        mm <- 0
        ff <- 0
        hh <- 0
        mm.a <- 0
        ff.a <- 0
        hh.a <- 0
        for(tt in 1:max(temp$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- temp$tmar-temp$bd+tt-1 >= temp$tms & temp$tmar-temp$bd+tt-1 < temp$tmar
            f.sex <- temp$tmar-temp$bd+tt-1 >= temp$tfs & temp$tmar-temp$bd+tt-1 < temp$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- 0
            p.f.bef <- 0    
            p.m.bef.a <- 0
            p.f.bef.a <- 0    
            p.m.bef[m.sex] <- (1 - exp(-bmb * epicf[cbind(temp$tmar[m.sex]-temp$bd[m.sex]+tt-1, temp$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * epicm[cbind(temp$tmar[f.sex]-temp$bd[f.sex]+tt-1, temp$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(temp$mage[m.sex]-temp$cd[m.sex]-temp$bd[m.sex]+tt-1, temp$cd[m.sex]+temp$bd[m.sex]-tt+1)]
            p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(temp$fage[f.sex]-temp$cd[f.sex]-temp$bd[f.sex]+tt-1, temp$cd[f.sex]+temp$bd[f.sex]-tt+1)]
            ## iterate probabilities based on previous values for only cases where it needs uptemping
            ss[e.sex] <- ssL[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mm[e.sex] <- mmL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            ff[e.sex] <- ffL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            hh[e.sex] <- hhL[e.sex] + ssL[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
              mmL[e.sex]*p.f.bef[e.sex] +
                ffL[e.sex]*p.m.bef[e.sex]
######################################################################
            ## iterate joint probabilities with survival
            mm.a[e.sex] <- mm.aL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
            ff.a[e.sex] <- ff.aL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
            hh.a[e.sex] <- hh.aL[e.sex] + ssL[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
              mm.aL[e.sex]*p.f.bef.a[e.sex] +
                ff.aL[e.sex]*p.m.bef.a[e.sex]
            ## uptempe last month states
            ssL[e.sex] <- ss[e.sex]
            mmL[e.sex] <- mm[e.sex]
            ffL[e.sex] <- ff[e.sex]
            hhL[e.sex] <- hh[e.sex]
            mm.aL[e.sex] <- mm.a[e.sex]
            ff.aL[e.sex] <- ff.a[e.sex]
            hh.aL[e.sex] <- hh.a[e.sex]
            out <- rbind(out, c(cmc = temp$tmar-temp$bd+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        ## probability of infection before partnership
        pi.m.bef <- mm+hh
        pi.f.bef <- ff+hh
        ## probability of infection before partnership & both individuals living to interview
        pi.m.bef.a <- mm.a+hh.a
        pi.f.bef.a <- ff.a+hh.a
        ## probability of infection by partner
        pi.m.part <- 0
        pi.f.part <- 0
        ## probability of infection by partner & both individuals living to interview
        pi.m.part.a <- 0
        pi.f.part.a <- 0
        ## probability of infection extracouply
        pi.m.exc <- 0
        pi.f.exc <- 0
        ## probability of infection extracouply & both individuals living to interview
        pi.m.exc.a <- 0
        pi.f.exc.a <- 0
        ## Track serodiscordant couples by route of infection for bookkeeping later
        mm.bef.a <- mm.a
        ff.bef.a <- ff.a    
        mm.exc.a <- 0
        ff.exc.a <- 0
        mm.bef.aL <- mm.a
        ff.bef.aL <- ff.a    
        mm.exc.aL <- 0
        ff.exc.aL <- 0
        ## probability of being infected by partner (constant, used inside loop)
        p.m.part <- 1 - exp(-bmp)
        p.f.part <- 1 - exp(-bfp)
        ## Now loop through marriage
        for(tt in 1:max(temp$cd-1))
          {
            ## are partners formed in a couple?
            fmd <- temp$cd >= tt              
######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## Survival probabilities
            s.p.m <- csurv[cbind(temp$mage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            s.p.f <- csurv[cbind(temp$fage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            ## Transmission probabilities from partner (jointly with survival)
            p.m.part.a <- p.m.part * s.p.m
            p.f.part.a <- p.f.part * s.p.f
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme*epicf[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe*epicm[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.m.exc.a <- p.m.exc * s.p.m
            p.f.exc.a <- p.f.exc * s.p.f
######################################################################
            ## iterate probabilities
            ss[fmd] <- ssL[fmd]*(1-p.m.exc)*(1-p.f.exc)
            mm[fmd] <- mmL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc*(1-p.f.exc)
            ff[fmd] <- ffL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc*(1-p.m.exc)
            hh[fmd] <- hhL[fmd] + ssL[fmd]* p.m.exc*p.f.exc +
              mmL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc) +
                ffL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)
######################################################################
            ## iterate probabilities jointly with survival until survey
            ## Note for probabilities of not being infected, we don't use the joint probability with being alive at sampling
            mm.a[fmd] <- mm.aL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc.a*(1-p.f.exc)
            ff.a[fmd] <- ff.aL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc.a*(1-p.m.exc)
            hh.a[fmd] <- hh.aL[fmd] + ssL[fmd]  *p.m.exc.a*p.f.exc.a +
              mm.aL[fmd]*(p.f.part.a + (1-p.f.part)*p.f.exc.a) +
                ff.aL[fmd]*(p.m.part.a + (1-p.m.part)*p.m.exc.a)
            ## Track cumulative incidence probabilities for partner infections
            pi.m.part[fmd] <- pi.m.part[fmd] + ffL[fmd]*p.m.part
            pi.f.part[fmd] <- pi.f.part[fmd] + mmL[fmd]*p.f.part        
            pi.m.part.a[fmd] <- pi.m.part.a[fmd] + ff.aL[fmd]*p.m.part.a
            pi.f.part.a[fmd] <- pi.f.part.a[fmd] + mm.aL[fmd]*p.f.part.a
            ## Track how many serodiscordant couples there are with
            ## infections from before partnership 
            mm.bef.a[fmd] <- mm.bef.aL[fmd]*(1-p.f.exc)*(1-p.f.part)
            ff.bef.a[fmd] <- ff.bef.aL[fmd]*(1-p.m.exc)*(1-p.m.part)
            ## Track cumulative incidence probabilities for infections
            ## prior to partnership just subtracting off partners with
            ## incident infections that die and remove couple from
            ## sample-able couples
            pi.m.bef.a[fmd] <- pi.m.bef.a[fmd] - mm.bef.aL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc)*(1-s.p.f)
            pi.f.bef.a[fmd] <- pi.f.bef.a[fmd] - ff.bef.aL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)*(1-s.p.m)
            ## Track cumulative incidence probabilties for extracouple infections
            mm.exc.a[fmd] <- mm.aL[fmd]-mm.bef.aL[fmd]
            ff.exc.a[fmd] <- ff.aL[fmd]-ff.bef.aL[fmd]        
            pi.m.exc[fmd] <- pi.m.exc[fmd] + (ssL[fmd] + ffL[fmd]*(1-p.m.part))*p.m.exc
            pi.f.exc[fmd] <- pi.f.exc[fmd] + (ssL[fmd] + mmL[fmd]*(1-p.f.part))*p.f.exc
            pi.m.exc.a[fmd] <- pi.m.exc.a[fmd] + (ssL[fmd] + ff.aL[fmd]*(1-p.m.part))*p.m.exc.a -
              mm.exc.aL[fmd]*(1-p.f.part)*p.f.exc*(1-s.p.f)
            pi.f.exc.a[fmd] <- pi.f.exc.a[fmd] + (ssL[fmd] + mm.aL[fmd]*(1-p.f.part))*p.f.exc.a -
              ff.exc.aL[fmd]*(1-p.m.part)*p.m.exc*(1-s.p.m)       
            ## update last month states
            ssL[fmd] <- ss[fmd]
            mmL[fmd] <- mm[fmd]
            ffL[fmd] <- ff[fmd]
            hhL[fmd] <- hh[fmd]       
            mm.aL[fmd] <- mm.a[fmd]
            ff.aL[fmd] <- ff.a[fmd]
            hh.aL[fmd] <- hh.a[fmd]
            mm.exc.aL[fmd] <- mm.exc.a[fmd]
            ff.exc.aL[fmd] <- ff.exc.a[fmd]            
            mm.bef.aL[fmd] <- mm.bef.a[fmd]
            ff.bef.aL[fmd] <- ff.bef.a[fmd]
            out <- rbind(out, c(cmc = temp$tmar+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        pser.a <- cbind(hh.a, mm.a, ff.a, ss)
        pser <- cbind(hh, mm, ff, ss)
        if(cprob)
          {
            scalar <- out$mm.a + out$ff.a + out$hh.a + out$ss
          }else{
            scalar <- 1
          }
        if(show.susc)   lines(out$cmc, out$ss/scalar, col = cols[4], lty = 1, lwd = 2)
        if(nsurv)
          {
            lines(out$cmc, out$mm, col = cols[2], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$ff, col = cols[3], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$hh, col = cols[1], lty = lty.nsurv, lwd = 2)
            if(!survive & blob) points(out$cmc[nrow(out)], pser[temp$ser], pch = 19,
                               col = cols[temp$ser], cex = 3)
          }
        if(surv)
          {
            lines(out$cmc, out$mm.a/scalar, col = cols[2], lty = lty.surv, lwd = 2)
            lines(out$cmc, out$ff.a/scalar, col = cols[3], lty = lty.surv, lwd = 2)
            lines(out$cmc, out$hh.a/scalar, col = cols[1], lty = lty.surv, lwd = 2)
            if(survive & blob) points(out$cmc[nrow(out)], pser.a[temp$ser]/scalar[nrow(out)], pch = 19,
                   col = cols[temp$ser], cex = blob.cex)
          }
        if(dead)
          {
             lines(out$cmc, out$mm + out$ff + out$hh - out$mm.a - out$ff.a - out$hh.a, col = col.dead, lty = lty.dead, lwd = 2)
           }
        segments(temp$tmar,-.01,temp$tmar,-drp*.1, lwd = t.lwd)
        text(temp$tmar,-drp*.1,"couple formation",pos=1)
        segments(temp$tms,-.01,temp$tms,-drp*.2, col = cols[2], lwd = t.lwd)
        text(temp$tms,-drp*.2,"male first sex",pos=2, col = cols[2])
        segments(temp$tfs,-.01,temp$tfs,-drp*.3, col = cols[3], lwd = t.lwd)
        text(temp$tfs,-drp*.3,"female first sex",pos=2, col = cols[3])
        segments(temp$tint,-.01,temp$tint,-drp*.2, lwd = t.lwd)
        text(temp$tint,-drp*.2,"DHS survey",pos=4, col = "black")
        if(show.age)
          {
            mtext(paste("male is", round(temp$mage/12,0),"yrs"), adj = .2, side = 3, line = -3)
            mtext(paste("female is", round(temp$fage/12,0),"yrs"), adj = .2, side = 3, line = -5)
          }
      }
    if(do.pdf) dev.off()
  }


## plot trajectory for selected couples
ctraj.area <- function(pars, dat, browse = F, width = 11, height = 8.5, last.year = 2012,
                       plot.cpls = NULL, ylab = "probability in couple serostatus / pop prevalence",
                       surv = T,             # show curves joint with survival
                       nsurv = F,            # show curves not joint with survival
                       survive = T,          # show point on surv curve
                       dead = T,          # show curve for prob of couple dying
                       cprob = T,         # show conditional probabilities
                       cols = c("orange", "darkgreen", "purple","gray"), cex.t = .8,
                       col.dead = "black",
                       lty.surv = 1, lty.dead = 1,
                       lty.nsurv = 2, blob = T, blob.cex = 3, show.susc = T, cex.leg = 1,
                       show.age = F, do.leg = T, drp = 1, x.ticks = T, t.lwd = 1,
                       do.pdf = T, pdf.name = "surv traj.pdf")
  {
    if(browse) browser()
    K <- nrow(dat)
    if(do.pdf)
      {
        pdf(pdf.name, width = width, height = height)
        par(mar=c(0,4,.5,2))
      }
    for(cpl in plot.cpls)
      {
        temp <- dat[cpl,]
        xlim <- c(80*12, (last.year-1900)*12)
        if(show.susc)
          {
            ylim <- c(-.3,1)
          }else{
            ylim <- c(-.3,.5)
          }
        plot(0,0, type = "n", xlim = xlim, axes = F,
             ylim = ylim, xlab = "",  ylab = ylab, bty = "n")
        if(surv & !nsurv & dead & do.leg)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F- & alive", "dead before DHS sampling"),
                                                        col = c(cols,col.dead), lwd = 3,
                                                        lty = c(rep(lty.surv,3),1,lty.dead), bty = "n", cex = cex.leg)
        if(surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+","M+F+ & alive",
                                                             "M+F- & alive","M-F+ & alive","M-F-"),
                                                col = c(cols[1:3],cols), lwd = 3,
                                                lty = c(rep(lty.nsurv,3),rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(surv & !nsurv & do.leg & !dead)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F-"),
                                                         col = cols, lwd = 3,
                                                         lty = c(rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(!surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+", "M-F-"),
                                                 col = cols, lwd = 3,
                                                 lty = c(rep(lty.nsurv,3),1), bty = "n", cex = cex.leg)
        if(do.leg) legend("left", c("HIV population prev in M", "HIV population prev in F"),
                          lty = 2:3, lwd = 3, col = "red", bty = "n", cex = cex.leg)
        if(x.ticks)
          {
            axis(1, seq(xlim[1],xlim[2], by = 12*10), seq(xlim[1],xlim[2], by = 12*10)/12 + 1900, pos = 0, las = 2, cex.axis = cex.leg)
            axis(1, seq(xlim[1],xlim[2], by = 60), labels = NA, pos = 0, las = 2)            
          }else{
            axis(1, seq(xlim[1],xlim[2], by = 60), labels = NA, pos = 0, las = 2)
          }
        axis(2, seq(0,1,l=3), las = 2, cex.axis = cex.leg)
        axis(2, seq(0,1,l=5), labels = NA)
        bmb <- as.numeric(pars["bmb"])
        bfb <- as.numeric(pars["bfb"])
        bme <- as.numeric(pars["bme"])
        bfe <- as.numeric(pars["bfe"])
        bmp <- as.numeric(pars["bmp"])
        rho <- exp(as.numeric(pars["lrho"]))
        bfp <- bmp * rho
        out <- data.frame(cmc = temp$tmar - temp$bd, ss = 1, mm = 0, ff = 0, hh = 0, mm.a = 0, ff.a = 0, hh.a = 0)
        ##                           pi.m.bef = 0, pi.f.bef = 0, pi.m.part = 0, pi.f.part = 0, pi.m.exc = 0, pi.f.exc = 0,
        ##                           pi.m.bef.a = 0, pi.f.bef.a = 0, pi.m.part.a = 0, pi.f.part.a = 0, pi.m.exc.a = 0, pi.f.exc.a = 0) 
        ssL <- 1
        mmL <- 0
        ffL <- 0
        hhL <- 0
        mm.aL <- 0
        ff.aL <- 0
        hh.aL <- 0
        ss <- 1
        mm <- 0
        ff <- 0
        hh <- 0
        mm.a <- 0
        ff.a <- 0
        hh.a <- 0
        for(tt in 1:max(temp$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- temp$tmar-temp$bd+tt-1 >= temp$tms & temp$tmar-temp$bd+tt-1 < temp$tmar
            f.sex <- temp$tmar-temp$bd+tt-1 >= temp$tfs & temp$tmar-temp$bd+tt-1 < temp$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- 0
            p.f.bef <- 0    
            p.m.bef.a <- 0
            p.f.bef.a <- 0    
            p.m.bef[m.sex] <- (1 - exp(-bmb * epicf[cbind(temp$tmar[m.sex]-temp$bd[m.sex]+tt-1, temp$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * epicm[cbind(temp$tmar[f.sex]-temp$bd[f.sex]+tt-1, temp$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(temp$mage[m.sex]-temp$cd[m.sex]-temp$bd[m.sex]+tt-1, temp$cd[m.sex]+temp$bd[m.sex]-tt+1)]
            p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(temp$fage[f.sex]-temp$cd[f.sex]-temp$bd[f.sex]+tt-1, temp$cd[f.sex]+temp$bd[f.sex]-tt+1)]
            ## iterate probabilities based on previous values for only cases where it needs uptemping
            ss[e.sex] <- ssL[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mm[e.sex] <- mmL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            ff[e.sex] <- ffL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            hh[e.sex] <- hhL[e.sex] + ssL[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
              mmL[e.sex]*p.f.bef[e.sex] +
                ffL[e.sex]*p.m.bef[e.sex]
######################################################################
            ## iterate joint probabilities with survival
            mm.a[e.sex] <- mm.aL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
            ff.a[e.sex] <- ff.aL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
            hh.a[e.sex] <- hh.aL[e.sex] + ssL[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
              mm.aL[e.sex]*p.f.bef.a[e.sex] +
                ff.aL[e.sex]*p.m.bef.a[e.sex]
            ## uptempe last month states
            ssL[e.sex] <- ss[e.sex]
            mmL[e.sex] <- mm[e.sex]
            ffL[e.sex] <- ff[e.sex]
            hhL[e.sex] <- hh[e.sex]
            mm.aL[e.sex] <- mm.a[e.sex]
            ff.aL[e.sex] <- ff.a[e.sex]
            hh.aL[e.sex] <- hh.a[e.sex]
            out <- rbind(out, c(cmc = temp$tmar-temp$bd+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        ## probability of infection before partnership
        pi.m.bef <- mm+hh
        pi.f.bef <- ff+hh
        ## probability of infection before partnership & both individuals living to interview
        pi.m.bef.a <- mm.a+hh.a
        pi.f.bef.a <- ff.a+hh.a
        ## probability of infection by partner
        pi.m.part <- 0
        pi.f.part <- 0
        ## probability of infection by partner & both individuals living to interview
        pi.m.part.a <- 0
        pi.f.part.a <- 0
        ## probability of infection extracouply
        pi.m.exc <- 0
        pi.f.exc <- 0
        ## probability of infection extracouply & both individuals living to interview
        pi.m.exc.a <- 0
        pi.f.exc.a <- 0
        ## Track serodiscordant couples by route of infection for bookkeeping later
        mm.bef.a <- mm.a
        ff.bef.a <- ff.a    
        mm.exc.a <- 0
        ff.exc.a <- 0
        mm.bef.aL <- mm.a
        ff.bef.aL <- ff.a    
        mm.exc.aL <- 0
        ff.exc.aL <- 0
        ## probability of being infected by partner (constant, used inside loop)
        p.m.part <- 1 - exp(-bmp)
        p.f.part <- 1 - exp(-bfp)
        ## Now loop through marriage
        for(tt in 1:max(temp$cd-1))
          {
            ## are partners formed in a couple?
            fmd <- temp$cd >= tt              
######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## Survival probabilities
            s.p.m <- csurv[cbind(temp$mage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            s.p.f <- csurv[cbind(temp$fage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            ## Transmission probabilities from partner (jointly with survival)
            p.m.part.a <- p.m.part * s.p.m
            p.f.part.a <- p.f.part * s.p.f
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme*epicf[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe*epicm[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.m.exc.a <- p.m.exc * s.p.m
            p.f.exc.a <- p.f.exc * s.p.f
######################################################################
            ## iterate probabilities
            ss[fmd] <- ssL[fmd]*(1-p.m.exc)*(1-p.f.exc)
            mm[fmd] <- mmL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc*(1-p.f.exc)
            ff[fmd] <- ffL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc*(1-p.m.exc)
            hh[fmd] <- hhL[fmd] + ssL[fmd]* p.m.exc*p.f.exc +
              mmL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc) +
                ffL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)
######################################################################
            ## iterate probabilities jointly with survival until survey
            ## Note for probabilities of not being infected, we don't use the joint probability with being alive at sampling
            mm.a[fmd] <- mm.aL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc.a*(1-p.f.exc)
            ff.a[fmd] <- ff.aL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc.a*(1-p.m.exc)
            hh.a[fmd] <- hh.aL[fmd] + ssL[fmd]  *p.m.exc.a*p.f.exc.a +
              mm.aL[fmd]*(p.f.part.a + (1-p.f.part)*p.f.exc.a) +
                ff.aL[fmd]*(p.m.part.a + (1-p.m.part)*p.m.exc.a)
            ## Track cumulative incidence probabilities for partner infections
            pi.m.part[fmd] <- pi.m.part[fmd] + ffL[fmd]*p.m.part
            pi.f.part[fmd] <- pi.f.part[fmd] + mmL[fmd]*p.f.part        
            pi.m.part.a[fmd] <- pi.m.part.a[fmd] + ff.aL[fmd]*p.m.part.a
            pi.f.part.a[fmd] <- pi.f.part.a[fmd] + mm.aL[fmd]*p.f.part.a
            ## Track how many serodiscordant couples there are with
            ## infections from before partnership 
            mm.bef.a[fmd] <- mm.bef.aL[fmd]*(1-p.f.exc)*(1-p.f.part)
            ff.bef.a[fmd] <- ff.bef.aL[fmd]*(1-p.m.exc)*(1-p.m.part)
            ## Track cumulative incidence probabilities for infections
            ## prior to partnership just subtracting off partners with
            ## incident infections that die and remove couple from
            ## sample-able couples
            pi.m.bef.a[fmd] <- pi.m.bef.a[fmd] - mm.bef.aL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc)*(1-s.p.f)
            pi.f.bef.a[fmd] <- pi.f.bef.a[fmd] - ff.bef.aL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)*(1-s.p.m)
            ## Track cumulative incidence probabilties for extracouple infections
            mm.exc.a[fmd] <- mm.aL[fmd]-mm.bef.aL[fmd]
            ff.exc.a[fmd] <- ff.aL[fmd]-ff.bef.aL[fmd]        
            pi.m.exc[fmd] <- pi.m.exc[fmd] + (ssL[fmd] + ffL[fmd]*(1-p.m.part))*p.m.exc
            pi.f.exc[fmd] <- pi.f.exc[fmd] + (ssL[fmd] + mmL[fmd]*(1-p.f.part))*p.f.exc
            pi.m.exc.a[fmd] <- pi.m.exc.a[fmd] + (ssL[fmd] + ff.aL[fmd]*(1-p.m.part))*p.m.exc.a -
              mm.exc.aL[fmd]*(1-p.f.part)*p.f.exc*(1-s.p.f)
            pi.f.exc.a[fmd] <- pi.f.exc.a[fmd] + (ssL[fmd] + mm.aL[fmd]*(1-p.f.part))*p.f.exc.a -
              ff.exc.aL[fmd]*(1-p.m.part)*p.m.exc*(1-s.p.m)       
            ## update last month states
            ssL[fmd] <- ss[fmd]
            mmL[fmd] <- mm[fmd]
            ffL[fmd] <- ff[fmd]
            hhL[fmd] <- hh[fmd]       
            mm.aL[fmd] <- mm.a[fmd]
            ff.aL[fmd] <- ff.a[fmd]
            hh.aL[fmd] <- hh.a[fmd]
            mm.exc.aL[fmd] <- mm.exc.a[fmd]
            ff.exc.aL[fmd] <- ff.exc.a[fmd]            
            mm.bef.aL[fmd] <- mm.bef.a[fmd]
            ff.bef.aL[fmd] <- ff.bef.a[fmd]
            out <- rbind(out, c(cmc = temp$tmar+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        pser.a <- cbind(hh.a, mm.a, ff.a, ss)
        pser <- cbind(hh, mm, ff, ss)
        if(cprob)
          {
            scalar <- out$mm.a + out$ff.a + out$hh.a + out$ss
          }else{
            scalar <- 1
          }
                                        #        if(show.susc)   lines(out$cmc, out$ss/scalar, col = cols[4], lty = 1, lwd = 2)
        if(nsurv)
          {
            lines(out$cmc, out$mm, col = cols[2], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$ff, col = cols[3], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$hh, col = cols[1], lty = lty.nsurv, lwd = 2)
            if(!survive & blob) points(out$cmc[nrow(out)], pser[temp$ser], pch = 19,
                                       col = cols[temp$ser], cex = 3)
          }
        if(surv)
          {
            if(out$cmc[1] > xlim[1])
              {
                out.begin <- out[rep(1, out$cmc[1]-xlim[1]),]
                out.begin$cmc <- xlim[1]:(out$cmc[1]-1)
                outpl <- rbind(out.begin,out)
              }else{
                outpl <- out[out$cmc >= xlim[1],]
              }
            scalar <- outpl$mm.a + outpl$ff.a + outpl$hh.a + outpl$ss
            polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar, rep(0, nrow(outpl))), col = cols[1])
            polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar, rev(outpl$hh.a/scalar + outpl$ff.a/scalar)), col = cols[3])
            polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar + outpl$ff.a/scalar, rev(outpl$hh.a/scalar + outpl$ff.a/scalar + outpl$mm.a/scalar)),  col = cols[2])
            if(show.susc)   polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar + outpl$ff.a/scalar + outpl$mm.a/scalar,
                                                                rev(outpl$hh.a/scalar + outpl$ff.a/scalar + outpl$mm.a/scalar + outpl$ss/scalar)), col = cols[4])
            if(survive & blob) points(outpl$cmc[nrow(outpl)]+12*2, pser.a[temp$ser]/scalar[nrow(outpl)], pch = 19, lwd = 1,
                   col = cols[temp$ser], cex = blob.cex)
            if(survive & blob) points(outpl$cmc[nrow(outpl)]+12*2, pser.a[temp$ser]/scalar[nrow(outpl)], pch = 21, lwd = 1,
                   col = "black", cex = blob.cex)
          }
        if(dead)
          {
             lines(outpl$cmc, outpl$mm + outpl$ff + outpl$hh - outpl$mm.a - outpl$ff.a - outpl$hh.a, col = col.dead, lty = lty.dead, lwd = 2)
           }
        segments(temp$tmar,-.01,temp$tmar,-drp, lwd = t.lwd)
        text(temp$tmar,-drp,expression(t),pos=1, cex = cex.t)
        segments(temp$tms,-.01,temp$tms,-drp, col = cols[2], lwd = t.lwd)
        text(temp$tms,-drp,expression(t),pos=1, col = cols[2], cex = cex.t)
        segments(temp$tfs,-.01,temp$tfs,-drp, col = cols[3], lwd = t.lwd)
        text(temp$tfs,-drp,expression(t),pos=1, col = cols[3], cex = cex.t)
        segments(temp$tint,-.01,temp$tint,-drp, lwd = t.lwd)
        text(temp$tint,-drp,expression(t),pos=1, col = "black", cex = cex.t)
#        lines(xlim[1]:xlim[2], epicm[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 2, lwd = 3)
 #       lines(xlim[1]:xlim[2], epicf[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 3, lwd = 3)
        if(show.age)
          {
            mtext(paste("male is", round(temp$mage/12,0),"yrs"), adj = .2, side = 3, line = -3)
            mtext(paste("female is", round(temp$fage/12,0),"yrs"), adj = .2, side = 3, line = -5)
          }
      }
    if(do.pdf) dev.off()
  }

## Plot posterior pairwise-correlations & histograms
sbpairs <- function(posts, truepars = NULL, file.nm, width = 10, height = 10, show.lines = T,
                    cex = 1, col = "black", nrpoints = 200, do.pdf = F, do.jpeg = T, big.ranges = T,
                    cex.axis = 1.5, cex.nm = 2.5, greek = T, show.cis = T, browse = F)
  {
    if(browse) browser()
    if(do.pdf) pdf(paste(file.nm, ".pdf", sep=""), width = width, height = height)
    if(do.jpeg) jpeg(paste(file.nm, ".jpeg", sep=""), width = width*100, height = height*100)
    ranges <- apply(rbind(posts,truepars), 2, range)
    if(big.ranges)
      {
        ranges[1,1:5] <- 0
        ranges[2,1:5] <- .06 #1.5* ranges[2,1:5]
        ranges[1,6] <- ranges[1,6] - .5
        ranges[2,6] <- ranges[2,6] + .5
      }
    cis <- apply(rbind(posts,truepars), 2, function(x) quantile(x, c(.025, .975)))
    par(mar=rep(3,4),oma=rep(2,4),mfrow=rep(ncol(posts),2))
    parnames <- c(expression(beta[Mb]), expression(beta[Fb]),
                  expression(beta[Me]), expression(beta[Fe]),
                  expression(beta[Mp]), expression(log(rho)))
    for(ii in 1:ncol(posts))
      {
        for(jj in 1:ncol(posts))
          {
            if(ii==jj)
              {
                hist(posts[,ii], breaks = 40, col = "black", main = "", xlab = "", ylab = "",
                     xlim = ranges[,ii], las = 2, cex.axis = cex.axis)
                if(show.lines) abline(v=truepars[ii], col = "red", lwd = 3)
                if(show.cis) abline(v=cis[,ii], col = "yellow", lwd = 3)
                mtext(parnames[ii], side = 3, line = -2, adj = .98, col = "red", cex = cex.nm)
              }else
            {
              smoothScatter(posts[,jj],posts[,ii], cex = cex, col = col, las = 2, cex.axis = cex.axis,
                            main = "", xlab = "", ylab = "", nrpoints = nrpoints,
                            xlim = ranges[,jj], ylim = ranges[,ii])
              if(show.lines) abline(v=truepars[jj], col = "red", lwd = 2)
              if(show.lines) abline(h=truepars[ii], col = "red", lwd = 2)              
            }
          }
      }
    if(do.pdf | do.jpeg) dev.off()
  }


## Make scatterplot with marginal histograms
    scatterhist = function(x, y, xlab="", ylab="", xb = 5, yb = .1, border = NA, cex = .4, xlim = c(0,60), ylim = c(0,1), ..., nn=F,
      col = "black", browse = F)
      {
        if(browse) browser()
        par(oma=c(1,1,0,0))        
        zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
        layout(zones, widths=c(3.5,2), heights=c(2,3.5))
        xbrks <- seq(min(xlim),max(xlim), by = xb)
        ybrks <- seq(min(ylim),max(ylim), by = yb)        
        xhist = hist(x, plot=FALSE, breaks = xbrks)
        yhist = hist(y, plot=FALSE, breaks = ybrks)
        top = max(c(xhist$counts, yhist$counts))
        if(nn) top <- max(xhist$counts)
        par(mar=c(2,2,1,1))
        if(!nn)
          {            
            plot(x,y, bty = "n", ylim = ylim, xlim = xlim, pch = 19, cex = cex, las = 2, axes=T, col = col)
          }else{
            plot(x,y, bty = "n", axes = F, type = "n")
          }
        par(mar=c(0,2,1,1))
        bp <- barplot(xhist$counts, axes=F, ylim=c(0, top), space=0, border = border, col = "black")
        show <- 1:length(xbrks) %% 2 == 1
        axis(1, (c(bp,max(bp)+1)-.5)[show], NA, las = 2)
        if(!nn)
          {
            par(mar=c(2,0,1,1))
            bp <- barplot(yhist$counts, axes=F, xlim=c(0, top), space=0, horiz=TRUE, border = border, col = "black")
            show <- 1:length(ybrks) %% 2 == 1
            axis(2, (c(bp,max(bp)+1)-.5)[show], NA, las = 2)
          }
      }
