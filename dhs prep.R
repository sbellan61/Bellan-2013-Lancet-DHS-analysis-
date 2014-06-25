## Steve Bellan, March 2012
## Prepare data sets for JAGS analysis. Going to have a single data
## frame for all countries having a DHS set and list with each object
## epidemic curve. Need to check for errors.
rm(list=ls())
dirnm <- "~/Dropbox/Inactive Projects/Finished Projects/Lancet MS/Couples Model Revision 121013/R scripts & Data Files/"
pth <- "~/Dropbox/Inactive Projects/Finished Projects/Lancet MS/Couples Model Revision 121013/R scripts & Data Files/DHS sets/"
fnms <- list.files(path = pth)
fnms
cnms <- fnms
cnms <- sub("_", " ", cnms);cnms <- sub("_", " ", cnms);cnms <- sub("_", " ", cnms);
cnms <- sub("_", " ", cnms)
cnms <- sub(".select.RData", "", cnms)
cnms

## assume that coverage is half of what's reported (due to adherence, etc...)
coverage.eff <- 1

setwd(dirnm)
## put all dhs data into the same data frame
for(ii in 1:length(fnms))
  {
    load(paste(pth,fnms[ii], sep = ""))
    ## for some reason these columns aren't in all of them, don't need em anyways
    rem.col <-names(Answers) %in% c("M_lives_with_spouse","How_previous_marriage_or_union_ended",
                                    "Dataset")
    temp <- data.frame(ds = cnms[ii], Answers[,!rem.col])
    if(ii==1) raw <- temp else raw <- rbind(raw, temp)
  }

names(raw) <- sub("_", "", names(raw))
## replace remaining with periods
while(sum(grepl("_", names(raw))>0))     names(raw) <- sub("_",".", names(raw))
names(raw)
xtabs(~raw$MmaritalHistory)
xtabs(~raw$FmaritalHistory)
xtabs(~raw$ds)
raw <- data.frame(uid = 1:nrow(raw),raw)
oraw <- raw

######################################################################
## get prevalence variables from UNAIDS epidemic curves
######################################################################
spectr <- read.csv("~/Dropbox/Inactive Projects/Finished Projects/Lancet MS/Couples Model Revision 121013/R scripts & Data Files/SSA prevalence.csv")
dim(spectr)
## remove empty rows
spectr <- na.exclude(spectr)
## remove south sudan (no data)
spectr <- spectr[!spectr$country=="Sudan South",]
dim(spectr)
names(spectr)
names(spectr) <- c("year","country","all.hiv.pop","m.hiv.pop","f.hiv.pop",
                   "all.pop","m.pop","f.pop",
                   "all.prev","mprev.all","fprev.all", "art.pop")
for(jj in 3:ncol(spectr)) spectr[,jj] <- as.numeric(levels(spectr[,jj])[spectr[,jj]])
## ART normalized prevalence, assume % of female on ART = % of ppl with HIV that are female * # on ARVs
spectr$perchiv.f <- spectr$f.hiv.pop / (spectr$m.hiv.pop + spectr$f.hiv.pop)
spectr$perchiv.m <- 1-spectr$perchiv.f
spectr$mprev <- (spectr$m.hiv.pop - spectr$perchiv.m * coverage.eff * spectr$art.pop) / spectr$m.pop
spectr$fprev <- (spectr$f.hiv.pop - spectr$perchiv.f * coverage.eff * spectr$art.pop) / spectr$f.pop
## art prevalences...
spectr$art.prev <- spectr$art.pop/spectr$all.hiv.pop
write.csv(spectr[spectr$year==2009, c("country","art.prev")], file = "2009 art prev.csv")
## only show columns that we're using
spectr[,c("mprev.all","fprev.all")] <- spectr[,c("mprev.all","fprev.all")]/100 # make proportion
spectr <- spectr[,c("year","country","mprev","fprev","mprev.all","fprev.all","art.prev")]

logit <- function(x) log(x/(1-x))
ilogit <- function(x) exp(x) / (1 + exp(x))
x.seq <- 1:((2013-1900)*12)
epicm <- data.frame(cmc = x.seq)
epicf <- data.frame(cmc = x.seq)
epicm.all <- data.frame(cmc = x.seq)
epicf.all <- data.frame(cmc = x.seq)
art.prev <- data.frame(cmc = x.seq)         # proportion not on ART (1-coverage)
pdf("ARV coverage model.pdf")
for(ii in 1:length(unique(spectr$country)))
  {
    cc <- unique(spectr$country)[ii]
    temp <- spectr[spectr$country==cc,]
    plot(temp$year, temp$art.prev, pch = 19, main = cc, xlim = c(1995,2011), ylim = c(0, .6))
    ## add 0 prevalence in 1980 to stabilize
    temp <- rbind(data.frame(year = 1980, country = cc, mprev = 10^-6, fprev = 10^-6, mprev.all = 10^-6, fprev.all = 10^-6,
                             art.prev = 10^-8), temp)
    ## make art prevalence 10^-8 to stabilize
    temp$art.prev[temp$art.prev==0] <- 10^-8
    temp$cmc <- (temp$year -1900)*12
    temp$logitm <- logit(temp$mprev)
    temp$logitf <- logit(temp$fprev)
    temp$logitm.all <- logit(temp$mprev.all)
    temp$logitf.all <- logit(temp$fprev.all)
    temp$art.prev <- logit(temp$art.prev)
    mmod <- smooth.spline(temp$logitm ~ temp$cmc)
    fmod <- smooth.spline(temp$logitf ~ temp$cmc)
    mmod.all <- smooth.spline(temp$logitm.all ~ temp$cmc)
    fmod.all <- smooth.spline(temp$logitf.all ~ temp$cmc)
    art.prev.mod <- smooth.spline(temp$art.prev ~ temp$cmc, df = 18)
    m.pred <- predict(mmod, x.seq)
    f.pred <- predict(fmod, x.seq)
    m.pred$y <- ilogit(m.pred$y)
    f.pred$y <- ilogit(f.pred$y)
    m.pred.all <- predict(mmod.all, x.seq)
    f.pred.all <- predict(fmod.all, x.seq)
    m.pred.all$y <- ilogit(m.pred.all$y)
    f.pred.all$y <- ilogit(f.pred.all$y)
    art.prev.pred <- predict(art.prev.mod, x.seq)
    art.prev.pred$y <- ilogit(art.prev.pred$y)
    art.prev.pred$y[art.prev.pred$x>109*12] <- art.prev.pred$y[109*12] * 1.01^(art.prev.pred$x[art.prev.pred$x>109*12] - 109*12)
    ## now use output of loess predicted at months to get normalized
    ## prevalence-months for each duration
    epicm <- data.frame(epicm, m.pred$y)
    epicf <- data.frame(epicf, f.pred$y)    
    epicm.all <- data.frame(epicm.all, m.pred.all$y)
    epicf.all <- data.frame(epicf.all, f.pred.all$y)
    art.prev <- data.frame(art.prev, art.prev.pred$y)
    lines(1900+ art.prev[,1]/12, art.prev[,ncol(art.prev)], col = "red")
  }
dev.off()
names(epicm)[-1] <- as.character(unique(spectr$country))
names(epicf)[-1] <- as.character(unique(spectr$country))
names(art.prev)[-1] <- as.character(unique(spectr$country))

## must match DHS countries to spectr countries
names(epicm)[names(epicm)=="Democratic Republic of the Congo"] <- "DRC"
names(epicf)[names(epicf)=="Democratic Republic of the Congo"] <- "DRC"
names(art.prev)[names(art.prev)=="Democratic Republic of the Congo"] <- "DRC"
names(epicm)[grepl('Tanzania', names(epicm))] <- "Tanzania"
names(epicf)[grepl('Tanzania', names(epicf))] <- "Tanzania"
names(art.prev)[grepl('Tanzania', names(art.prev))] <- "Tanzania"

unique(raw$ds)

names(epicm.all) <- names(epicm)
names(epicf.all) <- names(epicf)
raw$epic.ind <- NA
for(ii in 2:ncol(epicm))
  {
    cc <- names(epicm)[ii]
    ind <- grepl(cc, raw$ds)
    if(sum(ind)>0)
      {
        raw$epic.ind[ind] <- ii
        raw$epic.nm[ind] <- cc
      }
  }

## make sure the curves match the dhs data
unique(raw[,c("ds","epic.nm")])
###################################################################### 
## pool w africa
######################################################################
w.africa <- c("Burkina Faso 2003", "Cameroon 2004", "Ghana 2003", "Guinea 2005",
              "Liberia 2007", "Mali 2006", "Niger 2006", "Senegal 2005",
              "Sierra Leone 2008")
wa.ind <- unique(raw$epic.ind[raw$ds %in% w.africa])
wa.nm <- unique(raw$epic.nm[raw$ds %in% w.africa])
raw$group <- NA
raw$group[raw$ds %in% w.africa] <- "WA"
other.nm <- unique(raw$epic.nm)[!unique(raw$epic.nm) %in% wa.nm]
for(nn in other.nm)
  {
    raw$group[grepl(nn, raw$ds)] <- nn
  }
unique(raw[,c("ds","epic.nm","group")])
raw$wa <- F
raw$wa[raw$ds %in% w.africa] <- T
oraw <- raw                             # save original data set

######################################################################
## Make Zambia prev curves for Figure 1: Male beginning sex at 1980, married 1990
######################################################################
pdf("mprev zam fig 1.pdf", w = 3, h=1.7)
ymax <- .2
cex.nm <- .6
par(mar = c(1,3,0,0))
plot(0,0, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, axes = F, cex.axis = cex.nm)
axis(1, seq(1975,2010, by = 5), labels = NA)
axis(2, seq(0, .2, by = .05), las = 2, cex.axis = cex.nm)
ind <- 41
head(epicm[,c(1,ind)])
show <- epicm$cmc/12 + 1900 > 1988 & epicm$cmc/12 + 1900 <= 1990
col <- "dark gray"
polygon(c(epicm$cmc[show]/12+1900, rev(epicm$cmc[show]/12+1900)),
          c(epicm[show,ind], rep(0, sum(show))),
          col = col, border = NA)
show <- epicm$cmc/12 + 1900 >= 1990 & epicm$cmc/12 + 1900 <= 2007
col <- "red"
polygon(c(epicm$cmc[show]/12+1900, rev(epicm$cmc[show]/12+1900)),
          c(epicm[show,ind], rep(0, sum(show))),
          col = col, border = NA)
show <- epicm$cmc/12 + 1900 > 1975 & epicm$cmc/12 + 1900 < 2011
lines(epicm.all$cmc[show]/12+1900, epicm[show,ind], col = "black", lty = 2)
lines(epicm.all$cmc[show]/12+1900, epicm.all[show,ind], col = "black", lty = 1)
dev.off()
######################################################################


######################################################################
## Make Zambia prev curves for Figure 1: Female beginning sex at 1988, married 1990
######################################################################
pdf("fprev zam fig 1.pdf", w = 3, h=1.7)
ymax <- .2
cex.nm <- .6
par(mar = c(1,3,0,0))
plot(0,0, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, axes = F, cex.axis = cex.nm)
#axis(1, seq(1975,2010, by = 5))
axis(2, seq(0, .2, by = .05), las = 2, cex.axis = cex.nm)
ind <- 41
head(epicf[,c(1,ind)])
show <- epicf$cmc/12 + 1900 > 1980 & epicf$cmc/12 + 1900 <= 1990
col <- "dark gray"
polygon(c(epicf$cmc[show]/12+1900, rev(epicf$cmc[show]/12+1900)),
          c(epicf[show,ind], rep(0, sum(show))),
          col = col, border = NA)
show <- epicf$cmc/12 + 1900 >= 1990 & epicf$cmc/12 + 1900 <= 2007
col <- "red"
polygon(c(epicf$cmc[show]/12+1900, rev(epicf$cmc[show]/12+1900)),
          c(epicf[show,ind], rep(0, sum(show))),
          col = col, border = NA)
show <- epicf$cmc/12 + 1900 > 1975 & epicf$cmc/12 + 1900 < 2011
lines(epicf.all$cmc[show]/12+1900, epicf[show,ind], col = "black", lty = 2)
lines(epicf.all$cmc[show]/12+1900, epicf.all[show,ind], col = "black", lty = 1)
dev.off()
######################################################################


## make figure of epidemic curves
pdf("ssa epic.pdf", w = 4.5, h=4)
ymax <- .3
cex.nm <- .7
par(mfrow = c(2,1), mar = c(.5,2,1.5,0), oma = c(2,2,0,0))
plot(temp$year, temp$mprev, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
axis(1, seq(1975,2010, by = 5), NA)
countries <- unique(raw$epic.ind[!raw$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicm$cmc/12 + 1900 > 1975 & epicm$cmc/12 + 1900 < 2011
    lines(epicm$cmc[show]/12+1900, epicm[show,ind], col = cols[ii], lty = 2)
    lines(epicm.all$cmc[show]/12+1900, epicm.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("males", side = 3, line = 0, cex = 1)
legend("topleft", names(epicm)[countries], col=cols, lwd=1, bty = "n", ncol = 2, cex = cex.nm)
plot(temp$year, temp$mprev, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "year", ylab = "prevalence", main = "", las = 2, cex.axis = cex.nm)
countries <- unique(raw$epic.ind[!raw$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicf$cmc/12 + 1900 > 1975 & epicf$cmc/12 + 1900 < 2011
    lines(epicf$cmc[show]/12+1900, epicf[show,ind], col = cols[ii], lty = 2)
    lines(epicf.all$cmc[show]/12+1900, epicf.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("females", side = 3, line = 0, cex = 1)
legend("topleft", c("HIV prevalence","HIV prevalence * (1-ARV coverage)"), lty = 1:2, bty = "n", cex = cex.nm)
#mtext("population prevalence", side = 2, outer = T, line = 1, cex = cex.nm)
dev.off()

## make figure of epidemic curves
pdf("wa epic.pdf", w = 4.5, h=4)
ymax <- .08
cex.nm <- .7
par(mfrow = c(2,1), mar = c(.5,2,1.5,0), oma = c(2,2,0,0))
plot(temp$year, temp$mprev, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
axis(1, seq(1975,2010, by = 5), NA)
countries <- unique(raw$epic.ind[raw$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicm$cmc/12 + 1900 > 1975 & epicm$cmc/12 + 1900 < 2011
    lines(epicm$cmc[show]/12+1900, epicm[show,ind], col = cols[ii], lty = 2)
    lines(epicm.all$cmc[show]/12+1900, epicm.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("males", side = 3, line = 0, cex = 1)
legend("topleft", names(epicm)[countries], col=cols, lwd=1, bty = "n", ncol = 2, cex = cex.nm)
plot(temp$year, temp$mprev, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "year", ylab = "prevalence", main = "", las = 2, cex.axis = cex.nm)
countries <- unique(raw$epic.ind[raw$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicf$cmc/12 + 1900 > 1975 & epicf$cmc/12 + 1900 < 2011
    lines(epicf$cmc[show]/12+1900, epicf[show,ind], col = cols[ii], lty = 2)
    lines(epicf.all$cmc[show]/12+1900, epicf.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("females", side = 3, line = 0, cex = 1)
legend("topleft", c("HIV prevalence","HIV prevalence * (1-ARV coverage)"), lty = 1:2, bty = "n", cex = cex.nm)
mtext("population prevalence", side = 2, outer = T, line = 1, cex = cex.nm)
dev.off()

ltab <- data.frame(n = xtabs(~group, raw), nohivtest = NA, polygamy = NA,
                   no1scp = NA, errors = NA, analyzed = NA)
## HIV tested
xtabs(~FHIVresult + group, raw)
xtabs(~MHIVresult + group, raw)
levels(raw$FHIVresult)
nas <-  c("ERROR : V-, W+, M+", "ERROR : V-, W+, M-", "ERROR : V-, W-, M+",
          "Indeterminant", "Not enough samples to complete protocol")
noser <- is.na(raw$FHIVresult) | is.na(raw$MHIVresult)
noser <- noser | raw$FHIVresult %in% nas | raw$MHIVresult %in% nas
ltab$nohivtest <- xtabs(~group + noser, raw)[,2]
## create serostatuses
fp <- grepl("positive", raw$FHIVresult)
mp <- grepl("positive", raw$MHIVresult)
ser <- rep(NA, nrow(raw))
ser[fp & mp] <- 1 #"M+ F+"
ser[!fp & mp] <- 2 #"M+ F-"
ser[fp & !mp] <- 3 #"M- F+"
ser[!fp & !mp] <- 4 #"M- F-"
raw$ser <- ser
aggregate(noser, list(raw$ds), function(x) signif(mean(x),3)*100)
#write.csv(xtabs(~ds + ser, raw), file="~/Documents/R files/discordant couples/reldurmod/ser all breakdown.csv")

## One wife
xtabs(~ds + MnumWives, raw)
aggregate(wives, list(raw$ds), function(x) signif(mean(x),3)*100)
xtabs(~ds + is.na(MnumWives), raw)
wives <- raw$MnumWives > 1 | is.na(raw$MnumWives)
## losstab <- data.frame(losstab, onewife =xtabs(~ds + wives, raw)[,1])
## raw <- raw[!wives,]
ltab$polygamy <- xtabs(~group + wives, raw)[,2]
ltab
both <- !noser & !wives

## Married or Living together (everyone is by definition)
mmar <- raw$MmaritalStatus %in% c("Married", "Living together", "Living with partner")
fmar <- raw$FmaritalStatus %in% c("Married", "Living together", "Living with partner")
mar <- mmar & fmar
xtabs(~ds+mar,raw)
## Create variable indicating both are in their first & only marriage
## first need to make NA's unknowns to allow logical manipulation
levels(raw$Fnumber.of.unions) <- c(levels(raw$Fnumber.of.unions), "UNK")
levels(raw$Mnumber.of.unions) <- c(levels(raw$Mnumber.of.unions), "UNK")
raw$Mnumber.of.unions[is.na(raw$Mnumber.of.unions)] <- "UNK"
raw$Fnumber.of.unions[is.na(raw$Fnumber.of.unions)] <- "UNK"
## First marriage (one of them must be in 1st marriage to get mardur)
first.mar <- raw$Fnumber.of.unions=="Once" | raw$Mnumber.of.unions=="Once"
first.mar <- first.mar & mar
ltab$no1scp <- xtabs(~group + first.mar, raw)[,1]
ltab
aggregate(!first.mar, list(raw$ds), function(x) signif(mean(x),3)*100)


## interview of m & f was <= 1 month apart
tint.diff <- (raw$Minterview.cmc - raw$Finterview.cmc)
int.gr1mon <- abs(tint.diff) > 1 | is.na(tint.diff)
## missing age at first intercourse
no.aafi <- is.na(raw$FImp.Age.at.first.intercourse) | is.na(raw$MImp.Age.at.first.intercourse)
aggregate(no.aafi, list(raw$ds), function(x) signif(mean(x),3)*100)
aggregate(int.gr1mon, list(raw$ds), function(x) signif(mean(x),3)*100)

pdf("dhs problems.pdf")
                                        # they were interviewed about the same time, check for each country
raw$tint <- raw$Minterview.cmc
## calculaate ages using cmc
mage <- (raw$Minterview.cmc-raw$Mdob.cmc)/12
fage <- (raw$Finterview.cmc-raw$Fdob.cmc)/12
###################################################################### 
## calculate marriage durations
mmardur <- (raw$Minterview.cmc-raw$MDate.of.first.marriage)/12
fmardur <- (raw$Finterview.cmc-raw$FDate.of.first.marriage)/12
## should we be worried that for couples where both individuals are in
## their 1st and only marriage, the duration differs by couples a lot?
## removing Ethiopia because forced all males to be in the Once category above
mardiff <- mmardur-fmardur
sum(abs(mardiff[raw$both.first])>2.5)
hist(mmardur-fmardur, breaks = -50:50, col = "black", main = "male marriage duration - female marriage duration")
sum(both)

###################################################################### 
## marriage duration (female listed if she is the only 1 in the first
## marriage & vv or average if both in first marriage)
ffirst <- raw$Fnumber.of.unions=="Once" & !raw$Mnumber.of.unions=="Once"
mfirst <- raw$Mnumber.of.unions=="Once" & !raw$Fnumber.of.unions=="Once"
bfirst <- raw$Mnumber.of.unions=="Once" & raw$Fnumber.of.unions=="Once"
raw$mardur[ffirst] <- fmardur[ffirst]
raw$mardur[mfirst] <- mmardur[mfirst]
raw$mardur[bfirst] <- .5*(mmardur+fmardur)[bfirst]
raw$mardur.mon <- round(raw$mardur*12)  # must round months since averaging m & f month durations
raw$tmar <- raw$tint - raw$mardur.mon    # imputed cmc of marriage using mardur above
aggregate(bfirst, list(raw$ds), function(x) signif(mean(x),3)*100)
sum(both)
###################################################################### 
## now get yrs of sexual activity for each
###################################################################### 
## the 97 and 98 codes seem to indicate that first sexual intercourse
## was when started living with the partner, meaning that yrs of
## active sex = mardur. Note, the code in the survey says 95 not 97/98
## but doesn't seem to be any other option.
## female years of active sex
raw$fysa <- fage - raw$FImp.Age.at.first.intercourse
nna <- !is.na(raw$FImp.Age.at.first.intercourse)
raw$fysa[nna & raw$FImp.Age.at.first.intercourse>90] <- fmardur[nna & raw$FImp.Age.at.first.intercourse>90]
raw$mysa <- mage - raw$MImp.Age.at.first.intercourse
nna <- !is.na(raw$MImp.Age.at.first.intercourse)
raw$mysa[nna & raw$MImp.Age.at.first.intercourse>90] <- mmardur[nna & raw$MImp.Age.at.first.intercourse>90]


###################################################################### 
## time of first sex for m & f in CMC
raw$tfs <- raw$tint - raw$fysa*12
raw$tms <- raw$tint - raw$mysa*12


## % error in marriage duration for couples where both are in first marriage
hist(100*abs((mmardur-fmardur)[bfirst])/(.5*(mmardur+fmardur)[bfirst]), main ="% absolute difference in marriage duration listed by m & f", col = "black")
mardur.25 <- abs(mmardur-fmardur)/(.5*(mmardur+fmardur)) > .25
mardur.25[is.na(mardur.25)] <- F        # from 0 duration marriages
mardur.25 <- mardur.25 & bfirst         # only worried about error when both in first
mardur.25[is.na(mardur.25)] <- F
aggregate(mardur.25, list(raw$ds), function(x) signif(mean(x),3)*100)

## people married too young?
par(mfrow=c(2,1))
hist(fage-fmardur, breaks = 0:55, xlab = "age at first marriage", col = "black", main = "female")
hist(mage-mmardur, breaks = 0:85, xlab = "age at first marriage", col = "black", main = "male")
## Data to remove: people that say they married under 10 ys old if they are in their first marriage
## age in months at interview
raw$fage <- (raw$Finterview.cmc - raw$Fdob.cmc)
raw$mage <- (raw$Minterview.cmc - raw$Mdob.cmc)
mar.und.8 <- raw$mardur.mon >= raw$mage-8*12 |raw$mardur.mon >= raw$fage-8*12
mar.und.8[is.na(mar.und.8)] <- F
aggregate(mar.und.8, list(raw$ds), function(x) signif(mean(x),3)*100)

raw[mar.und.8,c("ds","mage","fage","tint","tmar")]
## Should we be worried that some people have been married longer than
## they've been having sex?
sum(raw$mardur > raw$fysa, na.rm=T)
sum(raw$mardur > raw$mysa, na.rm=T)
sum(raw$mardur > raw$fysa + 1, na.rm=T)
sum(raw$mardur > raw$mysa + 1, na.rm=T)
par(mfrow=c(1,2))
hist(raw$fysa - raw$mardur, main = "female", xlab = "", col="black")
abline(v=0, col = "red")
hist(raw$mysa - raw$mardur, main = "male", xlab="", col = "black")
abline(v=0, col = "red")
mtext("age at first marriage - age at first intercourse (yrs)", outer = T, side = 1, line = -2)
mar.bef.sex <- raw$fysa - fmardur < -1 | raw$mysa - mmardur < -1
mar.bef.sex[is.na(mar.bef.sex)] <- F
dev.off()
aggregate(mar.bef.sex, list(raw$ds), function(x) signif(mean(x),3)*100)

## no one should be having sex < 5 yrs old
early.sex <- raw$tint-raw$tms +60>= raw$mage
early.sex <- early.sex | raw$tint-raw$tfs +60>= raw$fage
early.sex[is.na(early.sex)] <- F
raw[early.sex,c("tint","tms","mage","tfs","fage")]
raw[early.sex,1:5]
aggregate(early.sex, list(raw$ds), function(x) signif(mean(x),3)*100)

## rows to remove
errs <- mar.und.8 | mardur.25 | mar.bef.sex | early.sex | no.aafi | int.gr1mon
ltab$errors <- xtabs(~group + errs, raw)[,2]
rem <- mar.und.8 | mardur.25 | mar.bef.sex | early.sex | no.aafi | int.gr1mon | noser | wives | !first.mar
## removed only because they didn't have serostatus
remser <- noser & !(mar.und.8 | mardur.25 | mar.bef.sex | early.sex | no.aafi | int.gr1mon | wives | !first.mar)
## chi sq for difference in  serostatus b/w those that were removed vs not
allrem <- rem
allraw <- data.frame(raw, remser, noser, rem)
allraw$ser[noser] <- NA
## Ethiopia 2005: "The survey was fielded from April 27 to August 30,
## 2005." (p. 11, Ethiopia DHS 2005 report)
## difference in dates is then
diff2005 <- 105*12 + 4 - min(allraw[allraw$ds=="Ethiopia 2005","tint"])
allraw[allraw$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] <- allraw[allraw$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] + diff2005
## Ethiopia 2011:"Allrawa collection took place over a five-month period
## from 27 December 2010 to 3 June 2011." (p. 10, Ethiopia DHS 2011
## report)
diff2011 <- 110*12 + 12 - min(allraw[allraw$ds=="Ethiopia 2011","tint"])
allraw[allraw$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] <- allraw[allraw$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] + diff2011

save(allraw, file = "alldhs raw.Rdata")   

xtabs(~ds + rem, raw)
aggregate(rem, list(raw$ds), function(x) signif(mean(x),3)*100)
ltab$analyzed <- xtabs(~group + rem, raw)[,1]
raw <- raw[!rem,]

names(raw)[names(raw)=="Mcircumcised"] <- "circ"
## remove any countries with less than 100 couples left
rem.count <- names(xtabs(~ds,raw))[xtabs(~ds,raw)<100]
rem <- raw$ds %in% rem.count | grepl("Sao Tome", raw$ds) # don't have epidemic curve for sao tome
oraw <- raw
raw <- raw[!rem,]

## Need to have ysa be >= mardur to avoid a nonnegative time b/w start
## of sex & start of marriage create serostatus of couple so for
## individuals where time of first sex was after marriage, set it to
## time of marriage.
apply(raw, 2, function(x) sum(is.na(x)))
aggregate(raw$tms>raw$tmar, list(raw$ds), function(x) signif(mean(x, na.rm=T),3)*100)
aggregate(raw$tfs>raw$tmar, list(raw$ds), function(x) signif(mean(x, na.rm=T),3)*100)
aggregate(raw$tms-raw$tmar, list(raw$ds), function(x) signif(range(x, na.rm=T),3)*100)
aggregate(raw$tfs-raw$tmar, list(raw$ds), function(x) signif(range(x, na.rm=T),3)*100)

nna <- !is.na(raw$tms) & !is.na(raw$tfs) & !is.na(raw$tmar)
raw$tms[nna & raw$tms>raw$tmar] <- raw$tmar[nna & raw$tms>raw$tmar]
raw$tfs[nna & raw$tfs>raw$tmar] <- raw$tmar[nna & raw$tfs>raw$tmar]
raw$fysa[nna & !raw$fysa>=raw$mardur] <- raw$mardur[nna & !raw$fysa>=raw$mardur]
raw$mysa[nna & !raw$mysa>=raw$mardur] <- raw$mardur[nna & !raw$mysa>=raw$mardur]

aggregate(raw$tms>raw$tmar, list(raw$ds), function(x) signif(mean(x, na.rm=T),3)*100)
aggregate(raw$tfs>raw$tmar, list(raw$ds), function(x) signif(mean(x, na.rm=T),3)*100)
aggregate(raw$tms-raw$tmar, list(raw$ds), function(x) signif(range(x, na.rm=T),3)*100)
aggregate(raw$tfs-raw$tmar, list(raw$ds), function(x) signif(range(x, na.rm=T),3)*100)

names(raw)[grepl('circum',names(raw))] <- 'circ'
## countries with circumcision are 2005/6 and later
rowSums(xtabs(~raw$ds + raw$circ))>0
tab8 <- xtabs(~raw$ds + raw$ser +raw$circ)[rowSums(xtabs(~raw$ds + raw$circ))>0,,]
tab8 <- rbind(tab8[,,2],tab8[,,3])
write.csv(tab8, file = "circum tab.csv")
######################################################################
names(raw)[names(raw)=="MHeard.of.drugs.to.help.infected.people.to.live.longer"] <- "m.k.arv"
names(raw)[names(raw)=="FHeard.of.drugs.to.help.infected.people.to.live.longer"] <- "f.k.arv"
names(raw)[names(raw)=="MLifetime.number.of.sexual.partners"] <- "mlsp"
names(raw)[names(raw)=="FLifetime.number.of.sexual.partners"] <- "flsp"
names(raw)[names(raw)=="MeverTestedHIVAIDS"] <- "mevtest"
names(raw)[names(raw)=="FeverTestedHIVAIDS"] <- "fevtest"
names(raw)[names(raw)=="Mnumber.of.unions"] <- "m.fun" # first union
names(raw)[names(raw)=="Fnumber.of.unions"] <- "f.fun"
levels(raw$m.fun) <- c(1,0,NA)           # 1 = first union, 0 not,
levels(raw$f.fun) <- c(1,0,NA)
## save loss filters

show <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","circ","mage","fage",
          "epic.ind","epic.nm","group", 'm.fun', 'f.fun', "m.k.arv","f.k.arv", "mlsp", "flsp", "mevtest","fevtest")
dat <- raw[,show]
dim(dat)

## add # of couples where both are in 1st cohabitation
both.first.mar <- raw$Fnumber.of.unions=="Once" & raw$Mnumber.of.unions=="Once"
ltab$bothfirstSCP <- xtabs(~group + both.first.mar,raw)[,2]
serstat <- xtabs(~group + ser, dat)[,c(4,2,3,1)]
ltab <- data.frame(ltab,serstat)
ltab
oltab <- ltab

for(ii in 1:nrow(ltab))
  {
    ltab[ii,8:12] <- paste(as.numeric(ltab[ii,8:12]),
                          " (", round(as.numeric(ltab[ii,8:12])/ltab[ii,7],3)*100,
                          "%)", sep="")
  }

for(ii in 1:nrow(ltab))
  {
    ltab[ii,3:7] <- paste(as.numeric(ltab[ii,3:7]),
                          " (", round(as.numeric(ltab[ii,3:7])/ltab[ii,2],3)*100,
                          "%)", sep="")
  }

write.csv(ltab, file = "loss table.csv")



pdf("wa epic.pdf", width = 8, height =6)
par(mfrow=c(1,2), mar = c(4,4,2,2))
plot(0,0, ylim = c(0,.35), xlim = c(1980,2012), type = "n", bty = "n",
     xlab = "year", ylab = "ART-normalized HIV prevalence in males")
for(ii in 1:length(wa.ind))
  {
    lines(1900 + 1/12*epicm$cmc, epicm[,wa.ind[ii]],
          col = rainbow(length(wa.ind))[ii])
  }
legend("topleft", names(epicm)[wa.ind], lwd = 1, col = rainbow(length(wa.ind)), bty = "n")
plot(0,0, ylim = c(0,.35), xlim = c(1980,2012), type = "n", bty = "n",
     xlab = "year", ylab = "ART-normalized HIV prevalence in females")
for(ii in 1:length(wa.ind))
  {
    lines(1900 + 1/12*epicm$cmc, epicf[,wa.ind[ii]],
          col = rainbow(length(wa.ind))[ii])
  }
dev.off()
######################################################################

######################################################################
## since it takes ~>1 month for antibodies to be detectable, for any
## couples that have been together for < 1 month, we set prev-yrs to 0
## for inside relationships
## dat[dat$tint-dat$tmar <= 1, c("m.dur.pm","f.dur.pm")] <- 0

dat$group <- factor(dat$group)

## Sexual partners per year of activity (s p acquisition rate: mspar, fspar)
dat$mspar <- dat$mlsp / (dat$tint - dat$tms + 1) * 12
dat$fspar <- dat$flsp / (dat$tint - dat$tfs + 1) * 12

pdf("SPAR plot.pdf")
par(mfrow = c(2,2))
breaks = seq(0, 200, by = .1)
hist(dat$mspar[dat$ser %in% c(3,4) & !is.na(dat$mspar)], breaks = breaks, main = "m-", col = "black", xlab = "SPAR",
     xlim = c(0,5))
hist(dat$fspar[dat$ser %in% c(2,4) & !is.na(dat$fspar)], breaks = breaks, main = "f-", col = "black", xlab = "SPAR",
     xlim = c(0,2))
hist(dat$mspar[dat$ser %in% c(1:2) & !is.na(dat$mspar)], breaks = breaks, main = "m+", col = "black", xlab = "SPAR",
     xlim = c(0,5))
hist(dat$fspar[dat$ser %in% c(1,3) & !is.na(dat$fspar)], breaks = breaks, main = "f+", col = "black", xlab = "SPAR",
     xlim = c(0, 2))
dev.off()

######################################################################

######################################################################
## Fix Ethiopian calendar which is off by 7-8 years from Gregorian calendar.

## Ethiopia 2005: "The survey was fielded from April 27 to August 30,
## 2005." (p. 11, Ethiopia DHS 2005 report)
## difference in dates is then
diff2005 <- 105*12 + 4 - min(dat[dat$ds=="Ethiopia 2005","tint"])
dat[dat$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] <- dat[dat$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] + diff2005

## Ethiopia 2011:"Data collection took place over a five-month period
## from 27 December 2010 to 3 June 2011." (p. 10, Ethiopia DHS 2011
## report)
diff2011 <- 110*12 + 12 - min(dat[dat$ds=="Ethiopia 2011","tint"])
dat[dat$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] <- dat[dat$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] + diff2011

save(dat,file= "alldhs.Rdata")
epicm <- as.matrix(epicm)
epicf <- as.matrix(epicf)
epicm.all <- as.matrix(epicm.all)
epicf.all <- as.matrix(epicf.all)
art.prev <- as.matrix(art.prev)

## if any epic index is negative set it to last value

save(epicm,file="allepicm.Rdata")
save(epicf,file="allepicf.Rdata")
save(epicm.all,file="epicmnart.Rdata")
save(epicf.all,file="epicfnart.Rdata")
save(art.prev,file="art.prev.Rdata")

######################################################################
write.csv(xtabs(~ds + ser, dat), file = "ser breakdown.csv")
######################################################################


######################################################################
## Prepare survival curve on monthly time intervals
######################################################################
pdf("surv by age for doc.pdf", width = 4.5, height = 3)
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
    csurv <- 1-cmort
    if(ii==1) plot(xseq/12, csurv, type = "l", xlim = c(0,25), col=cols[ii], lwd = 3, yaxt = "n",
         xlab = "years since seroconversion", ylab="probability of survival", bty = "n")
    if(ii!=1) lines(xseq/12, csurv, type = "l", col=cols[ii], lwd = 3)
  }
legend("topright",paste(age.seq/12,"yrs old"), col = cols, pch = 15, bty = "n", cex = 1, title = "age at seroconversion")
axis(2, seq(0,1,l=5), las = 2)
dev.off()

## Create a matrix of every month-age, cumulative survival probability to use in model
shp <- 2.3
ages <- 1:800                           # age in months
for(aa in ages)
  {
    if(aa==1) csurv <- 1 - pweibull(xseq, shape = shp, scale = 2000/shp/rep((aa/12)^.53, length(xseq)))
    if(aa>1)  csurv <- rbind(csurv, 1 - pweibull(xseq, shape = shp, scale = 2000/shp/rep((aa/12)^.53, length(xseq))))
  }

## So the cumulative probability a newly infected 20 yr old survives 8 yrs is
csurv[12*20,8*12]
## So the cumulative probability a newly infected 50 yr old survives 8 yrs is
csurv[12*50,8*12]

dim(csurv)

save(csurv, file = "csurv.Rdata")

seros <- dat$ser!=4
xtabs(~m.k.arv + group + seros,dat)
xtabs(~f.k.arv + group + seros,dat)
xtabs(~mevtest + group + ser,dat) / xtabs(~group + ser,dat)
xtabs(~fevtest + group + ser,dat)
xtabs(~circ + group + seros,dat)
names(dat)

## % ever tested and negative
test.tab <- round(xtabs(~mevtest + group,dat, subset = ser==4)[2,] / colSums(xtabs(~mevtest + group,dat, subset = ser %in% 4)),2)
test.tab <- rbind(test.tab, round(xtabs(~fevtest + group,dat, subset = ser==4)[2,] / colSums(xtabs(~fevtest + group,dat, subset = ser==4)),2))
## % ever tested given positive in discodant couple
test.tab <- rbind(test.tab, round(xtabs(~mevtest + group,dat, subset = ser==2)[2,] / colSums(xtabs(~mevtest + group,dat, subset = ser==2)),2))
test.tab <- rbind(test.tab, round(xtabs(~fevtest + group,dat, subset = ser==3)[2,] / colSums(xtabs(~fevtest + group,dat, subset = ser==3)),2))
## % ever tested given negative in discordant couple
test.tab <- rbind(test.tab, round(xtabs(~mevtest + group,dat, subset = ser==3)[2,] / colSums(xtabs(~mevtest + group,dat, subset = ser==3)),2))
test.tab <- rbind(test.tab, round(xtabs(~fevtest + group,dat, subset = ser==2)[2,] / colSums(xtabs(~fevtest + group,dat, subset = ser==2)),2))
## % ever tested given concordant positive
test.tab <- rbind(test.tab, round(xtabs(~mevtest + group,dat, subset = ser==1)[2,] / colSums(xtabs(~mevtest + group,dat, subset = ser==1)),2))
test.tab <- rbind(test.tab, round(xtabs(~fevtest + group,dat, subset = ser==1)[2,] / colSums(xtabs(~fevtest + group,dat, subset = ser==1)),2))

rownames(test.tab) <- c("M--", "F--", "M+-", "F-+", "M-+", "F+-", "M++", "F++")
test.tab <- t(test.tab)

test.tab <- test.tab[,c("M--", "M+-", "M-+", "M++", "F--",  "F-+",  "F+-",  "F++")]
write.csv(t(test.tab), file = "testing table.csv")

test.tab


mtest <- xtabs(~group + mevtest, dat) / rowSums(xtabs(~group + mevtest, dat))
ftest <- xtabs(~group + fevtest, dat) / rowSums(xtabs(~group + fevtest, dat))
gentest <- cbind(mtest[,2],ftest[,2])
colnames(gentest) <- c("male", "female")
write.csv(gentest, file = "gentest.csv")
