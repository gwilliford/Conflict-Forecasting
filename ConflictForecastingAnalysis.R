setwd("C:/Users/gwill/Dropbox/Research/Bayes Forecasting/Analysis - Cluster")
library(readstata13)
library(caret)
library(R2jags)
library(ROCR)
library(doSNOW)
library(foreach)
library(xtable)
library(countrycode)
source("C:/Users/gwill/Dropbox/QL Folder/Bayes - MCMCTab.R")

# Load data and subset to relevant variables------------------------------------
dat.load <- read.dta13("ClusterDataMerged.dta")
dat.subs <- dat.load[dat.load$west == 1, ]

# Standardize
demean <- function (x) {
 x <- x - mean(x, na.rm = T)
}
standardize <- function(x) {
  demean <- x - mean(x, na.rm = T)
  std <- demean / (2 * sd(demean, na.rm = T))
  x <- std
}
dat.stand <- dat.subs
attach(dat.subs)
dat.stand$caprat <- standardize(caprat)
dat.stand$caprat2 <- standardize(caprat2)
dat.stand$pchcap <- standardize(pchcap)
dat.stand$mpdy <- demean(mpdy)
dat.stand$nukdy <- demean(nukdy)
dat.stand$polity21 <- standardize(polity21)
dat.stand$polity22 <- standardize(polity22)
dat.stand$politym <- standardize(politym)
dat.stand$defense <- demean(defense)
dat.stand$igosum <- standardize(igosum)
dat.stand$prevmids <- standardize(prevmids)
dat.stand$ccdistance <- standardize(ccdistance)
dat.stand$trival <- demean(trival)
dat.stand$terr <- demean(terr)
detach(dat.subs)
dat.stand <- na.omit(dat.stand)

# Split into training and test sets
dat.stand$test <- ifelse(dat.stand$year < 1990, 0, 1)
dat.train <- dat.stand[dat.stand$test == 0, ]
dat.test  <- dat.stand[dat.stand$test == 1, ]

# Undersample training data
set.seed(19910321)
dat.undersamp <- downSample(dat.train, as.factor(dat.train$midonset))

# Set up simulated predicted probabilities
capexp <- seq(-.6, 1, .1)
pchexp <- seq(-3, 7, 1)
attach(dat.stand)
ppmat.caprat <- expand.grid(1, capexp, 1, median(pchcap), median(mpdy), median(polity21), median(polity22), median(politym), median(defense), median(igosum), median(prevmids), median(nukdy))#, median(trival), median(terr))
ppmat.pchcap <- expand.grid(1, median(caprat), median(caprat2), pchexp, median(mpdy), median(polity21), median(polity22), median(politym), median(defense), median(igosum),  median(prevmids), median(nukdy))#, median(trival), median(terr))
varnamevec <- c("Constant", "Capability Ratio", "Capability Ratio Squared", "Percent Change Capabilities", "Major Power", "Polity A", "Polity B", "Polity A*B", "Defensive Alliance", "Intergovernmental Organizations", "Previous MIDs", "Nuclear Power")#, "Rivalry", "Territorial Claim")
colnames(ppmat.caprat) <- varnamevec
colnames(ppmat.pchcap) <- varnamevec
ppmat.caprat[, 3] <- ppmat.caprat[, 2]^2
detach(dat.stand)

# Cauchy function
logit.cauchy <- function () {
   for (i in 1:nrow) {
     y[i] ~ dbern(p[i])
     logit(p[i]) <- inprod(X[i,],b)
   }
   b[1] ~ dt(0, pow(10, -2), 1)
   for (i in 2:ncov) {
     b[i] ~ dt(0, pow(2.5, -2), 1)
   }
}

# Function to adjust predicted probabilities for undersampling
ppcalib <- function(pp, cgamma) {
  (cgamma * pp) / (cgamma * pp - pp + 1)
}
cgamma = sum(dat.undersamp$midonset == 0)/sum(dat.train$midonset == 0)

# Create X matrix for estimation
mform <- as.formula(midonset ~ caprat + caprat2 + pchcap + mpdy + polity21 + polity22 + politym + defense + igosum + prevmids + nukdy)
Xtemp1 <- model.matrix(mform, data = dat.undersamp)
colnames(Xtemp1) <- varnamevec
Xtrain <- rbind(Xtemp1, ppmat.caprat, ppmat.pchcap)
y.us <- dat.undersamp[, "midonset"]
y.full <- c(y.us, rep(NA, nrow(ppmat.caprat)), rep(NA, nrow(ppmat.pchcap)))
datfin <- list(X = Xtrain, nrow = nrow(Xtrain), ncov = ncol(Xtrain), y = y.full)

# # X Test
# Xtest <- model.matrix(mform, data = dat.train)
# colnames(Xtrain) <- varnamevec
# colnames(Xtest) <- varnamevec

# Run models
m.cauchy <- jags.parallel(data = datfin,
  parameters.to.save = c("b", "p"),
  model = logit.cauchy,
  n.chains = 2,
  n.iter = 100000,
  n.thin = 1
)

mtab <- mcmctab(m.cauchy$BUGSoutput)
mtab <- as.matrix(mtab[1:datfin$ncov, 1:4]); mtab
rownames(mtab) <- varnamevec
colnames(mtab) <- c("Mean","Std. Dev.","2.5 Percent","97.5 Percent")
print(xtable(mtab, align=c("l","c","c","c","c"),label="tab:results"),type="latex",file="./tab_results.tex",sanitize.text.function=function(x){x})


#quantfun <- function(x) quantile(x, c(.05, .95))
#pp.results.med <- m.cauchy$BUGSoutput$median$p
m.cauchy.mcmc <- as.mcmc(m.cauchy)
m.cauchy.mcmclist <- as.mcmc.list(m.cauchy.mcmc)
pp.res.full <- as.data.frame(summary(m.cauchy.mcmc)$quantiles[(datfin$ncov + 2):(datfin$ncov + datfin$nrow + 1), c(1, 3, 5)])
pp.res.full$order <- substring(row.names(pp.res.full), 3, nchar(row.names(pp.res.full))-1)
pp.res.full <- pp.res.full[order(as.numeric(pp.res.full$order)), ]
pp.res.full <- pp.res.full[, 1:3]
pp.cauchy.us <- pp.res.full[1:nrow(dat.undersamp), ]
pp.caprat.st <- nrow(dat.undersamp) + 1
pp.caprat.en <- pp.caprat.st + nrow(ppmat.caprat) - 1
pp.pchcap.st <- pp.caprat.en + 1
pp.pchcap.en <- pp.pchcap.st + nrow(ppmat.pchcap) - 1
pp.res.caprat <- pp.res.full[pp.caprat.st:pp.caprat.en, ]
pp.res.pchcap <- pp.res.full[pp.pchcap.st:pp.pchcap.en, ]

matplot(capexp, pp.res.caprat[1:3], type = "lll", lty = c(2, 1, 2), col = 1)
matplot(pchexp, pp.res.pchcap[1:3], type = "lll", lty = c(2, 1, 2), col = 1)

# Prediction performance in test sample
b.cauchy <- m.cauchy$BUGSoutput$median$b
X.test <- model.matrix(mform, data = dat.test)
X.train.full <- model.matrix(mform, data = dat.train)

# Calculate os pps
xb.cauchy.test <- as.matrix(X.test) %*% b.cauchy
pp.cauchy.test.temp <- exp(xb.cauchy.test) / (1 + exp(xb.cauchy.test))
pp.cauchy.test <- pp.cauchy.test.temp # ppcalib(pp.cauchy.test.temp, cgamma)
#xb.cauchy.train.full <- as.matrix(X.train.full) %*% b.cauchy
#pp.cauchy.train.full.temp <- exp(xb.cauchy.train.full) / (1 + exp(xb.cauchy.train.full))
#pp.cauchy.train.full <- pp.cauchy.test.temp #ppcalib(pp.cauchy.train.full.temp, cgamma)

# PVs
pv.cauchy.us <- ifelse(pp.cauchy.us[2] > .5, 1, 0)
pv.cauchy.test <- ifelse(pp.cauchy.test > .5, 1, 0)

# Confusion matrices
cmat.test <- table(dat.test$midonset, pv.cauchy.test)
cmat.test.perc <- round(prop.table(cmat.test, 1) * 100, 2)
for (i in 1:2) {
  cmat.test[i, 1] <- paste(cmat.test[i, 1], " (", cmat.test.perc[i, 1], "%)", sep = "")
  cmat.test[i, 2] <- paste(cmat.test[i, 2], " (", cmat.test.perc[i, 2], "%)", sep = "")
}
rownames(cmat.test) <- c("No Dispute Observed"," Dispute Observed")
colnames(cmat.test) <- c("No Dispute Predicted", "Dispute Predicted")
print(xtable(cmat.test, digits = 0), type = "latex", file = "./tab_insamp2x2.tex", booktabs = T)
cmat.test

cmat.us <- table(y.us, pv.cauchy.us)
cmat.us.perc <- round(prop.table(cmat.us, 1) * 100, 2)
for (i in 1:2) {
  cmat.us[i, 1] <- paste(cmat.us[i, 1], " (", cmat.us.perc[i, 1], "%)", sep = "")
  cmat.us[i, 2] <- paste(cmat.us[i, 2], " (", cmat.us.perc[i, 2], "%)", sep = "")
}
rownames(cmat.us)<-c("No Dispute Observed","Dispute Observed")
colnames(cmat.us)<-c("No Dispute Predicted","Dispute Predicted")
print(xtable(cmat.us, digits = 0), type = "latex", file = "./tab_os2x2.tex", booktabs = T)
cmat.us

# Prediction stats
brier.cauchy.test <- mean((pp.cauchy.test - dat.test$midonset)^2)
pred.cauchy.test <- ROCR::prediction(pp.cauchy.test, dat.test$midonset)
roc.cauchy.test <- ROCR::performance(pred.cauchy.test, "tpr", "fpr")
auc.cauchy.test <- unlist(ROCR::performance(pred.cauchy.test, "auc")@y.values)
auc.cauchy.test
brier.cauchy.test

brier.cauchy.us <- mean((pp.cauchy.us[, 2] - y.us))
pred.cauchy.us <- ROCR::prediction(pp.cauchy.us[2], y.us)
roc.cauchy.us <- ROCR::performance(pred.cauchy.us, "tpr", "fpr")
auc.cauchy.us <- unlist(ROCR::performance(pred.cauchy.us, "auc")@y.values)
auc.cauchy.us
brier.cauchy.us

# Plot ROC Curves
isleg <- c(paste("AUC =", round(auc.cauchy.us, 2)), paste("Brier Score =", round(brier.cauchy.us, 2)))
osleg <- c(paste("AUC =", round(auc.cauchy.test, 2)), paste("Brier Score =", round(brier.cauchy.test, 2)))

png("ROCis.png", width = 400, height = 400)
par(mar = c(4,4,1,1))
plot(roc.cauchy.us, asp = 1)
temp <- legend("bottomright", legend = c(" ", " "), text.width = strwidth("Brier Score = X"), xjust = 1, yjust = 1)
text(temp$rect$w + temp$rect$left, temp$text$y, isleg, pos = 2)
dev.off()

png("ROCos.png", width = 400, height = 400)
par(mar = c(4,4,1,1))
plot(roc.cauchy.test)
temp <- legend("bottomright", legend = c(" ", " "), text.width = strwidth("Brier Score = X"), xjust = 1, yjust = 1)
text(temp$rect$w + temp$rect$left, temp$text$y, osleg, pos = 2)
dev.off()

# Separation plots
sepplot <- function(obsdv, pp, filename, fill = "red") {
  df <- data.frame(obsdv = obsdv, pp = pp, rank =  rank(pp))
  df2 <- df[df$obsdv == 1, ]
  minx <- min(df2$rank)
  maxx <- max(df2$rank)
  ggsave(
    filename = filename,
    ggplot(aes(x = rank, y = obsdv), data = df) +
    geom_rect(aes(fill = obsdv, xmin = rank - 0.5, xmax = rank + .5, ymin = 0,
                    ymax = 1, group = rank), fill = fill, data = df2) +
    xlim (0, maxx) + geom_line(aes(x = rank, y = pp)) +
    labs(x = "Predicted Probability Rank", y = "Predicted Probability") +
    theme(plot.margin = unit(c(0,0,0,0), "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()),
    width = 6.5,
    height = 1.5,
    dpi = 1200
  )
}

sepplot(dat.undersamp$midonset, pp.caprat.us[, 2], "NDUinsampSepPlot.png")
sepplot(dat.train$midonset, pp.cauchy.train, "NDUoutsampSepPlot.png")
sepplot(dat.test$midonset, pp.cauchy, "NDUSepPlot.png")

# False negatives analysis
fnmat <- cbind(dat.test, pv.cauchy.test)
fnmat <- fnmat[fnmat[, "midonset"] == 1 & fnmat[, "pv.cauchy.test"] == 0, ]
dymid <-read.csv("gml-ndy-disputes-2.0.csv")
merger <- merge(fnmat, dymid, by.x = c("ccode1", "ccode2", "year"), by.y = c("ccode1", "ccode2", "year"), all.x = T)
fntab <- merger[, c("ccode1", "ccode2", "year", "dispnum", "fatality", "maxdur")]
fntab$ccode1 <- countrycode(fntab$ccode1, "cown", "country.name")
fntab$ccode2 <- countrycode(fntab$ccode2, "cown", "country.name")
fntab[1, 1] <- "Venezuela"
fntab[4, 2] <- "Venezuela"
fntab[5, 2] <- "Venezuela"
fntab <- fntab[order(fntab$ccode1, fntab$ccode2, fntab$year, fntab$dispnum, fntab$fatality, fntab$maxdur), ]
colnames(fntab)<-c("Country A", "Country B", "Year", "Mid Number", "Fatality", "Max Duration")
#fntab[c(1, 2, 4, 5), "Hostility Level"] <- "Use of Force"
#fntab[c(3, 6), "Hostility Level"] <- "Display of Force"
print(xtable(fntab, digits = 0),type = "latex", label = "tab_outsamp_fn", file = "tab_outsamp_fn.tex", include.rownames = F, booktabs = T, align = "r")
fntab
beep(8)

