# Calculate the predicted values for the mortality model, at the modelled levels of the random effects,
# and create a binned plot of the residuals.
# Sean Tuck
# 05/09/2013

library(lme4)
library(arm)

pred.lmer.binom <- function(model, X, Zt){
  if(missing(X)){
    X <- model@X
  }
  b <- fixef(model)
  if(missing(Zt)){
    Zt <- as.matrix(model@Zt)
  }
  z <- unlist(ranef(model))
  plogis(X %*% b + t(Zt) %*% z)
}

mort.dat<- read.table('PhD/Borneo Mortality/MaluaDataI_EditMar13.txt', header=TRUE, sep='\t') # Edited data file, see query 1 below.

mort.dat$pl.li<- factor(paste(mort.dat$pl,'.',mort.dat$li,sep=''))
mort.dat$pl<- as.factor(mort.dat$pl)
mort.dat$bl<- as.factor(mort.dat$bl)

mort.model <- lmer(survival ~ species + sd + plantingdate + (1|bl:pl) + (1|pl/pl.li), data=mort.dat, family=binomial, REML=FALSE)

preds <- pred.lmer.binom(mort.model)
resids <- resid(mort.model)
binnedplot(preds,resids)