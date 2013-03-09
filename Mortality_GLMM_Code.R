# Code to explore and begin modelling Borneo mortality analysis (MaluaDataI.txt)
# Binomial GLMMs will be fitted to data, to estimate mortality of species and assess the relative importance
# of measured local environmental conditions in determining mortality.
# Read metadata (Malua_MetadataDataI.xlsx) file and Hector et al. 2011 for info on data
#
# Sean Tuck - PhD student, University of Oxford
# Last updated: 07/03/2012
##########
# Install packages for modelling
library(lme4)
library(arm)
##########
# Set working directory
setwd('/Users/sean/Dropbox')   # University computer
setwd('~/Dropbox')             # Home mac

# Read in mortality data
MortDat<- read.table('./SBE_survival/MaluaDataI.txt', header=TRUE, sep="\t")
MortDat<- read.table('PhD/Borneo Mortality/MaluaDataI_EDIT.txt', header=TRUE, sep='\t') # Edited data file, see query 1 below.

########## IMPORTANT
# 18 unique genus species combinations found (not incl. 'Dont plant'). Plots contain max. of 16 species.
# Unique species found when combining genus and species info that are not included in list on MaluaPlantingDesign.doc, are:
# 1. Dryobalanops tomentella; 2. Hopea faguetiana
# If this should be corrected, which component of species info is incorrect - species or genus?
# The following code would preserve erroneous lines's species data and make the genera data conform to relevant genera in 
# MaluaPlantingDesign.doc, and then create a new file MaluaDataI_EDIT.txt (rename MaluaDataI.txt if overwriting old file is needed):
# Uncomment section below to use correction
#
# Create genus sp. variable to give unique name for each species (rep of sp names in different genera)
# MortDat$genusspecies<- as.factor(paste(MortDat$genus,MortDat$species,sep=" "))
# MortDat$genus[MortDat$genusspecies == 'Dryobalanops tomentella']<- factor('Parashorea')
# MortDat$genus[MortDat$genusspecies == 'Hopea faguetiana']<- factor('Shorea')
# MortDat<- MortDat[,1:32] # Remove newly created genusspecies variable
# write.table(MortDat, file='MaluaDataI_EDIT.txt', sep='\t', row.names=FALSE)
#
# Also, the metadata file says "NA for all real missing values and unplanted points (1-9)". But the following code...
# unique(MortDat$survival[MortDat$unpntdpoint != 0])
# shows that there are 0s, not just NAs, at some unplanted points.
# Malua_MetadataDataI.xls shows 4 canopycover classes. In dataset, NA, 1-4, and 5 and 7 are found. Are 5 and 7 erroneous?
# NB: MaluaPlantingDesign.doc: gd:fd matrix shows one class in 4 gd with 2 fd (tall & medium III). The species listed underneath are all
# tall and the data show no 4 gd with 2 fd (see table(MortDat$gd,MortDat$fd)). Make change in PlantingDesign document?
#########################

#########################
#########################
# Data exploration
# What proportion of rows with NA survival?
sum(!is.na(MortDat$survival))/length(MortDat$survival) # 42% NAs
sum(!is.na(MortDat$survival)) # 62350 rows with survival data
# Of rows with survival data, what proportion are dead?
sum(MortDat$survival[which(!is.na(MortDat$survival))])/sum(!is.na(MortDat$survival)) # 36% recorded dead
sum(MortDat$survival[which(!is.na(MortDat$survival))]) # 22333 dead counts

# Quick look at survival by species
spsurvival<- table(MortDat$species, MortDat$survival)
barplot(c(spsurvival[,1],spsurvival[,2]), beside=TRUE, legend.text=c("Died","Survived"), ylab="Frequency", 
        main="Numbers by species of survived and dead specimens", axis.lty="solid", cex.names=0.4)
# Any species have higher frequency of survived than dead?
which(spsurvival[,2] > spsurvival[,1]) # Hopea sangal is the only plant with more observed survivals than deaths
spsurvival<- cbind(spsurvival,Proportion=spsurvival[,2]/(spsurvival[,1]+spsurvival[,2]))
summary(spsurvival[,3]) # Observed species-averaged survival range from 15-60%, with mean 35%

##########
# How many species under different planting conditions
table(MortDat$species,MortDat$sd) # TYPOs: Only 1 D. tomentella recorded, in 16sp plot (& it's NA); H. faguetiana only found in 1sp plots
par(mfrow=c(1,3))
for(i in 1:length(unique(MortDat$sd))) hist(table(MortDat$species,MortDat$sd)[,i],xlab=NULL,main=sort(unique(MortDat$sd))[i])
par(mfrow=c(1,1)) 
# Once typos and Dont plant are accounted for, 1sd and 16sd conditions have even spread of records for different species.

# Observed survival under different planting conditions (sd)
spsdsurvival<- table(MortDat$species,MortDat$sd,MortDat$survival)
spsdproportion<- spsdsurvival[,,2]/(spsdsurvival[,,1]+spsdsurvival[,,2])
boxplot(spsdproportion) # Looks to be no great change between sd classes, maybe other than slightly less variation at high sd

plot(as.factor(sort(unique(MortDat$sd))),NULL,ylim=c(0.1,0.7),type='n')
cols<- rainbow(nrow(spsdproportion))
for(i in 1:nrow(spsdproportion)) lines(as.factor(sort(unique(MortDat$sd))),spsdproportion[i,], col=cols[i])
# No real trend in by-species proportion survived with sd classes. Some increase with sd some decrease.
# If anything, many are lowest at intermediate sd.

##########
# Survival vs. damage types
table(MortDat$insectdamage, MortDat$survival)
table(MortDat$mammaldamage, MortDat$survival)
table(MortDat$treefalldamage, MortDat$survival)
# Not much increase in deaths as any type of damage increases in intensity.

##########
# Records for by-species by-plot survival
table(MortDat$species,MortDat$pl)
tapply(MortDat$pl,MortDat$species,summary)
# Max records found in a plot for all species is 124, and min for most is 1. Distribution similar for all species.

##########
# Planting times
table(MortDat$species,MortDat$timelag) # Amount of plants, for each species, for each timelag length (time between planting and survey)
par(mfrow=c(1,2))
hist(MortDat$plantdatenum)
hist(MortDat$surveydatenum)
par(mfrow=c(1,1))
# Boxplot of timelag (surveydate - plantingdate) conditional on survival
boxplot(MortDat$timelag ~ MortDat$survival, ylab="timelag",xlab="survival", horizontal=TRUE)
# Distribution of planting times for each genus
par(mfrow=c(2,3))
for(i in 1:length(levels(MortDat$genus))) hist(MortDat$plantdatenum[MortDat$genus == levels(MortDat$genus)[i]],xlab=NULL,main=levels(MortDat$genus)[i])
par(mfrow=c(1,1))
# Distribution of planting times for each species
tapply(MortDat$plantdatenum,MortDat$species, summary)
  
# NOTES:
# Based on data exploration, I would expect diversity (sd) to have little effect on estimated species mortality rates. Some species
# seem to be more likely to survive in mixture, and others when grown independently, but few responses seem dramatic. Observed
# proportion of survival differs greatly between species, so I expect estimated mortality rates to be different between species.
# Looks to be no significant impact of timelag (time between planting and survey) on survival.
# The number of species found under 1sd and 16sd conditions is relatively evenly spread. 4sd shows some greater variation, 
# mainly because 4 species (D. conformis, D. lanceolata, H. sangal, and H. spp.) occur much more regularly than the rest. At 4sd,
# the least occurring species is S. leprosula, at 638 records. So, at 4sd at least, there is some variation in how well each species
# is represented.
# The distributions of planting times between species are very similar - data were collected at similar times, so there is minimal
# species-associated bias or skew in collection dates (and by extension in any corresponding temporal-specific conditions).


#########################
#########################
# Modelling
# REML=FALSE used so models can be compared by AIC when appropriate.

# Null model - just random effects. Random effects specification, before adding fixed effects.
ranefmodpl<- lmer(survival ~ 1 + (1|pl), data=MortDat, family=binomial, REML=FALSE)
# It is important to include line level as a random effect as this is the funademental physical level of the experimental design.
# Line IDs are replicated within plots, however, so a plot.line ID needs to be created to incorporate line level in a relevant way.
MortDat$pl.li<- factor(paste(MortDat$pl,'.',MortDat$li,sep=''))
ranefmodli<- lmer(survival ~ 1 + (1|pl.li), data=MortDat, family=binomial, REML=FALSE)
AIC(ranefmodli)-AIC(ranefmodpl)
# Line level a much better fit.
# Line level included as a separate random effect. Not ideally how we'd want to specify it, as lines are physically
# nested within plots, nested within blocks.
# Need plot level as factor class for line nested within plot to work.
MortDat$pl<- as.factor(MortDat$pl)
ranefmodplli<- lmer(survival ~ 1 + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
# Include block level into the random effects structure as well, to see how much variance at the block level? How would this be
# specified - block/plot/line? Three nested grouping factor possible? Would it be easier/sensible to include block as a fixed effect?
# Keep just pl/pl.li for now.
# As specified by ranefmodplli, plot level and line level grouping factors capture a similar amount of variation in the data.

ranefmod2<- lmer(survival ~ 1 + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
anova(ranefmod,ranefmod2)
# Random intercept for forest types explains a lot of variation in the data and dramatically reduces deviance. Keep in model to focus
# on differences between species and effects of diversity on mortality rates.
##########

# Is foresttype random effect necessary, or is an appropriate fixed effect capturing similar information useful? Compare effects of
# canopycover fixed effect and foresttype random effect.
MortDat$canopycover<- as.factor(MortDat$canopycover)
spcnpymod<- lmer(survival ~ canopycover + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
spcnpymod2<- lmer(survival ~ canopycover + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
anova(ranefmodplli,spcnpymod)
anova(ranefmod2,spcnpymod2)
# Inclusion of canopycover categorical fixed effect is significant and improves AIC, but the coefficient values are rather wild.
# While they show large effects, when compared to the intercept as treatment constrasts, the effect is not large. Do show some small
# trend increasing with canopycover density however. Even with foresttype random effect, canopycover is improving the model fit.
# A problem, however, is correlation values of 1 and -1 between canopycover classes. This probably has something to do with the
# large standard error in parameter estimates, and rubbish p values. So, canopycover should not be included in the model, despite
# showing a better fit to the data and being a significant addition according to a likelihood ratio test.
##########

spmortmod<- lmer(survival ~ species + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
anova(ranefmod2,spmortmod)
# Likelihood Ratio Test shows inclusion of species fixed effect is justified. Species with most positive coefficients
# corroborates thoughts after data exploration on which species had highest proportion of survival. Likewise for those
# with lowest proportion of survival. Species coefficients within genera are not similar relative to among genera, so
# I think keeping data at resolution of species level is informative (want to estimate species mortality rates anyway).
##########

spsdmortmod<- lmer(survival ~ species + sd + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
anova(spmortmod2,spsdmortmod)
# Species diversity classes included as a fixed effect to see if increasing diversity significantly effects mortality.
# This does not improve the models explanatory power and shows diversity to have no effect on the mortality of a species.
spgdfdmod<- lmer(survival ~ species + gd*fd + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
spgdmortmod<- lmer(survival ~ species + gd + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
spfdmortmod<- lmer(survival ~ species + fd + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)
# No sign of significant effect of increasing generic diversity, functional diversity in terms of canopy mixture, or including the
# interaction between increasing functional diversity as generic diversity increases.
##########

# Include block as a fixed effect. Different response for different blocks?
spmortmodbl<- lmer(survival ~ species + bl + (1|foresttype) + (1|pl/pl.li), data=MortDat, family=binomial, REML=FALSE)

##########

# Models that include timelag and types of damage are not converging on meaningful values. Changing the approximation technique from
# Laplacian to Gauss-Hermite Quadrature (greater accuracy) does not solve false convergence. Leads to wild model estimates.
# Models that include any size covariates (continuous variables) also fail to converge...

##########

finalmod<- spmortmod # Current minimal adequate model is spmortmod
# Estimate mortality rates
fixef(finalmod)
mortpars<- vector(length=length(fixef(finalmod)))
mortpars[1]<- fixef(finalmod)[1]
mortpars[2:length(mortpars)]<- fixef(finalmod)[2:length(mortpars)]+mortpars[1] # Extract par estimates from treatment contrasts
mortrates<- matrix(c(levels(MortDat$species)[c(1,3:17)], 1/(1+exp(-mortpars))), ncol=2) 
# Back-transform logit values on to the original scale
# Object 'mortrates' stores the estimated mortality rates (percentage of survival) for all species.
# Proportion of survival ranges from 11% to 49% among species.

##########
# Model checking

# Use `binnedplot()` to look at residual diagnostics. Need to work out a way to predict from lmer GLMMs. Apparently examples 
# of methods to do this in Gelman & Hill, which Andy has.
# Use Nakagawa & Schielzeth's (2012) method to obtain Rsq value for models.
# Use cross-validation to see how well the par estimates generalise.
# Fit model to set of data randomly drawn from binomial distribution, instead of survival data, to make sure the
# estimates represent the real world, given the data.
