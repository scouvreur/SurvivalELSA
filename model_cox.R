##############################################################
## load packages                                            ##
##############################################################
# Clear workspace variables
rm(list = ls())
cat("\014")

# options(stringsAsFactors=F)
# install.packages("survival")
# install.packages("Formula")
# install.packages("ggplot2")
# install.packages("Hmisc")
# install.packages("lattice")
# install.packages("SparseM")
# install.packages("tinytex")
# tinytex::install_tinytex(force = TRUE)
library(survival)
library(Formula)
library(ggplot2)
library(Hmisc)
library(SparseM)
library(lattice)

# install.packages("sandwich")
library(sandwich)
# install.packages("rms")
library(rms)

#install.packages("devtools")
#library(devtools)
#install.version("rmarkdown",version=1.8)


##############################################################
## 1. Set working directory and get data                    ##
##############################################################
setwd("~/Dropbox/Documents/Projects/DataScience/SurvivalELSA")
# elsa_cf <- read.csv("ELSA_CF_TRUE.csv", sep=",")
# save(elsa_cf,file='elsa_cf.rdata')

load("elsa.RData")
attach(elsa_cf)

##############################################################
## 2. Cox regression models                                 ##
##############################################################

##############################################################
## 2.1. Simple models                                       ##
##############################################################

# Sex only
cox1 <- coxph(Surv(time, death)~ factor(sex), data=elsa_cf)
summary(cox1)
# Sex + age
cox2 <- coxph(Surv(time, death)~ age1 + factor(sex), data=elsa_cf)
summary(cox2)


# LR test: simply use anova of the two model fits
# Example to test the effect of adding CHD into the model
cox3 <- coxph(Surv(time, death)~ age1 + factor(sex) + factor(chd1), data=elsa_cf)
summary(cox3)
# likelihood ratio test:
anova(cox2,cox3)
# This returns a very small p-value, meaning that adding CHD to the model
# improves the likelihood. We can also simply see this in the p-value of the
# effect estimate of CHD in model 3.


##############################################################
## 2.2. COVARIATES : EXPLORATION AND SELECTION              ##
##############################################################
# Model 4 : all potential covariates included  #
cigst1_=factor(cigst1)
educ1_=factor(educ1)
cox4 <- coxph(Surv(time, death) ~ age1 + sex + chd1 + cancer1 +
              educ1_ + cigst1_ + physinact1 +
              alcohol1, data=elsa_cf, method="breslow")
summary(cox4)

### ALL COVARIATES ARE ASSOCIATED EXCEPT ALCOHOL ###

# Now that we have selected our set of covariates to include in the
# model, we want to test the crude and adjusted effect of cognitive
# function on mortality risk

#### MODEL 5 : COGNITIVE FUNCTION SCORE (CONTINUOUS)

cox5a <- coxph(Surv(time, death) ~ cf1 + age1 + sex,
               data=elsa_cf, method="breslow") # only adjusted for age and sex
cox5b <- coxph(Surv(time, death) ~ cf1 + age1 + sex +
               chd1 + cancer1+ educ1_ + cigst1_ +
               physinact1, data=elsa_cf, method="breslow") # adjusted for covariates
summary(cox5a)
summary(cox5b)

## In the multivariable-adjusted model (5b), the increase in 1 point of
## cognitive function score was associated with a 3% decrease in mortality risk
## Explanation: a hazard ratio of 0.97 means a 0.97-1=-0.03 *100=-3% 

## Note: all covariates remain significantly associated except for education

## We create a variable cf1/10 to interpret the HR as the reduction in risk
## for an increase of 10 points of score (instead of 1)
elsa_cf$cf1_10<-elsa_cf$cf1/10
summary(elsa_cf$cf1)
summary(elsa_cf$cf1_10)

cox5c <- coxph(Surv(time, death) ~ cf1_10 + age1 + sex +
               chd1 + cancer1+ educ1_ + cigst1_ + physinact1,
               data=elsa_cf, method="breslow") # adjusted for covariates
summary(cox5c)

## --> In the multivariable-adjusted model (5b), the increase in 10 points
## of cognitive function score was associated with a 28 % decrease in
## mortality risk

# Finally, we also look at it across quintiles
# Divide cognitive function score (quintiles)
quantile(elsa_cf$cf1, prob=c(0.20, 0.40, 0.60, 0.80))
elsa_cf$q_cf1 <-cut(elsa_cf$cf1, breaks=c(0, 38, 45, 50, 57, 194))
cox5c <- coxph(Surv(time, death) ~ factor(q_cf1) + age1 + sex +
               chd1 + cancer1+ educ1_ + cigst1_  +  physinact1,
               data=elsa_cf, method="breslow") # adjusted for covariates
summary(cox5c)

##############################################################
## 2.3. TESTING FOR INTERACTIONS                            ##
##############################################################

# LR test: simply use anova of the two model fits
# Example to test the effect of an interaction between age and cognitive function
## !!! THE MODELS HAVE TO BE NESTED
cox6a <- coxph(Surv(time, death)~ age1 + factor(sex) + cf1, data=elsa_cf)
cox6b <- coxph(Surv(time, death)~ age1 + factor(sex) + cf1 + age1*cf1, data=elsa_cf)
summary(cox6a)
summary(cox6b)
# likelihood ratio test:
anova(cox6a,cox6b)

## Very small p-value: evidence of a potential interaction between cognitive
## function and age, i.e. the effect of cognitive function on mortality
## differs according to ages

# Create agegr (age in categories) and stratify analysis
elsa_cf$agegr[elsa_cf$age1 < 60] <- 1
elsa_cf$agegr[elsa_cf$age1>=60 & elsa_cf$age1 < 70] <- 2
elsa_cf$agegr[elsa_cf$age1 >= 70] <- 3

# Create subsets #§
elsaage1 <-subset(elsa_cf, elsa_cf$agegr==1)
elsaage2 <-subset(elsa_cf, elsa_cf$agegr==2)
elsaage3 <-subset(elsa_cf, elsa_cf$agegr==3)

# Run models in each strata#
coxf.age1 <- coxph(Surv(time, death) ~ cf1_10 + age1 + sex +
                   chd1 + cancer1+ factor(educ1) + factor(cigst1) +
                   physinact1, data=elsaage1, method="breslow")
coxf.age2 <- coxph(Surv(time, death) ~ cf1_10 + age1 + sex +
                   chd1 + cancer1+ factor(educ1) + factor(cigst1) +
                   physinact1, data=elsaage2, method="breslow")
coxf.age3 <- coxph(Surv(time, death) ~ cf1_10 + age1 + sex +
                   chd1 + cancer1+ factor(educ1) + factor(cigst1)  +
                   physinact1, data=elsaage3, method="breslow")
summary(coxf.age1)
summary(coxf.age2)
summary(coxf.age3)

##############################################################
## 3. CHECKING PROPORTIONAL HAZARDS ASSUMPTION              ##
##############################################################

# Plot 1 - log (-log survival)
coxa <- cph(Surv(time, death) ~ cf1 + age1 + sex + chd1, 
            data=elsa_cf, x=TRUE,y=TRUE, method="breslow")

survplot(coxa, cf1=c(38,45,50,57), age1=mean(age1), sex=0, chd1=0,
         logt=TRUE, loglog=TRUE, xlim=c(-4, 2.5), ylim=c(-10, 0))

## Create a time-dependent model and test for interactions with time
## for all covariates
time.dep <- coxph (Surv(time, death)~ cf1_10 + age1  + sex + chd1 +
                   cancer1 + physinact1 + educ1_ + cigst1_ +  
                   strata(agegr), data=elsa_cf, method="breslow",
                   na.action=na.exclude)
time.dep.zph <- cox.zph(time.dep, transform = 'log')
time.dep.zph

# plots of residuals for all predictors (use arrows to go from one
# plot to the other)
plot(time.dep.zph)

# plots for the predictor number in bracket. In the model, cf1_10 is
# the first predictor so we are plotting cf1_10
plot(time.dep.zph[1])
abline(h=0, lty=3)
