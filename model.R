##############################################################
## Clear workspace and load packages                        ##
##############################################################
# Clear workspace variables
rm(list = ls())
cat("\014")

#install.packages("survival")
#install.packages("ggplot2")
library(survival)
library(ggplot2)

##############################################################
## 1. Set working directory and get data                    ##
##############################################################
setwd("~/Dropbox/Documents/Projects/DataScience/SurvivalELSA")
load("elsa.RData")
attach(elsa_cf)

##############################################################
## 2. Getting a sense of the data                           ##
##############################################################
# Get list of variables in the dataset 'elsa_cf'
names(elsa_cf)

# Get summary statistic for 'age1' and 'time'
summary(elsa_cf$age1)
summary(elsa_cf$time)

# Histogram separately for people with an event and survivors #

#Create a subset#
elsadead <-subset(elsa_cf, elsa_cf$death==1)
elsalive <-subset(elsa_cf, elsa_cf$death==0)

summary(elsadead$time)
summary(elsalive$time)
# Histogram in each subset#
qplot(elsadead$time, geom = 'histogram', bins = 40)
qplot(elsalive$time, geom = 'histogram', bins = 40)

# Tabulate the outcome 'death'
table(elsa_cf$death)
prop.table(table(elsa_cf$death))

##############################################################
## 3. Survival and KM curve                                 ##
##############################################################
# Calculate person-years and incidence rate #
pyear <- sum(elsa_cf$time)
incidencerate <- 1802/pyear
print(pyear)
print(incidencerate)

# This gives a summary of survival time again
summary(Surv(time, death))

# Create a survival object and return summary #
km <- survfit(Surv(time, death)~ 1)
km 
# It doesn't give any value for the median survival


# Let's assess this graphically
# Draw the KM plot with CI
plot (km, lty=1, lwd=2, xlab="Time", ylab="Survival Probability", 
      col=rainbow(1))
# If we don't want Confidence limits
# We have to specify it in survfit
kmnoci <- survfit(Surv(time, death)~ 1, conf.type="none")
plot (kmnoci,lty=1, lwd=2, xlab="Time", ylab="Survival Probability", 
      col=rainbow(1))
# Quantiles of survival #
quantile(kmnoci)
# No value for the 25% either

# We can return the survival table at every time
# or specify at a specific time point, here 11y
# which is the end of follow-up
summary(kmnoci, times=11)
# We see here that at the end of follow-up
# the survival was ~80% hence less than 25%
# developed the event, which is why R can't estimate
# a 25% survival
##############################################################
## 4. Group comparisons                                     ##
##############################################################
# Draw the KM plot by sex
kmsex <- survfit( Surv(time, death)~ strata(sex), data=elsa_cf, 
                  conf.type="none")
plot(kmsex, lty=1, lwd=1, xlab="Time", ylab="Survival Probability", 
     col=rainbow(2)) 
legend("bottomleft", c("Men", "Women") , lty=1, lwd=1, 
       col=rainbow(2) )

#This function implements the G-rho family of Harrington and Fleming (1982), 
# with weights on each death of S(t)^rho, where S is the Kaplan-Meier 
#estimate of survival. 
# With rho = 0 this is the log-rank or Mantel-Haenszel test, 
# and with rho = 1 it is equivalent to the Peto & Peto modification
# of the Gehan-Wilcoxon test.
survdiff(Surv(time, death) ~ sex, rho=0)
# It is not significant: no difference between men and women
##############################################################
## 5. Exercise                                              ##
##############################################################
# Divide cognitive function score (tertiles)
quantile(elsa_cf$cf1, prob=c(0.33, 0.66))
elsa_cf$t_cf1 <-cut(elsa_cf$cf1, breaks=c(0, 42, 52, 194))
# Divide cognitive function score (quintiles)
quantile(elsa_cf$cf1, prob=c(0.20, 0.40, 0.60, 0.80))
elsa_cf$q_cf1 <-cut(elsa_cf$cf1, breaks=c(0, 38, 45, 50, 57, 194))

# Draw the KM plot by tertiles t_cf1
kmcft <- survfit( Surv(time, death)~ strata(t_cf1), data=elsa_cf, 
                  conf.type="none")
plot(kmcft, lty=1, lwd=2, xlab="Time", ylab="Survival Probability", 
     col=rainbow(3, start=0.1, end=0.9) ) 
legend("bottomleft", inset=.05, title="KM per cognitive function tertiles",
       c("T1", "T2", "T3"), lty=1, lwd=2 ,col = rainbow(3, start=0.1, end=0.9) )
#logrank
survdiff(Surv(elsa_cf$time, elsa_cf$death) ~ elsa_cf$t_cf1, rho=0)
# Significant difference

# Draw the KM plot by quintiles q_cf1
kmcfq <- survfit( Surv(time, death)~ strata(q_cf1), data=elsa_cf, 
                  conf.type="none")
plot(kmcfq, lty=1,lwd=2, xlab="Time", ylab="Survival Probability", 
     col=rainbow(5,start=0.1, end=0.9) ) 
legend("bottomleft", inset=.05, title="KM per cognitive function quintiles",
       c("Q1", "Q2", "Q3", "Q4", "Q5"), lty=1, lwd=2, 
       col = rainbow(5,start=0.1, end=0.9) )
#logrank
survdiff(Surv(elsa_cf$time, elsa_cf$death) ~ elsa_cf$q_cf1, rho=0)
# Significant difference
