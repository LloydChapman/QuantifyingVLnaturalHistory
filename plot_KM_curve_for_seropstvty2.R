rm(list=ls())
library(survival)
library(ggplot2)
source("hzd_ratio_with_conf_intvls.R")
source("createSurvivalFrame.R")

# Read in data
data<-read.csv("BNVL2004_augmntd_apprxmtd_death_dates2.csv",header=T,stringsAsFactors=F)

# Set ELISA cut-off for seropositivity and strong seropositivity
ELISA_pstve<-20
ELISA_strng_pstve<-61

# Convert dates into numeric serial dates
data[,c(10:12,16:18,42:43,53)] <- apply(data[,c(10:12,16:18,42:43,53)],2,function(x){as.numeric(as.Date(x,"%d/%m/%y"))})

# Filter out KA cases with onset before 2002 survey and individuals without an ELISA reading in 2002
no_ELISA02<-is.na(data$EIACU02)
data<-data[!no_ELISA02 & (data$ONSYR>=2002 | is.na(data$ONSYR)) & (data$KA==0 | data$FEV_ONS>data$INTDT02 | is.na(data$FEV_ONS>data$INTDT02)),]

# Calculate the time to onset for KA patients
time_to_onset<-data$FEV_ONS-data$INTDT02
# Calculate the maximum follow-up time for individuals without KA
time_to_follow_up<-apply(data[,16:18],1,function(x){max(x,na.rm=T)})-data$INTDT02

# Make a time vector of the times to events for KA and non-KA individuals
time<-numeric(length=nrow(data))
time<-apply(cbind(time_to_onset,time_to_follow_up),1,function(x){min(x,na.rm=T)})

# Make a vector of baseline (2002) serology status 
group<-numeric(length=nrow(data))
# Seronegative
group[data$EIACU02<ELISA_pstve]<-1
# Seropositive
group[data$EIACU02>ELISA_pstve & data$EIACU02<ELISA_strng_pstve]<-2
# Strong seropositive
group[data$EIACU02>ELISA_strng_pstve]<-3

# Combine data for Kaplan-Meier analysis
KMdata<-data.frame(IDNUM=data[,2],time,time_to_onset,time_to_follow_up,KA=data$KA)

# Create survival object
surv.obj<-Surv(time,KMdata$KA)

# Fit Cox prop hazard model with indicator variables for different groups to estimate hazard ratios with confidence intervals
ind2<-numeric(nrow(KMdata))
ind3<-numeric(nrow(KMdata))
ind2[group==2]<-1
ind3[group==3]<-1

coxph.obj<-coxph(surv.obj~ind2+ind3)
summary(coxph.obj)

# Perform log-rank test to determine if survival differs between any of serogroups
survdiff.obj<-survdiff(surv.obj~group,rho=0)
print(survdiff.obj)

# # Do log-rank tests for seropositive and strong seropositive groups against seronegative group
# gp12idx<-(group==1 | group==2)
# survdiff_ngtve_pstve.obj<-survdiff(Surv(time[gp12idx],KMdata$KA[gp12idx])~group[gp12idx],rho=0)
# print(survdiff_ngtve_pstve.obj)
# HR_and_95CI_ngtve_pstve<-hzd_ratio_with_conf_intvls(survdiff_ngtve_pstve.obj)
# print(HR_and_95CI_ngtve_pstve)
# 
# gp13idx<-(group==1 | group==3)
# survdiff_ngtve_strng_pstve.obj<-survdiff(Surv(time[gp13idx],KMdata$KA[gp13idx])~group[gp13idx],rho=0)
# print(survdiff_ngtve_strng_pstve.obj)
# HR_and_95CI_strng_pstve<-hzd_ratio_with_conf_intvls(survdiff_ngtve_strng_pstve.obj)
# print(HR_and_95CI_strng_pstve)

# coxph.obj<-coxph(surv.obj~group)

survfit.obj<-survfit(surv.obj~group,data=KMdata)
print(survfit.obj,print.rmean=T)
# plot(SurvFit,conf.int=TRUE,col=c("green","blue","red"),xlab="Follow-up time (days)",ylab="Proportion not progressing to KA")
# legend("bottom",c("Seronegative","Seropositive","Strong seropositive"),col=c("green","blue","red"),lty=1,inset=c(0.4,0.05))
# title("Kaplan-Meier curve for progression to KA \naccording to seropositivity")

surv.frame<-createSurvivalFrame(survfit.obj)
# ggplot(data=surv.frame, aes(colour = strata, group = strata)) + geom_step(aes(x = time, y = surv), direction = "hv") + geom_ribbon(aes(x = time,ymax = upper, ymin = lower, fill = strata), directions = "hv", linetype = 0, alpha = 0.25) + geom_point(data = subset(surv.frame, n.censor == 1), aes(x = time, y = surv), shape = 20) + theme_bw()
ggplot(data=surv.frame, aes(colour = strata, group = strata)) + geom_step(aes(x = time, y = surv), direction = "hv") + geom_step(aes(x = time, y = upper), directions = "hv", linetype = 2, alpha = 0.5) + geom_step(aes(x = time, y = lower), direction = "hv", linetype = 2, alpha = 0.5) + geom_point(data = subset(surv.frame, n.censor == 1), aes(x = time, y = surv), shape = 20) + theme_bw() + scale_colour_discrete(name="Serology status",breaks=c("group=1","group=2","group=3"),labels=c("negative","positive","strongly positive")) + xlab("Time since first survey (days)") + ylab("Proportion not progressing to KA") + xlim(0,800)+ylim(0,1)
