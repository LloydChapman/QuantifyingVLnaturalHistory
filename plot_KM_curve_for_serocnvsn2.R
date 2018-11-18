rm(list=ls())
library(survival)
library(ggplot2)
source("pair_log_rank_test.R")
source("hzd_ratio_with_conf_intvls.R")
source("createSurvivalFrame.R")

# Read in data
data<-read.csv("BNVL2004_augmntd_apprxmtd_death_dates2.csv",header=T,stringsAsFactors=F)

# Set ELISA cut-off for seropositivity and strong seropositivity
ELISA_pstve<-20
ELISA_strng_pstve<-61

# Convert dates into numeric serial dates
data[,c(10:12,16:18,42:43,53)] <- apply(data[,c(10:12,16:18,42:43,53)],2,function(x){as.numeric(as.Date(x,"%d/%m/%y"))})

# Filter out KA cases with onset before 2002 survey and individuals without an ELISA reading in 2002 or 2003
no_ELISA02_or_no_ELISA03<-(is.na(data$EIACU02) | is.na(data$EIACU03))
no_KA<-(data$KA==0)
onset_after_srvy_st<-(data$ONSYR>=2002 & (data$FEV_ONS>data$INTDT02 | is.na(data$FEV_ONS>data$INTDT02)))
data<-data[!no_ELISA02_or_no_ELISA03 & (no_KA | onset_after_srvy_st),]

# Calculate the time to onset for KA patients
time_to_onset<-data$FEV_ONS-data$INTDT02
# Calculate the maximum follow-up time for individuals without KA
time_to_follow_up<-apply(data[,16:18],1,function(x){max(x,na.rm=T)})-data$INTDT02

# Make a time vector of the times to events for KA and non-KA individuals
time<-numeric(length=nrow(data))
time<-apply(cbind(time_to_onset,time_to_follow_up),1,function(x){min(x,na.rm=T)})

##########
# Make a vector for seroconversion 
group<-numeric(length=nrow(data))

seronegatives02<-(data$EIACU02<ELISA_pstve)
seropositives02<-(data$EIACU02>=ELISA_pstve & data$EIACU02<ELISA_strng_pstve)
serostrongpositives02<-(data$EIACU02>=ELISA_strng_pstve)
seronegatives03<-(data$EIACU03<ELISA_pstve)
seropositives03<-(data$EIACU03>=ELISA_pstve & data$EIACU03<ELISA_strng_pstve)
serostrongpositives03<-(data$EIACU03>=ELISA_strng_pstve)

# no seroconversion (remained seronegative or seropositive)
group[(seronegatives02 & seronegatives03) | (seropositives02 & seropositives03)]<-1
# remained strongly seropositive
group[serostrongpositives02 & serostrongpositives03]<-2
# from seropositive or strongly seropositive to seronegative (or strong seropositive to seropositive)
group[((seropositives02 | serostrongpositives02) & seronegatives03) | (serostrongpositives02 & seropositives03)]<-3
# from seronegative to seropositive
group[seronegatives02 & seropositives03]<-4
# from seronegative or seropositive to strong seropositive
group[(seronegatives02 | seropositives02) & serostrongpositives03]<-5
##########

# Combine data for Kaplan-Meier analysis
KMdata<-data.frame(IDNUM=data[,2],time,time_to_onset,time_to_follow_up,KA=data$KA)
# Check there are no individuals without a seroconversion status
print(summary(!group==0))
# KMdata<-KMdata[group!=0,]
# group<-group[group!=0]

# Create survival object
surv.obj<-Surv(KMdata$time,KMdata$KA)
survfit.obj<-survfit(surv.obj~group,data=KMdata)
# plot(survfit.obj,conf.int=TRUE,col=c("green","blue","red"),xlab="Follow-up time (days)",ylab="Proportion not progressing to KA")
# legend("bottom",c("Seronegative","Seropositive","Strong seropositive"),col=c("green","blue","red"),lty=1,inset=c(0.4,0.05))
# title("Kaplan-Meier curve for progression to KA \naccording to seropositivity")

# Fit Cox prop hazard model with indicator variables for different groups to estimate hazard ratios with confidence intervals
ind2<-numeric(nrow(KMdata))
ind3<-numeric(nrow(KMdata))
ind4<-numeric(nrow(KMdata))
ind5<-numeric(nrow(KMdata))
ind2[group==2]<-1
ind3[group==3]<-1
ind4[group==4]<-1
ind5[group==5]<-1

coxph.obj<-coxph(surv.obj~ind2+ind3+ind4+ind5)
summary(coxph.obj)

# # Perform log-rank test to determine if survival differs between any of serogroups
# survdiff.obj<-survdiff(surv.obj~group,rho=0)
# print(survdiff.obj)
# 
# # Do log-rank tests for different seroconversion groups against no seroconversion group
# # Remained strongly seropositive vs no conversion
# survdiff_none_vs_rmnd_strng_pstve.obj<-pair_log_rank_test(KMdata,group,1,2)
# print(survdiff_none_vs_rmnd_strng_pstve.obj)
# HR_and_95CI_none_vs_rmnd_strng_pstve<-hzd_ratio_with_conf_intvls(survdiff_none_vs_rmnd_strng_pstve.obj)
# print(HR_and_95CI_none_vs_rmnd_strng_pstve)
# 
# # Serodeconversion vs no conversion
# survdiff_none_vs_decnvsn.obj<-pair_log_rank_test(KMdata,group,1,3)
# print(survdiff_none_vs_decnvsn.obj)
# HR_and_95CI_none_vs_decnvsn<-hzd_ratio_with_conf_intvls(survdiff_none_vs_decnvsn.obj)
# print(HR_and_95CI_none_vs_decnvsn)
# 
# # Seroconversion against no conversion
# survdiff_none_vs_cnvsn.obj<-pair_log_rank_test(KMdata,group,1,4)
# print(survdiff_none_vs_cnvsn.obj)
# HR_and_95CI_none_vs_cnvsn<-hzd_ratio_with_conf_intvls(survdiff_none_vs_cnvsn.obj)
# print(HR_and_95CI_none_vs_cnvsn)
# 
# # Strong serocnvsn against no conversion
# survdiff_none_vs_strng_cnvsn.obj<-pair_log_rank_test(KMdata,group,1,5)
# print(survdiff_none_vs_strng_cnvsn.obj)
# HR_and_95CI_none_vs_strng_cnvsn<-hzd_ratio_with_conf_intvls(survdiff_none_vs_strng_cnvsn.obj)
# print(HR_and_95CI_none_vs_strng_cnvsn)
# 
# chisq <- matrix(0,5,5) 
# for (i in 1:5) 
# { 
#   for (j in (1:5)[-i]) 
#     { 
#     temp <- survdiff(Surv(time,KA)~group,data=KMdata,subset=(group %in% unique(group[order(group)])[c(i,j)])) 
#     chisq[i,j] <- temp$chisq
#     # print(temp)
#     } 
# }
# print(chisq)

surv.frame<-createSurvivalFrame(survfit.obj)
# ggplot(data=surv.frame, aes(colour = strata, group = strata)) + geom_step(aes(x = time, y = surv), direction = "hv") + geom_ribbon(aes(x = time,ymax = upper, ymin = lower, fill = strata), directions = "hv", linetype = 0, alpha = 0.25) + geom_point(data = subset(surv.frame, n.censor == 1), aes(x = time, y = surv), shape = 20) + theme_bw()
ggplot(data=surv.frame, aes(colour = strata, group = strata)) + geom_step(aes(x = time, y = surv), direction = "hv") + geom_step(aes(x = time, y = upper), directions = "hv", linetype = 2, alpha = 0.5) + geom_step(aes(x = time, y = lower), direction = "hv", linetype = 2, alpha = 0.5) + geom_point(data = subset(surv.frame, n.censor == 1), aes(x = time, y = surv), shape = 20) + theme_bw() + scale_colour_discrete(name="Seroconversion",breaks=c("group=1","group=3","group=4","group=2","group=5"),labels=c("positive/strongly \npositive to negative","none","negative to \npositive","strongly positive to \nstrongly positive","negative/positive to \nstrongly positive")) + xlab("Time since first survey (days)") + ylab("Proportion not progressing to KA")
