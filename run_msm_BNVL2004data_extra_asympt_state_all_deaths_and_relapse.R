######## RUN CODE TO PROCESS DATA INTO RIGHT FORMAT AND LOAD MSM LIBRARY ########
source('process_BNVL2004data_extra_asympt_state_all_deaths_and_relapse.R')
source('plot_survival_prob.R')
library(msm)

# Make table of total numbers of transitions observed between each pair of states
transtns_table<-statetable.msm(state,IDNUM,data_long)
print(transtns_table)
dim(data_long)

######## CALCULATE ROUGH INITIAL GUESS FOR TRANSITION INTENSITY MATRIX FROM OBSERVED NUMBERS OF TRANSITIONS ########
# Specify structure of transition intensity matrix
Qinit<-rbind(c(1,1,0,0,1,1),c(0,1,1,0,0,0),c(0,0,1,1,1,0),c(1,0,1,1,1,0),c(0,0,0,0,0,0),c(0,0,0,1,1,1))
# Calculate an initial guess for Q by assuming that times in data are exact transition times
Qinit<-crudeinits.msm(state~time,IDNUM,qmatrix=Qinit,data=data_long,censor=c(91,92,93),censor.states=list(c(1,4),c(2,4,6),c(1,2,6)))
print(Qinit)

######## RUN MSM METHOD TO FIT MODEL TO DATA (FIND MAX LIKELIHOOD) AND ESTIMATE TRANSITION INTENSITY MATRIX ########
# Change optimisation function using opt.method="...". Options are "optim", "nlm" or "bobyqa"
# optim options are: Nelder-Mead (default, gives different answer), BFGS, CG, L-BFGS-B, SANN (doesn't work)
# Change optimisation algorithm using method="...". See page 30 of msm 2014 manual for options. 
# Without covariates
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,censor=c(91,92,93),method="BFGS",censor.states=list(c(1,4),c(2,4,6),c(1,2,6)),control=list(maxit=100,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))

# With covariates
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~data_long$SEX,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~data_long$SEX,censor=c(91,92,93), method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~data_long$age_groups,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~data_long$age_groups,censor=c(91,92,93), method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~data_long$Netsummer,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates= list("1-2"=~data_long$age_groups,"2-3"=~data_long$age_groups,"2-4"=~data_long$age_groups,"3-4"=~data_long$age_groups,"4-1"=~data_long$age_groups),censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# 3-5 seems to be the problematic transition for an age covariate
# system.time(BNVL2004_extra_asympt_state.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates= list("1-2"=~data_long$age_m),censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))

# Print out transition intensities
print(BNVL2004_extra_asympt_state.msm)

# Make survival plots for different states
# plot_survival_prob(BNVL2004_extra_asympt_state.msm,legend.pos=c(2.4,0.98),xlab="Time (yrs)",ylim=c(0.9,1),lwd=2)

# Calculate transition probability matrix for time interval of 3 years
# P<-pmatrix.msm(BNVL2004_extra_asympt_state.msm,t=3,ci="normal",cl=0.95)
# print(P)

# Calculate mean times spent in each state and convert into days
sjn_times<-sojourn.msm(BNVL2004_extra_asympt_state.msm)
sjn_times_days<-365*sjn_times
print(sjn_times_days)

######## CALCULATE PROBABILITY OF DEVELOPING KA FROM INFECTION ########
## Probability assuming survival
prob_symptms_no_death_fn<-function(x,covrtes="mean")
{
  qmatrix.obj<-qmatrix.msm(x,covariates=covrtes)
  Q<-qmatrix.obj$estimates
  prob_symptms_no_death<-Q[1,2]/(Q[1,2]+Q[1,6])
  return(prob_symptms_no_death)
}

Q<-BNVL2004_extra_asympt_state.msm$Qmatrices$baseline
print(Q)
prob_symptms_no_death<-prob_symptms_no_death_fn(BNVL2004_extra_asympt_state.msm)
print(prob_symptms_no_death)

# # Calculate confidence intervals for probability of developing KA assuming survival by bootstrapping
# prob_symptms_no_death.list<-boot.msm(BNVL2004_extra_asympt_state.msm,stat=prob_symptms_no_death_fn,B=1000)
# prob_symptms_no_death.vec<-unlist(prob_symptms_no_death.list)
# prob_symptms_no_death.vecdble<-as.double(prob_symtpms_no_death.vec)
# std_dev_prob_symptms_no_death<-sd(prob_symptms_no_death.vecdble,na.rm=T)
# CI_prob_symptms_no_death<-quantile(prob_symptms_no_death.vecdble,c(0.025,0.975),na.rm=T)
# names(std_dev_prob_symptms_no_death)<-"se"
# print(c(std_dev_prob_symptms_no_death,CI_prob_symptms_no_death))
# 
# ## Probability accounting for mortality
# prob_symptms<--qratio.msm(BNVL2004_extra_asympt_state.msm,c(1,2),c(1,1),ci="none")
# print(prob_symptms)
# 
# # Calculate confidence intervals for probability of developing KA using different methods 
# # delta method
# qr1<-qratio.msm(BNVL2004_extra_asympt_state.msm,c(1,2),c(1,1))
# names(qr1)[3:4]<-c("U","L")
# prob_symptms_with_CI_delta_mthd<-c(-qr1[1],qr1[2],-qr1[4],-qr1[3])
# print(prob_symptms_with_CI_delta_mthd)
# # multivariate normal method
# qr2<-qratio.msm(BNVL2004_extra_asympt_state.msm,c(1,2),c(1,1),ci="normal")
# names(qr2)[3:4]<-c("U","L")
# prob_symptms_with_CI_norm_mthd<-c(-qr2[1],qr2[2],-qr2[4],-qr2[3])
# print(prob_symptms_with_CI_norm_mthd)
# # # bootstrap method
# # qr3<-qratio.msm(BNVL2004_extra_asympt_state.msm,c(1,2),c(1,1),ci="bootstrap")
# # names(qr3)[3:4]<-c("U","L")
# # prob_symptms_with_CI_btstrp<-c(-qr3[1],qr3[2],-qr3[4],-qr3[3])
# # prob_symptms.list<-boot.msm(BNVL2004_extra_asympt_state.msm,stat=function(x){-qratio.msm(x,c(2,3),c(2,2),ci="none")},B=100)
# # prob_symptms_vec<-unlist(prob_symptms.list)
# # std_dev_prob_symptms<-sd(prob_symptms_vec)
# # CI_prob_symptms<-quantile(prob_symptms_vec,c(0.025,0.975))
# # names(std_dev_prob_symptms)<-"se"
# # print(c(std_dev_prob_symptms,CI_prob_symptms))

######## PERFORM MODEL ASSESSMENT ########
# N.B. Can't do likelihood ratio test with 5-state model as models aren't nested - 
# some transitions in 5-state model are missing from 6-state model
AIC_extra_asympt_state<-AIC(BNVL2004_extra_asympt_state.msm)
print(AIC_extra_asympt_state)

# Uncomment to save variables in workspace to msmBNVL2004_extra_asympt_state_all_deaths_and_relapse.RData
# save.image("msmBNVL2004_extra_asympt_state_all_deaths_and_relapse.RData")
