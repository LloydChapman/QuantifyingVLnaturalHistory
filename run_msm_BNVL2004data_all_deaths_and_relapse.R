######## RUN CODE TO PROCESS DATA INTO RIGHT FORMAT AND LOAD MSM LIBRARY ########
source('process_BNVL2004data_all_deaths_and_relapse.R')
source('plot_survival_prob.R')
source('plot_obs_num_in_each_state.R')
source('plot_prevalence.R')
library(msm)

## Make table of total numbers of transitions observed between each pair of states
transtns_table<-statetable.msm(state,IDNUM,data_long)
print(transtns_table)

######## CALCULATE ROUGH INITIAL GUESS FOR TRANSITION INTENSITY MATRIX FROM OBSERVED NUMBERS OF TRANSITIONS ########
# Specify structure of transition intensity matrix
Qinit<-rbind(c(1,1,0,0,1),c(0,1,1,1,1),c(0,0,1,1,1),c(1,0,1,1,1),c(0,0,0,0,0))
# Calculate an initial guess for Q by assuming that times in data are exact transition times
Qinit<-crudeinits.msm(state~time,IDNUM,qmatrix=Qinit,data=data_long,censor=c(91,92,93),censor.states=list(c(1,4),c(2,4),c(1,2)))
 
# Guess based on trnstn intnsties for max lklhd from big parameter sweep
# Qinit<-rbind(c(-0.21,0.21,0,0,0.0055),c(0,-2.5,0.9,1.6,0.017),c(0,0,-2.9,2.9,0.09),c(0.32,0,0.01,-0.32,0.006),c(0,0,0,0,0))

print(Qinit)

######## RUN MSM METHOD TO FIT MODEL TO DATA (FIND MAX LIKELIHOOD) AND ESTIMATE TRANSITION INTENSITY MATRIX ########
# Change optimisation function using opt.method="...". Options are "optim", "nlm" or "bobyqa"
# optim options are: Nelder-Mead (default, gives different answer), BFGS, CG, L-BFGS-B, SANN (doesn't work)
# Change optimisation algorithm using method="...". See page 30 of msm 2014 manual for options. 
# Without covariates
# with BFGS method
system.time(BNVL2004.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,censor=c(91,92,93),method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=100,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))#,fixedpars=c(2,5,10)
# BNVL2004.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,censor=c(91,92,93),method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=200,fnscale=3500,reltol=1e-20))
# with CG method
# system.time(BNVL2004.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))

# With covariates
# with BFGS method
system.time(BNVL2004.age.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~age_groups,constraint=list(age_groupsB=c(1,2,3,4,2,5,6,7,8,2),age_groupsC=c(1,2,3,4,2,5,6,7,8,2)),censor=c(91,92,93),method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004.age.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=list("1-2"=~age_groups,"2-3"=~age_groups,"2-4"=~age_groups,"3-4"=~age_groups,"4-1"=~age_groups,"4-3"=~age_groups),censor=c(91,92,93), method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
system.time(BNVL2004.sex.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~SEX,censor=c(91,92,93), method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004.sex.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=list("1-2"=~SEX,"2-3"=~SEX,"2-4"=~SEX,"3-4"=~SEX,"4-1"=~SEX,"4-3"=~SEX),censor=c(91,92,93), method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
system.time(BNVL2004.bednet.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=list("1-2"=~Netsummer,"2-3"=~Netsummer,"2-4"=~Netsummer,"3-4"=~Netsummer,"4-1"=~Netsummer,"4-3"=~Netsummer),censor=c(91,92,93), method="BFGS",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))

# with CG method
# system.time(BNVL2004.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~age_groups,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3430,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=~SEX,censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))
# system.time(BNVL2004.msm<-msm(state~time,subject=IDNUM,data=data_long,qmatrix=Qinit,obstype=obstype,covariates=list("1-2"=~Netsummer,"2-3"=~Netsummer,"2-4"=~Netsummer,"3-4"=~Netsummer,"4-1"=~Netsummer,"4-3"=~Netsummer),censor=c(91,92,93), method="CG",censor.states=list(c(1,4),c(2,4),c(1,2)),control=list(maxit=1000,fnscale=3500,reltol=1e-20,trace=1,REPORT=1)))

# Print out transition intensities
# print(BNVL2004.msm)
# print(BNVL2004.age.msm)
# print(BNVL2004.sex.msm)
# print(BNVL2004.bednet.msm)

# Make survival plots for different states
# plot_survival_prob(BNVL2004.msm,legend.pos=c(2.4,0.98),xlab="Time (yrs)",ylim=c(0.9,1),lwd=2)

# Calculate transition probability matrix for time interval of 3 years
# P<-pmatrix.msm(BNVL2004.msm,t=3,ci="normal",cl=0.95)
# print(P)

## Calculate mean times spent in each state and convert into days
sjn_times<-sojourn.msm(BNVL2004.msm)
sjn_times_days<-365*sjn_times
print(sjn_times_days)

######## CALCULATE PROBABILITY OF DEVELOPING KA FROM ASYMPTOMATIC INFECTION ########
## Probability assuming survival
prob_symptms_no_death_fn<-function(x,covrtes="mean")
{
  qmatrix.obj<-qmatrix.msm(x,covariates=covrtes)
  Q<-qmatrix.obj$estimates
  prob_symptms_no_death<-Q[2,3]/(Q[2,3]+Q[2,4])
  return(prob_symptms_no_death)
}

# Without covariates
Q<-BNVL2004.msm$Qmatrices$baseline
print(Q)
prob_symptms_no_death<-prob_symptms_no_death_fn(BNVL2004.msm)
print(prob_symptms_no_death)

# With covariates
# Age
prob_symptms_0to14<-prob_symptms_no_death_fn(BNVL2004.age.msm,list(age_groupsB=0,age_groupsC=0))
prob_symptms_15to45<-prob_symptms_no_death_fn(BNVL2004.age.msm,list(age_groupsB=1,age_groupsC=0))
prob_symtpms_over45<-prob_symptms_no_death_fn(BNVL2004.age.msm,list(age_groupsB=0,age_groupsC=1))
print(prob_symptms_0to14)
print(prob_symptms_15to45)
print(prob_symtpms_over45)

# Sex
prob_symptms_male<-prob_symptms_no_death_fn(BNVL2004.sex.msm,list(SEX=1))
prob_symptms_female<-prob_symptms_no_death_fn(BNVL2004.sex.msm,list(SEX=2))
print(prob_symptms_male)
print(prob_symptms_female)

# Bed net use
prob_symptms_no_nets<-prob_symptms_no_death_fn(BNVL2004.bednet.msm,list(Netsummer=0))
prob_symptms_nets<-prob_symptms_no_death_fn(BNVL2004.bednet.msm,list(Netsummer=1))
print(prob_symptms_no_nets)
print(prob_symptms_nets)

# # Calculate confidence intervals for probability of developing KA assuming survival by bootstrapping
# prob_symptms_no_death.list<-boot.msm(BNVL2004.msm,stat=prob_symptms_no_death_fn,B=1000)
# prob_symptms_no_death_vec<-unlist(prob_symptms_no_death.list)
# std_dev_prob_symptms_no_death<-sd(prob_symptms_no_death_vec)
# CI_prob_symptms_no_death<-quantile(prob_symptms_no_death_vec,c(0.025,0.975))
# names(std_dev_prob_symptms_no_death)<-"se"
# print(c(std_dev_prob_symptms_no_death,CI_prob_symptms_no_death))

# ## Probability accounting for mortality
# prob_symptms<--qratio.msm(BNVL2004.msm,c(2,3),c(2,2),ci="none")
# print(prob_symptms)
# 
# # Calculate confidence intervals for probability of developing KA from asymptomatic infection using different methods 
# # delta method
# qr1<-qratio.msm(BNVL2004.msm,c(2,3),c(2,2))
# names(qr1)[3:4]<-c("U","L")
# prob_symptms_with_CI_delta_mthd<-c(-qr1[1],qr1[2],-qr1[4],-qr1[3])
# print(prob_symptms_with_CI_delta_mthd)
# # multivariate normal method
# qr2<-qratio.msm(BNVL2004.msm,c(2,3),c(2,2),ci="normal")
# names(qr2)[3:4]<-c("U","L")
# prob_symptms_with_CI_norm_mthd<-c(-qr2[1],qr2[2],-qr2[4],-qr2[3])
# print(prob_symptms_with_CI_norm_mthd)
# # bootstrap method
# qr3<-qratio.msm(BNVL2004.msm,c(2,3),c(2,2),ci="bootstrap")
# names(qr3)[3:4]<-c("U","L")
# prob_symptms_with_CI_btstrp<-c(-qr3[1],qr3[2],-qr3[4],-qr3[3])
# prob_symptms.list<-boot.msm(BNVL2004.msm,stat=function(x){-qratio.msm(x,c(2,3),c(2,2),ci="none")},B=100)
# prob_symptms_vec<-unlist(prob_symptms.list)
# std_dev_prob_symptms<-sd(prob_symptms_vec)
# CI_prob_symptms<-quantile(prob_symptms_vec,c(0.025,0.975))
# names(std_dev_prob_symptms)<-"se"
# print(c(std_dev_prob_symptms,CI_prob_symptms))

######## PERFORM MODEL ASSESSMENT ########
# Do likelihood ratio tests for adding each of the covariates to the model
LR1<-lrtest.msm(BNVL2004.msm,BNVL2004.age.msm)
LR2<-lrtest.msm(BNVL2004.msm,BNVL2004.sex.msm)
LR3<-lrtest.msm(BNVL2004.msm,BNVL2004.bednet.msm)
print(LR1)
print(LR2)
print(LR3)

# Calculate Akaike information criterion for each model with covariates
AICs<-AIC(BNVL2004.msm,BNVL2004.age.msm,BNVL2004.sex.msm,BNVL2004.bednet.msm)
print(AICs)

# Compare observed and expected prevalence 
# (N.B.: NOT APPROPRIATE FOR ASSESSING MODEL FIT AS WE DON'T KNOW NUMBERS INITIALLY IN EACH STATE)
BNVL2004prevalence.msm<-prevalence.msm(BNVL2004.msm,times=seq(3.2,5,0.25),initstates=c(858,149,5,591,0))
print(BNVL2004prevalence.msm)
# plot.prevalence.msm(BNVL2004.msm,mintime=3.2,maxtime=5,initstates=c(858,300,5,591,0),covariates="mean")
plot_prevalence(BNVL2004.msm,mintime=3.2,maxtime=5,initstates=c(858,149,5,591,0),covariates="mean",legend.pos=c(2002.5,10),startyr=1999)

# Plot observed numbers of people in each state
frame()
par(mfrow=c(1,1))
plot_obs_num_in_each_state(BNVL2004.msm,mintime=3.2,maxtime=5,startyr=1999)


# Uncomment to save variables in workspace to msmBNVL2004_all_deaths_and_relapse.RData
# save.image("msmBNVL2004_all_deaths_and_relapse.RData")
