rm(list=ls()) # clear all
## Set treatment period in days
trtmnt_prd<-20
symptm_prd<-365

## Set cut-offs for ELISA and LST positivity
ELISA_thr<-20
LST_thr<-5

## Read in data from csv file
data<-read.csv("BNVL2004_all_deaths_and_relapse.csv",header=T,stringsAsFactors=F)
## Remove people who have no ELISA or LST readings and no KA and no death
hasnodata<-(is.na(data$EIACU02) & is.na(data$EIACU03) & is.na(data$EIACU04 ) & is.na(data$MEANLST02) & is.na(data$MEANLST03) & is.na(data$MEANLST04) & data$KA==0 & is.na(data$NATDEATH))
data<-data[!hasnodata,]

## Convert dates into time since start of study (taken as 01/01/1999)
STARTDT="1999-01-01"
data[,c(10:12,16:18,42:43,53,57)] <- apply(data[,c(10:12,16:18,42:43,53,57)],2,function(x){as.numeric(as.Date(x,"%d/%m/%y")-as.Date(STARTDT))})

## Work out the age of each individual at time zero
age<-(data$AGE02-(data$INTDT02/365))
age[is.na(age)]<-(data$AGE03[is.na(age)]-(data$INTDT03[is.na(age)]/365))
age[is.na(age)]<-(data$AGE04[is.na(age)]-(data$INTDT04[is.na(age)]/365))
sum(is.na(age)) ### there is 1 person with 1 ELISA and 1 LST reading but no age information

######## CONVERT DATA TO LONG DATA FORMAT ########
## Make a data frame with the first survey observations
data1<-data[!is.na(data$INTDT02),] # select people who have a first survey date
data1<-data1[is.na(data1$KADEATH_DT) | (!is.na(data1$KADEATH_DT) & data1$INTDT02<=data1$KADEATH_DT),] # people who don't have a death date, or it is after the first survey

## Add first survey data
data1$timeindys<-data1$INTDT02
data1$ELISA_m<-data1$EIACU02
data1$LST_m<-data1$MEANLST02
data1$KA_m<-data1$KA_02
data1$death_m<-0
data1$timepoint<-"S1"
data1$age_m<-data1$AGE02
data_long<-data1
dim(data_long)

## Make a data frame with the second survey observations
data1<-data[!is.na(data$INTDT03),] # select people who have a second survey date
data1<-data1[is.na(data1$KADEATH_DT) | (!is.na(data1$KADEATH_DT) & data1$INTDT03<=data1$KADEATH_DT),] # people who don't have a death date, or it is after the second survey
data1$timeindys<-data1$INTDT03
data1$ELISA_m<-data1$EIACU03
data1$LST_m<-data1$MEANLST03
data1$KA_m<-data1$KA_03
data1$age_m<-data1$AGE03
data1$death_m<-0
data1$timepoint<-"S2"
data_long<-rbind(data_long,data1)
dim(data_long)

## Make a data frame with the third survey observations, ignore interviews which have a date but are after someone died
data1<-data[!is.na(data$INTDT04),] # those people who have a third survey date
data1<-data1[is.na(data1$KADEATH_DT) | (!is.na(data1$KADEATH_DT) & data1$INTDT04<=data1$KADEATH_DT),]
data1$timeindys<-data1$INTDT04
data1$ELISA_m<-data1$EIACU04
data1$LST_m<-data1$MEANLST04
data1$age_m<-data1$AGE04
data1$KA_m<-data1$KA_04
data1$death_m<-0
data1$timepoint<-"S3"
data_long<-rbind(data_long,data1)
dim(data_long)

######## CLEAN ########
## Remove survey points when no LST and no ELISA - don't move this or will lose some of the KA observations
no_LST_ELISA_idx<-is.na(data_long$LST_m) & is.na(data_long$ELISA_m)
data_long<-data_long[!no_LST_ELISA_idx,]
dim(data_long)

######## ADD OBSERVATIONS OF KA, TREATMENT AND DEATH ########
# KA observations
data1<-data[data$KA==1,]
dim(data1)
data1$timeindys<-data1$FEV_ONS
data1$ELISA_m<-data1$EIA_KA
data1$LST_m<-data1$LST_KA
data1$age_m<-(age[data$KA==1]+data1$timeindys/365)
data1$KA_m<-1
data1$death_m<-0
data1$timepoint[data1$FEV_ONS<data1$INTDT02]<-"K01"
data1$timepoint[is.na(data1$FEV_ONS) & !is.na(data1$ONSYR) & data1$ONSYR<=2002]<-"K01"
data1$timepoint[!is.na(data1$FEV_ONS) & data1$FEV_ONS==data1$INTDT02]<-"KS1"
data1$timepoint[!is.na(data1$FEV_ONS) & data1$FEV_ONS>data1$INTDT02 & data1$FEV_ONS<data1$INTDT03]<-"K12"
data1$timepoint[!is.na(data1$FEV_ONS) & data1$FEV_ONS==data1$INTDT03]<-"KS2"
data1$timepoint[!is.na(data1$FEV_ONS) & data1$FEV_ONS>data1$INTDT03 & data1$FEV_ONS<data1$INTDT04]<-"K23"
data1$timepoint[!is.na(data1$FEV_ONS) & data1$FEV_ONS==data1$INTDT04]<-"KS3"
data1$timepoint[!is.na(data1$FEV_ONS) & data1$FEV_ONS>data1$INTDT04]<-"K30"
# Assume any remaining unlabelled observations of KA onset for individuals without a 3rd survey date occurred between 2nd and 3rd surveys
data1$timepoint[is.na(data1$timepoint) & is.na(data1$INTDT04)]<-"K23"
# UNCOMMENTING THIS LINE SETS SYMPTOM ONSET DATE TO START OF ONSET YEAR FOR KA PATIENTS WITHOUT A SYMPTOM ONSET DATE
# data1$timeindys[data1$timepoint=="K01"]<-(data1$ONSYR[data1$timepoint=="K01"]-as.numeric(format(as.Date(STARTDT),"%Y")))*365

data_long<-rbind(data_long,data1)
dim(data_long)

# Treatment - observation at the end of the treatment period (RXDT+trtmnt_prd)
data1<-data[data$KA==1 & !is.na(data$RXDT) & (is.na(data$KADEATH_DT) | (data$RXDT+trtmnt_prd)!=data$KADEATH_DT),]
dim(data1)
data1$timeindys<-(data1$RXDT+trtmnt_prd)
data1$ELISA_m<-NA
data1$LST_m<-NA
data1$age_m<-(age[data$KA==1 & !is.na(data$RXDT) & (is.na(data$KADEATH_DT) | (data$RXDT+trtmnt_prd)!=data$KADEATH_DT)]+data1$timeindys/365)
data1$death_m<-0  
data1$KA_m<-0
data1$timepoint<-"TR"
dim(data1)
data_long<-rbind(data_long,data1)
dim(data_long)

# Relapse to KA
data1<-data[!is.na(data$RELAPSE),]
dim(data1)
data1$timeindys<-data1$RELAPSE
data1$ELISA_m<-NA
data1$LST_m<-NA
data1$age_m<-(age[!is.na(data$RELAPSE)]+data1$timeindys/365)
data1$death_m<-0  
data1$KA_m<-1
data1$timepoint<-"R"
dim(data1)
data_long<-rbind(data_long,data1)
dim(data_long)

# Re-treatment
data1<-data[!is.na(data$RELP_RX),]
dim(data1)
data1$timeindys<-data1$RELP_RX
data1$ELISA_m<-NA
data1$LST_m<-NA
data1$age_m<-(age[!is.na(data$RELP_RX)]+data1$timeindys/365)
data1$death_m<-0  
data1$KA_m<-0
data1$timepoint<-"TR_R"
dim(data1)
data_long<-rbind(data_long,data1)
dim(data_long)

# KA deaths
data1<-data[!is.na(data$KADEATH_DT),]
data1$timeindys<-data1$KADEATH_DT
data1$ELISA_m<-NA
data1$LST_m<-NA
data1$KA_m<-0
data1$death_m<-1
data1$age_m<-data1$AGE_KADEATH
data1$timepoint<-"D_KA"
data_long<-rbind(data_long,data1)
dim(data_long)

# Natural deaths
data1<-data[!is.na(data$NATDEATH_DT),]
data1$timeindys<-data1$NATDEATH_DT
data1$ELISA_m<-NA
data1$LST_m<-NA
data1$KA_m<-0
data1$death_m<-1
data1$age_m<-age[!is.na(data$NATDEATH_DT)]+data1$timeindys/365
data1$timepoint<-"D_NAT"
data_long<-rbind(data_long,data1)
dim(data_long)

######## CLEAN ########
## Remove any observations without a time
data_long<-data_long[!is.na(data_long$timeindys),]
dim(data_long)

## Remove data points after death
predeath_idx<-((is.na(data_long$KADEATH_DT) & is.na(data_long$NATDEATH_DT)) | (!is.na(data_long$KADEATH_DT) & data_long$timeindys<=data_long$KADEATH_DT) | (!is.na(data_long$NATDEATH_DT) & data_long$timeindys<=data_long$NATDEATH_DT))
data_long<-data_long[predeath_idx,]

## Re-order data by the patient ID, IDNUM, and time of observation
data_long<-data_long[order(data_long$IDNUM,data_long$timeindys),]

######## CLASSIFY DATA (ASSIGN STATES TO EACH OBSERVATION) ########
## Make extra markers
# for having ELISA reading
hasE<-!is.na(data_long$ELISA_m)
data_long$hasE<-hasE
# for being ELISA positive
posE<-(!is.na(data_long$ELISA_m) & data_long$ELISA_m>=ELISA_thr)
data_long$posE<-posE
# for being ELISA negative
negE<-(!is.na(data_long$ELISA_m) & data_long$ELISA_m<ELISA_thr)
data_long$negE<-negE
# for having LST reading
hasL<-!is.na(data_long$LST_m)
data_long$hasL<-hasL
# for being LST positive
posL<-(!is.na(data_long$LST_m) & data_long$LST_m>=LST_thr)
data_long$posL<-posL
# for being LST negative
negL<-(!is.na(data_long$LST_m) & data_long$LST_m<LST_thr)
data_long$negL<-negL
# for active KA
active_KA<-(data_long$KA_m==1)
data_long$active_KA<-active_KA
# for being pre_KA (i.e. having an onset date in the future)
pre_KA<-logical(nrow(data_long))
pre_KA[!is.na(data_long$FEV_ONS) & data_long$timeindys<data_long$FEV_ONS]<-TRUE
data_long$pre_KA<-pre_KA
# for being post KA (i.e. having been treated for KA or more than a year after symptom onset)
post_KA<-logical(nrow(data_long))
# for end of treatment
end_of_trtmnt<-(data_long$timepoint=="TR" | data_long$timepoint=="TR_R")
# post treatment
post_KA[data_long$KA==1 & !is.na(data_long$RXDT) & is.na(data_long$KADEATH_DT) & data_long$timeindys>=data_long$RXDT+trtmnt_prd]<-TRUE 
# no treatment date, but more than a year after onset
post_KA[data_long$KA==1 & is.na(data_long$RXDT) & is.na(data_long$KADEATH_DT) & ((!is.na(data_long$FEV_ONS) & data_long$timeindys>=data_long$FEV_ONS+symptm_prd) | (is.na(data_long$FEV_ONS) & data_long$timeindys>(data_long$ONSYR-as.numeric(format(as.Date(STARTDT),"%Y")))*365+symptm_prd))]<-TRUE
data_long$post_KA<-post_KA

## Make a function to test whether the individual was LST positive at any time in the past
post_LST_pstve<-function(x)
{
  post_LST_pstve_idx<-numeric(nrow(x))
  for (i in 1:nrow(x))
  {
    LST_results<-x[i,c("MEANLST02", "MEANLST03", "MEANLST04")]
    test_times<-x[i,c("INTDT02", "INTDT03", "INTDT04")]
    crrnt_time<-x$timeindys
    post_LST_pstve_idx[i]<-any(LST_results[test_times<=crrnt_time]>=LST_thr,na.rm=TRUE) # has the individual ever had positive LST?
  }
  return(as.logical(post_LST_pstve_idx))
}

## Make a function to test whether the individual had positive ELISA on the previous test
last_ELISA_pstve<-function(x)
{
  last_ELISA_pstve_idx<-numeric(nrow(x))
  for (i in 1:nrow(x)) #
  {
    ELISA<-x[i,c("EIACU02", "EIACU03", "EIACU04")]
    test_times<-x[i,c("INTDT02", "INTDT03", "INTDT04")]
    
    crrnt_time<-x$timeindys[i]
    ttimes<-test_times[test_times<crrnt_time]
    
    if ((length(ttimes)>0)&(sum(is.na(ttimes))!=length(ttimes)))    
    {
      previous_time<-which.max(ttimes)
      
      last_ELISA_pstve_idx[i]<-any(ELISA[previous_time]>=ELISA_thr) # did the individual have positive ELISA on the previous test?
    }
    else
    {     
      last_ELISA_pstve_idx[i]<-FALSE # no previous test
    }
    
  }
  return(as.logical(last_ELISA_pstve_idx))
}

## Apply last_ELISA_pstve and post_LST_pstve functions to dataframe
last_ELISA_pstve_idx<-last_ELISA_pstve(data_long) # ELISA+ on previous test
data_long$last_ELISA_pstve<-last_ELISA_pstve_idx
post_LST_pstve_idx<-post_LST_pstve(data_long) # positive LST in the past

# print(table(post_LST_pstve_idx,pre_KA))

## Convert units for time from days to years
data_long$time<-data_long$timeindys/365

## Classify ages into groups
data_long$age_groups[data_long$age_m>=0 & data_long$age_m<15]<-"A"
data_long$age_groups[data_long$age_m>=15 & data_long$age_m<=45]<-"B"
data_long$age_groups[data_long$age_m>45]<-"C"

## Classify state for each observation
# Susceptible: KA-, LST never positive, no previous KA, ELISA- and ELISA- on previous survey
state1_idx<-(!active_KA & negE & negL & !last_ELISA_pstve_idx & !end_of_trtmnt)

# Asymptomatic: KA-, ELISA+ and LST-
state2_idx<-(!active_KA & posE & negL & !end_of_trtmnt & !post_KA)

# Active KA: KA+ & not finished treatment
state3_idx<-(active_KA & !end_of_trtmnt)

# Recovered/dormant: LST+ or (ELISA+ & post KA) or (missing ELISA & LST- & post KA) or serodeconverted or at end of treatment
state4_idx<-(!active_KA & (posL
                           | (posE & post_KA)
                           | (!hasE & negL & post_KA & last_ELISA_pstve_idx)
                           | (last_ELISA_pstve_idx & negE)
                           | end_of_trtmnt)) 

# Dead
state5_idx<-(data_long$death_m==1)

## Censored observations: for missing LST or ELISA or previous survey ELISA readings
# State 1 or 4: (E- & missing LST & prev ELISA -/missing) or (missing ELISA & L- & post KA & prev ELISA -/missing) or (E- & L- & missing prev ELISA)
state91_idx<-(!active_KA & !end_of_trtmnt & ((negE & !hasL & (!last_ELISA_pstve_idx | is.na(last_ELISA_pstve_idx)))
                          | (!hasE & negL & post_KA & (!last_ELISA_pstve_idx | is.na(last_ELISA_pstve_idx)))
                          | (negE & negL & is.na(last_ELISA_pstve_idx))))
# State 2 or 4: for pre KA or no KA: either (E+ & missing LST) or (missing ELISA & L- & possibility of serodeconversion)
state92_idx<-(!active_KA & !end_of_trtmnt & ((posE & !hasL & !post_KA) 
                          | (!hasE & negL & last_ELISA_pstve_idx & !post_KA )))
# State 1 or 2: for preK or no K: missing ELISA, L- & no possibility of serodeconversion
state93_idx<-(!active_KA & !end_of_trtmnt & !hasE & negL & !last_ELISA_pstve_idx & !post_KA)

# Make vector of states for observations
state<-numeric(nrow(data_long))
state[state1_idx]<-1
state[state2_idx]<-2
state[state3_idx]<-3
state[state4_idx]<-4
state[state5_idx]<-5
state[state91_idx]<-91
state[state92_idx]<-92
state[state93_idx]<-93

# print(summary(as.factor(state)))

## Classify type of each observation (see "obstype" on page 38 of 2015 msm manual):
# 1 = observation at arbitrary time
# 2 = exact transition time with state at previous observation retained until current observation (can be 
# transition to different state or repeated observation of same state)
# 3 = exact transition time but previous state is unknown (e.g. KA onset without previous observation)
obstype<-numeric(nrow(data_long))
# Classify observations of states 1, 2 and 4 (apart from at treatment end date) and censored observations as at arbitrary times
obstype[state1_idx|state2_idx|(state4_idx & !end_of_trtmnt)|state91_idx|state92_idx|state93_idx]<-1
# Classify observations of active KA for which previous state is unknown as type 3 exact (see below)
obstype[state3_idx & !is.na(data_long$FEV_ONS)]<-3
# obstype[state3_idx & is.na(data_long$FEV_ONS)]<-1
# Classify observations of state 4 at end of treatment as type 2 exact
obstype[state4_idx & end_of_trtmnt]<-2
# Classify observations of dead KA patients as type 2 exact if death date was not set to end of year and observations of natural deaths as type 2 exact if month of death was recorded
obstype[state5_idx & (format(as.Date(data_long$KADEATH_DT,origin=STARTDT),"%m-%d")!="12-31" | !is.na(data_long$NATDEATH_MO))]<-2
# Classify observations of dead KA patients as at arbitrary times if death date was set to end of year and observations of natural deaths as at arbitrary times if month of death was not recorded
obstype[state5_idx & (format(as.Date(data_long$KADEATH_DT,origin=STARTDT),"%m-%d")=="12-31" | is.na(data_long$NATDEATH_MO))]<-1

# Overwrite observation type for exactly observed transitions to KA where previous state is known (i.e. is either state 2 or state 3)
prev_obsv_IDNUM<-c(0,data_long$IDNUM[1:nrow(data_long)-1])
prev_state2_idx<-c(logical(1),state2_idx[1:length(state2_idx)-1])
prev_state3_idx<-c(logical(1),state3_idx[1:length(state3_idx)-1])
obstype[data_long$IDNUM==prev_obsv_IDNUM & state3_idx & (prev_state2_idx | prev_state3_idx)]<-2
# for relapse
prev_state4_idx<-c(logical(1),state4_idx[1:length(state4_idx)-1])
obstype[data_long$IDNUM==prev_obsv_IDNUM & state3_idx & prev_state4_idx & !is.na(data_long$RELAPSE)]<-2

# Append vectors of states and observation types to dataframe
data_long<-cbind(data_long,state,obstype)

# Uncomment to write output to data_long_all_deaths_and_relapse.csv
# write.csv(data_long,"data_long_all_deaths_and_relapse.csv")

