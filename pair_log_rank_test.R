pair_log_rank_test<-function(data,group,gp1,gp2)
{
  gp12idx<-(group==gp1 | group==gp2)
  survdiff.obj<-survdiff(Surv(data$time[gp12idx],data$KA[gp12idx])~group[gp12idx],rho=0)
  return(survdiff.obj)
}