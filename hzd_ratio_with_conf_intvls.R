hzd_ratio_with_conf_intvls<-function(x)
{
HR<-(x$obs[2]/x$exp[2])/(x$obs[1]/x$exp[1])
# UL95CI<-exp(log(HR)+qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
# LL95CI<-exp(log(HR)-qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
UL95CI<-exp((x$obs[2]-x$exp[2])/x$var[2,2]+qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
LL95CI<-exp((x$obs[2]-x$exp[2])/x$var[2,2]-qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
return(c(HR,LL95CI,UL95CI))
}