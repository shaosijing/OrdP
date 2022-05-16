#‘ calculate the mean of estimated cross-lag coefficient and percentage of significant effects from number of replicated simulations
#' This function allows you to simulate longitudinal ordinal data and estimate the parameters with multilevel models for a nubmer of times
#' @param n number of categories for the ordinal outcome variable
#' @param numSample number of participants
#' @param numAssess number of assessments per participant
#' @param thresh thresholds for ordinal outcome
#' @param autoreg_coeff autoregressive coefficient
#' @param crosslag_coeff cross-lag coefficient
#' @param gamma_00 fixed intercept
#' @param gamma_00_sd random intercept
#' @param gamma_01_sd random autoregressive cofficient
#' @param gamma_02_sd random cross-lag coefficient sd
#' @param Compliance compliance rate in percentage
#' @param reps number of replications
#' @return the mean of estimated cross-lag coefficient and percentage of significant effects from number of replicated simulations
#' @export



ord_power<-function(n,numSample,numAssess,thresh,autoreg_coeff,crosslag_coeff,gamma_00,gamma_00_sd,gamma_01_sd,gamma_02_sd,Compliance, reps){
  cl <- parallel::makePSOCKcluster(2)
  doParallel::registerDoParallel(cl)
  X <- 1:reps

  Y=foreach::foreach(x = X, .packages=c('mnormt','tidyverse','DataCombine','EMAtools',
                               'ordinal'), .export = c("get_param","extract_modparams")) %dopar% {
                                 obs<-get_param(n,numSample,numAssess,thresh,autoreg_coeff,crosslag_coeff,gamma_00,gamma_00_sd,gamma_01_sd, gamma_02_sd, Compliance)
                               }

  parallel::stopCluster(cl)

  out<-extract_modparams(Y, reps)
  crosslag_est<-out[,1]
  crosslag_p<-out[,2]

  crosslag = round(mean(crosslag_est),2)
  count<- sum(crosslag_p<0.05)
  perc<-round(count/reps,2)
  res<-list(crosslag,perc)
  return(res)
}
