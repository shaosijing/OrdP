#‘ calculate the mean of estimated cross-lag coefficient and percentage of significant effects from number of replicated simulations
#' This function allows you to simulate longitudinal ordinal data and estimate the parameters with multilevel models for a nubmer of times
#' @param n number of categories for the ordinal outcome variable
#' @param numSample number of participants
#' @param numAssess number of assessments per participant
#' @param thresh_CON condition for thresholds: normal or skewed
#' @param autoreg_coeff autoregressive coefficient
#' @param crosslag_coeff cross-lag coefficient
#' @param crosslag_sk cross-lag skewness and kuertosis
#' @param gamma_00 fixed intercept
#' @param gamma_00_sd random intercept
#' @param gamma_01_sd random autoregressive cofficient
#' @param gamma_02_sd random cross-lag coefficient sd
#' @param Compliance compliance rate in percentage
#' @param crosslag_prior prior for crosslag
#' @param ar_sk ar skewness and kuertosis
#' @param corr correlation between random ar and random crosslag
#' @param reps number of replications
#' @return the mean of estimated cross-lag coefficient and percentage of significant effects from number of replicated simulations
#' @import EMAtools
#' @import Rcpp
#' @import Matrix
#' @import brms
#' @import tidyverse
#' @import dbplyr
#' @import ordinal
#' @import rstan
#' @import Brobdingnag
#' @importFrom stats rnorm runif
#' @export



ord_power<-function(n,numSample,numAssess,thresh_CON,autoreg_coeff,crosslag_coeff,crosslag_sk,gamma_00,gamma_00_sd, gamma_01_sd,gamma_02_sd,Compliance,crosslag_prior,ar_sk, corr, reps){

   cl <- parallel::makePSOCKcluster(2)
  doParallel::registerDoParallel(cl)
  X <- 1:reps

  Y=foreach::foreach(x = X, .packages=c('tidyverse','DataCombine','EMAtools',
                               'ordinal', 'covsim', 'brms', 'Matrix','rstan','dbplyr','Rcpp'), .export = c("get_param","extract_modparams")) %dopar% {
                                 obs<-get_param(n,numSample,numAssess,thresh_CON,autoreg_coeff,crosslag_coeff,crosslag_sk,gamma_00,gamma_00_sd, gamma_01_sd,gamma_02_sd,Compliance,crosslag_prior,ar_sk, corr)
                               }

  parallel::stopCluster(cl)

  out<-extract_modparams(Y, reps)
  crosslag_est<-out[,1]
  crosslag_p<-out[,2]

  crosslag = round(mean(crosslag_est),2)
  count1<- sum(crosslag_p<0.05)
  count2<- sum(crosslag_p>=0.05)

  if (nrow(out) > 0){
  perc_sig<-round(count1/nrow(out),2)
  perc_notsig<-round(count2/nrow(out),2)
  } else if (nrow(out) == 0){
  perc_sig<-NA
  perc_notsig<-NA
  }
  res<-list(crosslag,perc_sig, perc_notsig, nrow(out))
  return(res)
}
