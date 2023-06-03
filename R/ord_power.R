#‘ calculate the mean of estimated cross-lag coefficient and percentage of significant effects from number of replicated simulations
#' This function allows you to simulate longitudinal ordinal data and estimate the parameters with multilevel models for a nubmer of times
#' @param n number of categories for the ordinal outcome variable
#' @param numSample number of participants
#' @param numAssess number of assessments
#' @param thresh_CON threshold conditions, 1 represents normal thresholds, 2 represents right-skewed thresholds
#' @param autoreg_coeff autoregressive coefficient
#' @param crosslag_coeff cross-lag coefficient
#' @param crosslag_skew cross-lag skewness
#' @param crosslag_kurt cross-lag kurtosis
#' @param gamma_00 fixed intercept
#' @param gamma_00_sd random intercept
#' @param gamma_01_sd random autoregressive cofficient sd
#' @param gamma_02_sd random cross-lag coefficient sd
#' @param Compliance compliance rate in percentage
#' @param ar_skew ar skewness
#' @param ar_kurt ar kurtosis
#' @param corr correlation between random ar and random crosslag, 1 indicating correlation of 0.5, 2 indicating no corrrelation
#' @param reps number of replications
#' @return the mean of estimated cross-lag coefficient, the percentage of significant effects, the percentage of non-significant effects, total replication number, percentage of replications that produced errors, and percentage of replications that produced warnings
#' @import EMAtools
#' @import Rcpp
#' @import Matrix
#' @import ordinal
#' @importFrom stats rnorm runif
#' @export
#' @examples
#' \donttest{
#' ord_power(n = 5,
#' numSample = 35,
#' numAssess = 73,
#' thresh_CON = 2,
#' autoreg_coeff = 0.919,
#' crosslag_coeff = 0.4,
#' gamma_00 = 0,
#' gamma_00_sd = 1.043,
#' gamma_01_sd = 0.591,
#' gamma_02_sd = 0.13,
#' Compliance = 1,
#' reps=10)
#' }
ord_power <- function(n,
                      numSample,
                      numAssess,
                      thresh_CON,
                      autoreg_coeff,
                      crosslag_coeff,
                      crosslag_skew = 0,
                      crosslag_kurt = 3,
                      gamma_00,
                      gamma_00_sd,
                      gamma_01_sd,
                      gamma_02_sd,
                      Compliance = 1,
                      ar_skew = 0,
                      ar_kurt = 3,
                      corr = 1,
                      reps = 100){

  # cl <- parallel::makePSOCKcluster(1)
  #doParallel::registerDoParallel(cl)
  #X <- 1:reps

  #Y=foreach::foreach(x = X, .packages=c('tidyverse','DataCombine','EMAtools',
  #                             'ordinal', 'covsim', 'brms', 'Matrix','rstan','dbplyr','Rcpp'), .export = c("get_param","extract_modparams")) %dopar% {


  # n = 7
  # numSample = 80
  # numAssess = 28
  # thresh_CON = 1
  # autoreg_coeff = 0.2
  # crosslag_coeff = 0.2
  # crosslag_skew = 0
  # crosslag_kurt = 3
  # gamma_00 = 1
  # gamma_00_sd = 0.2
  # gamma_01_sd = 0
  # gamma_02_sd = 0
  # Compliance = 1
  # ar_skew = 0
  # ar_kurt = 3
  # corr = 1
  # reps=20


  obs <- replicate(reps, try(get_param(n,
                                     numSample,
                                     numAssess,
                                     thresh_CON,
                                     autoreg_coeff,
                                     crosslag_coeff,
                                     crosslag_skew = crosslag_skew,
                                     crosslag_kurt = crosslag_kurt,
                                     gamma_00,
                                     gamma_00_sd,
                                     gamma_01_sd,
                                     gamma_02_sd,
                                     Compliance,
                                     ar_skew = ar_skew,
                                     ar_kurt = ar_kurt,
                                     corr = corr)))




  # out <- do.call(rbind, obs)
  out <- obs


  crosslag_est <- out[1,]
  crosslag_p <- out[2,]
  error<- out[3,]

  # mean(crosslag_est[crosslag_est< quantile(crosslag_est, 0.9, na.rm=T)], na.rm=T)

  crosslag = round(mean(crosslag_est, na.rm=T),2)
  count1 <- sum(crosslag_p < 0.05, na.rm = T) # detect as significant
  count2 <- sum(crosslag_p >= 0.05, na.rm = T) # detect as nonsignificant


  perc_sig <- round(count1/ncol(out),2)
  perc_notsig <- round(count2/ncol(out),2)
  errorPerc <- round(sum(error == 2)/ncol(out),2)
  warningPerc <- round(sum(error == 1)/ncol(out),2)
  res <- list(crosslag,
              perc_sig,
              perc_notsig,
              nrow(out),
              errorPerc,
              warningPerc)

  return(res)
}



