#‘ simulate longitudinal ordinal data and estimate the parameters with multilevel models
#' This function allows you to simulate longitudinal ordinal data and estimate the parameters with multilevel models
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
#' @return the estimated cross-lag coefficient, its corresponding p-value and whether an error (2) or warning (1) message is produced
#' @export

get_param <- function(n,
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
                    Compliance,
                    ar_skew = 0,
                    ar_kurt = 3,
                    corr = 1){

  N = 1:numSample
  assess = 1:numAssess

  datt = data.frame(expand.grid(assess,N))
  colnames(datt) = c("assessment","N")
  datt$si_cat = NA
  datt$si_star = NA
  datt$pred = rnorm(nrow(datt))

  if (thresh_CON == 1){ #normal
    thresh = c(1:(n-1))
    thresh = thresh - mean(thresh)
  } else if (thresh_CON == 2){ #skewed
    thresh = c(1:(n-1))
  }

  # random intercept and random slope in autoregressive
  #gam <- c(gamma_00, autoreg_coeff,crosslag_coeff)
  if (gamma_01_sd == 0 & gamma_02_sd == 0){
    int<-gamma_00 + rnorm(max(N), 0, gamma_00_sd)
    autoreg<- rep(autoreg_coeff, max(N))
    crosslag <- rep(crosslag_coeff,max(N))
  } else if (gamma_01_sd != 0 & gamma_02_sd == 0){ # no variablity in cross-lag coefficient
    if (corr ==1 ){
      a<-matrix(c(1,0.5,0.5, 1),nrow = 2)
    } else if (corr ==2){
      a<-matrix(c(1,0,0,1),nrow = 2)
    }
    stdevs <- c(gamma_00_sd,gamma_01_sd)
    b <- stdevs %*% t(stdevs)
    G <- b * a

    gam <- c(gamma_00, autoreg_coeff)

    uj <- covsim::rIG(max(N), sigma = G, skewness = c(0,ar_skew),
                      excesskurt = c(3,ar_kurt), reps = 1)[[1]]

    int<-uj[,1]
    autoreg<- uj[,2]
    crosslag <- rep(crosslag_coeff,max(N))
  } else if (gamma_01_sd == 0 & gamma_02_sd != 0){
    if (corr ==1 ){
      a<-matrix(c(1,0.5,0.5, 1),nrow = 2)
    } else if (corr ==2){
      a<-matrix(c(1,0,0,1),nrow = 2)
    }
    stdevs <- c(gamma_00_sd,gamma_02_sd)
    b <- stdevs %*% t(stdevs)
    G <- b * a

    gam <- c(gamma_00, crosslag_coeff)

    uj <- covsim::rIG(max(N), sigma = G, skewness = c(0,crosslag_skew),
                      excesskurt = c(3,crosslag_kurt), reps = 1)[[1]]
    int<-uj[,1]
    autoreg<- rep(autoreg_coeff,max(N))
    crosslag <- uj[,2]
  } else if (gamma_01_sd != 0 & gamma_02_sd != 0){
    gam <- c(gamma_00, autoreg_coeff, crosslag_coeff)

    if (corr ==1 ){
      a<-matrix(c(1,-0.5,0.5,-0.5, 1, -0.5, 0.5,-0.5, 1),nrow = 3)
    } else if (corr ==2){
      a<-matrix(c(1,0,0,0, 1, 0, 0, 0, 1),nrow = 3)
    }
    stdevs <- c(gamma_00_sd,gamma_01_sd,gamma_02_sd)
    b <- stdevs %*% t(stdevs)
    G <- b * a

    uj <- covsim::rIG(max(N), sigma = G, skewness = c(0,ar_skew, crosslag_skew),
                      excesskurt = c(3,ar_kurt, crosslag_kurt), reps = 1)[[1]]

    int<-uj[,1]
    autoreg<- uj[,2] +autoreg_coeff
    crosslag<-uj[,3] + crosslag_coeff
  }

  count = 0

  thresh[n]=100
  thresh1=append(thresh,-100,0)
  datt$si_cat = NULL

  for(i in 1:nrow(datt)){

    if(datt[i,"assessment"]==1){
      count = count + 1
      datt[i,"si_star"] = int[count] + rnorm(1,0,1)
    }else{
      datt[i,"si_star"] = int[count] + autoreg[count]*datt[i-1,"si_cat"] + crosslag[count]*datt[i-1,"pred"] + rnorm(1,0,1)
    }

    for (j in (1:n)){
      if (datt[i,"si_star"] >= thresh1[j] & datt[i,"si_star"] < thresh1[j+1] || datt[i,"si_star"] >= thresh1[j+1] ){
        datt[i,"si_cat"] = j
      }
    }


  }

  datt2 <- DataCombine::slide(datt,Var="si_cat",GroupVar="N",
                              NewVar="si_cat_lead",slideBy=1,TimeVar="assessment")


  prop.m = 1-Compliance/100  # 7% missingness
  mcar   = runif(nrow(datt2), min=0, max=1)
  datt2$si_cat_lead = ifelse(mcar<prop.m, NA, datt2$si_cat_lead)

  datt2$si_cat_lead<-as.factor(datt2$si_cat_lead)
  #datt2$si_cat<-as.factor(datt2$si_cat)

  thresh_length = length(table(datt2$si_cat_lead))
  datt2$N <- factor(datt2$N)
  #randnum<-floor(runif(1,0,1)*1000)




  if (gamma_01_sd == 0){
    factory1 <- function() {
      warn <- err <- NULL
      res <- withCallingHandlers(
        tryCatch(
          # fun(x),
          ordinal::clmm(si_cat_lead ~ si_cat + pred + (1|N), data=datt2, link = "probit"),
          error = function(e) {
            err <<- conditionMessage(e)
            NULL
          }),
        warning = function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      list(res=res, warn = warn, err = err)
    }
    mod1 <- factory1()
  } else if (gamma_02_sd == 0){
    factory2 <- function() {
      warn <- err <- NULL
      res <- withCallingHandlers(
        tryCatch(
          # fun(x),
          ordinal::clmm(si_cat_lead ~ si_cat + pred + (1+si_cat|N), data=datt2,link = "probit"),
          error = function(e) {
            err <<- conditionMessage(e)
            NULL
          }),
        warning = function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      list(res=res, warn = warn, err = err)
    }


    mod2<-factory2()
  }else{
    factory3 <- function() {
      warn <- err <- NULL
      res <- withCallingHandlers(
        tryCatch(
          # fun(x),
          ordinal::clmm(si_cat_lead ~ si_cat+pred+(1+si_cat+pred|N), data = datt2, link = "probit"),
          error = function(e) {
            err <<- conditionMessage(e)
            NULL
          }),
        warning = function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      list(res=res, warn = warn, err = err)
    }

    mod3<-factory3()
  }







  if (gamma_01_sd == 0){
    sum.m <- summary(mod1$res)
    mod.m <- mod1
  } else if (gamma_02_sd == 0){
    sum.m <- summary(mod2$res)
    mod.m <- mod2
  }else{
    sum.m <- summary(mod3$res)
    mod.m <- mod3
  }



  if (!is.null(mod.m$err)){
    res <- c(NA, NA, 2)
  }else if (!is.null(mod.m$warn)){
    res <- c(sum.m$coefficients[thresh_length+1,1],sum.m$coefficients[thresh_length+1,4], 1)
  }else{
    res <- c(sum.m$coefficients[thresh_length+1,1],sum.m$coefficients[thresh_length+1,4], 0)
  }


  #coefficients[6,1]: est
  #coefficients[6,4]: p value
  return(res)

}


# get_param(n = 7,
#                     numSample = 80,
#                     numAssess = 28,
#                     thresh_CON = 1,
#                     autoreg_coeff = 0.2,
#                     crosslag_coeff = 0.2,
#                     crosslag_skew = 0,
#                     crosslag_kurt = 3,
#                     gamma_00 = 1,
#                     gamma_00_sd = 0.2,
#                     gamma_01_sd = 0,
#                     gamma_02_sd = 0,
#                     Compliance = 1,
#                     ar_skew = 0,
#                     ar_kurt = 3,
#                     corr = 1)

