#‘ simulate longitudinal ordinal data and estimate the parameters with multilevel models
#' This function allows you to simulate longitudinal ordinal data and estimate the parameters with multilevel models
#' @param n number of categories for the ordinal outcome variable
#' @param numSample number of participants
#' @param numDays number of days
#' @param assess four accessments per day
#' @param thresh thresholds for ordinal outcome
#' @param autoreg_coeff autoregressive coefficient
#' @param crosslag_coeff cross-lag coefficient
#' @param gamma_00 fixed intercept
#' @param gamma_00_sd random intercept
#' @param gamma_01_sd random autoregressive cofficient
#' @return the estimated cross-lag coefficient and its corresponding p-value
#' @export

get_para<-function(n,numSample,numDays,assess,thresh,autoreg_coeff,crosslag_coeff,gamma_00,gamma_00_sd, gamma_01_sd){
  # numSample: number of participants
  # numDays: number of EMA days
  # assess: number of measurements each day (same for everyone)
  # thresh: thresholds, fixed for everyone
  # autoreg: autoregressive coefficient
  # crosslag: crosslag coefficient
  # gamma_00: fixed intercept
  # gamma_00_sd: random effect for intercept


  N = 1:numSample
  days = 1:numDays

  datt = data.frame(expand.grid(assess,days,N))
  colnames(datt) = c("hour","day","N")
  datt$timepoint = rep(1:(max(days)*length(assess)),length(N))
  datt$si_cat = NA
  datt$si_star = NA
  datt$pred = rnorm(nrow(datt))

  # autoreg <- rep(autoreg_coeff,max(N)) #leave out random effect for now
  #int<-  rnorm(numSample,gamma_00,gamma_00_sd)
  #autoreg <- rnorm(max(N),autoreg_coeff, 0.569)
  crosslag <- rep(crosslag_coeff,max(N))

  # random intercept and random slope in autoregressive
  gam <- c(gamma_00, autoreg_coeff)
  G<-matrix(c(gamma_00_sd,-0.54,-0.54,gamma_01_sd),nrow = 2)
  uj <- mnormt::rmnorm(max(N), mean = rep(0, 2), varcov = G)
  betaj <- matrix(gam, nrow = max(N), ncol = 2, byrow = TRUE) + uj
  int<-betaj[,1]
  autoreg<- betaj[,2]

  count = 0

  thresh[5]=100
  thresh1=append(thresh,-100,0)
  datt$si_cat = NULL

  for(i in 1:nrow(datt)){

    if(datt[i,"timepoint"]==1){
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
                 NewVar="si_cat_lead",slideBy=1,TimeVar="timepoint")



  # brm_int_sub1 <- brm(si_cat_lead ~ 1+si_cat+pred+(1|N),
  #                      family = cumulative("probit"), data = datt2,
  #                     chains = 3, iter = 400, cores = 1)

  # sum=summary(brm_int_sub1)
  # auto_est_CI_brms<- c(sum$fixed$`l-95% CI`[5],sum$fixed$`u-95% CI`[5])

  # clmm
  datt2$si_cat_lead<-as.factor(datt2$si_cat_lead)
  #datt2$si_cat<-as.factor(datt2$si_cat)
  mod=ordinal::clmm2(si_cat_lead ~ si_cat+pred+(1|N), data = datt2, link = "probit")
  sum=summary(mod)
  res<-list(c(sum$coefficients[6,1],sum$coefficients[6,4]))
  #coefficients[6,1]: est
  #coefficients[6,4]: p value
  return(res)

}
