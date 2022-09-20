#‘ simulate longitudinal ordinal data and estimate the parameters with multilevel models
#' This function allows you to simulate longitudinal ordinal data and estimate the parameters with multilevel models
#' @param n number of categories for the ordinal outcome variable
#' @param numSample number of participants
#' @param numAssess number of assessments
#' @param thresh_CON condition for thresholds: normal or skewed
#' @param autoreg_coeff autoregressive coefficient
#' @param crosslag_coeff cross-lag coefficient
#' @param crosslag_sk cross-lag skewness and kuertosis
#' @param gamma_00 fixed intercept
#' @param gamma_00_sd random intercept
#' @param gamma_01_sd random autoregressive cofficient sd
#' @param gamma_02_sd random cross-lag coefficient sd
#' @param Compliance compliance rate in percentage
#' @param crosslag_prior prior for crosslag
#' @param ar_sk ar skewness and kuertosis
#' @param corr correlation between random ar and random crosslag
#' @return the estimated cross-lag coefficient and its corresponding p-value
#' @export

get_param<-function(n,numSample,numAssess,thresh_CON,autoreg_coeff,crosslag_coeff,crosslag_sk,gamma_00,gamma_00_sd, gamma_01_sd,gamma_02_sd,Compliance,crosslag_prior,ar_sk, corr){

  N = 1:numSample
  assess = 1:numAssess

  datt = data.frame(expand.grid(assess,N))
  colnames(datt) = c("assessment","N")
  datt$si_cat = NA
  datt$si_star = NA
  datt$pred = rnorm(nrow(datt))

  if (crosslag_sk == 1){
    crosslag_skew = 0
    crosslag_kurt = 3
  } else if (crosslag_sk == 2){
    crosslag_skew = 0.7
    crosslag_kurt = 3
  } else if (crosslag_sk == 3){
    crosslag_skew = 1.5
    crosslag_kurt = 3
  }

  if (ar_sk == 1){
    ar_skew = 0
    ar_kurt = 3
  } else if (ar_sk == 2){
    ar_skew = 0.7
    ar_kurt = 3
  } else if (ar_sk == 3){
    ar_skew = 1.5
    ar_kurt = 3
  }

  if (thresh_CON == 1){ #normal
    thresh = c(-1.5, -1, 1, 1.5)
  } else if (thresh_CON == 2){ #skewed
    thresh = c(0,1,2,3)
  }
  # autoreg <- rep(autoreg_coeff,max(N)) #leave out random effect for now
  #int<-  rnorm(numSample,gamma_00,gamma_00_sd)
  #autoreg <- rnorm(max(N),autoreg_coeff, 0.569)
  #crosslag <- rep(crosslag_coeff,max(N))

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

  thresh[5]=100
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
  randnum<-floor(runif(1,0,1)*10000)
  if (crosslag_prior == 1){
    set.seed(randnum)
    test<-tryCatch(mod<-{ordinal::clmm2(si_cat_lead ~ si_cat+pred+(1|N), data = datt2, link = "probit"); return(mod)}, warning = function(w) print("Warning"),error = function(e) {print("Error")})

    if (test == "Warning"){
      if (exists("mod") == T){
      sum = summary(mod)
      res<-list(c(sum$coefficients[thresh_length+1,1],sum$coefficients[thresh_length+1,4],0))
      } else {res<-list(c(NA, NA, 2))}
    } else if (test == "Error"){
      res<-list(c(NA, NA, 2))
    } else {
      if (exists("mod") == T){
      sum = summary(mod)
      res<-list(c(sum$coefficients[thresh_length+1,1],sum$coefficients[thresh_length+1,4],1))
      } else {
        res<-list(c(NA, NA, 2))
      }
    }

      #if (test == "Warning"){
       #set.seed(randnum)
       #mod <- ordinal::clmm2(si_cat_lead ~ si_cat+pred+(1|N), data = datt2, link = "probit")
      #  sum <- summary(mod)
       # res<-list(c(sum$coefficients[thresh_length+1,1],sum$coefficients[thresh_length+1,4], 1))
     #} #else if (test == "Error"){
      #  res<-list(c(NA, NA, 2))
       # }
      #  else {
     #     set.seed(randnum)
      #    mod <- ordinal::clmm2(si_cat_lead ~ si_cat+pred+(1|N), data = datt2, link = "probit")
     #   sum <- summary(mod)
     #   res<-list(c(sum$coefficients[thresh_length+1,1],sum$coefficients[thresh_length+1,4],0))
     # }

  } else if (crosslag_prior == 2){
    datt2$si_cat_lead <-as.ordered(datt2$si_cat_lead)
    mod = try(brm(si_cat_lead ~ si_cat+pred+(1|N), data = datt2,family = cumulative("probit", threshold="flexible")))
    sum=summary(mod)
    a<-sum$fixed[thresh_length+1,3]
    b<-sum$fixed[thresh_length+1,4]
    if (a>0 & b > 0) {
      res<-list(c(sum$fixed[thresh_length+1,1],0))
    } else if (a <0 & b<0) {
      res<-list(c(sum$fixed[thresh_length+1,1],0))
    } else {res<-list(c(sum$fixed[thresh_length+1,1],1))}

     }

#coefficients[6,1]: est
  #coefficients[6,4]: p value
  return(res)

}

