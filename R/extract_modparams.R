#‘ extract model parameters
#' This function allows you to extract model parameters from the multilevel models
#' @param result_list results from multilevel models
#' @param reps number of replications
#' @return the model parameters
#' @importFrom foreach %do%
#' @export


extract_modparams <- function(result_list, reps) {
  modparams <- matrix()
  modparams<- foreach::foreach (i=1:reps,.combine='rbind') %do% return(result_list[[i]][[1]])

  colnames(modparams) <- c("res", "res2", "res3")

  return(modparams)
}
