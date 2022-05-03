#‘ extract model parameters
#' This function allows you to extract model parameters from the multilevel models
#' @param result_list results from multilevel models
#' @param nreps number of replications
#' @return the model parameters
#' @importFrom foreach %do%
#' @export


extract_modparams <- function(result_list, nreps) {
  modparams <- matrix()
  modparams<- foreach::foreach (k=1:nreps,.combine='rbind') %do% return(result_list[[k]][[1]])

  colnames(modparams) <- c("res", "res2")

  return(modparams)
}
