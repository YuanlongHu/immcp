#' Create a new list from a data.frame of drug target and disease biomarker as input
#'
#'
#' @title get_list
#' @param disease_biomarker A character. This is disease biomarker.
#' @param drugtarget A data frame of drug target. The first column is drug name, and the second column is drug target.
#' @return list
#' @export
#' @author Yuanlong Hu



get_list <- function(dataframe){
  dataframe <- dataframe[,c(1,2)]
  names(dataframe) <- c("c1", "c2")
  target0 <- list()
  for (i in unique(dataframe[,1])) {
    target1 <- dataframe$c2[dataframe$c1 == i]
    target1 <- list(target1)
    target0 <- c(target0, target1)
  }
  names(target0) <- unique(dataframe[,1])
  return(target0)
}
