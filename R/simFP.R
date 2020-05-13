#' Calculate the similarity between Drug Fingerprints
#'
#'
#' @title simFP
#' @param FP A ScoreFP object
#' @return a matrix
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#' @export
#' @author Yuanlong Hu



simFP <- function(FP){
 FP <- FP@Fingerprint[-1]
 f1 <- as.data.frame(FP) %>%
   t() %>%
   Matrix::Matrix() %>%
   proxyC::simil(method = "jaccard") %>%
   as.matrix()

 return(f1)
}


#' Visualization of a drug fingerprints similarity matrix using ggplot2
#'
#'
#' @title Heatmap_simFP
#' @param data A drug fingerprints similarity matrix
#' @return a ggplot2
#' @importFrom ggcorrplot ggcorrplot
#' @export
#' @author Yuanlong Hu


Heatmap_simFP <- function(data, ...){
  ggcorrplot(data, ...)
}

