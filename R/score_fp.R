#' Calculate the similarity of immune fingerprint and the characteristics of network topology
#'
#'
#' @title score_immpc
#' @param disease_biomarker A character.
#' @param disease_network A data frame of disease network.This data frame containing a symbolic edge list in the first two columns.
#' @param target A data frame of drug target.
#' @param geneset A list.
#' @return ScoreResult object
#' @importFrom igraph graph.data.frame
#' @importFrom igraph centr_degree
#' @importFrom igraph mean_distance
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#' @export
#' @author Yuanlong Hu


score_fp <- function(FP,n=100){

  f <- FP@Fingerprint
  target <- FP@DrugTarget

  Tanimoto_f <- function(disease, drug){
    f <- rbind(disease, drug)
    res_f <- as.matrix(f) %>%
      Matrix::Matrix() %>%
      proxyC::simil(method = "jaccard") %>%
      as.matrix() %>%
      as.data.frame()
    res_f <- res_f[1,2]
    return(res_f)
  }

  cat("Calculating \n")

  f1 <- as.data.frame(f) %>%
    t() %>%
    Matrix::Matrix() %>%
    proxyC::simil(method = "jaccard") %>%
    as.matrix() %>%
    as.data.frame()
  f1 <- data.frame(Drug = names(f[-1]), Tanimoto = f1[-1,1], stringsAsFactors = F)
  cat("Done \n")

  cat("Adjusting \n")
  res2_list <- pbapply::pblapply(f[-1], function(x){
    set.seed(123)
    res_r <- replicate(n, Tanimoto_f(disease = sample(f$disease), drug = x))
    res_r <- c(mean(res_r), sd(res_r))
    names(res_r) <- c("mean", "sd")
    res_r
  })

  res2_list <- as.data.frame(res2_list) %>%
    t() %>%
    as.data.frame()

  res2_list <- data.frame(Drug = names(f)[-1],
                          mean= res2_list$mean,
                          sd = res2_list$sd)


  res2_list <- merge(res2_list, f1, by="Drug")
  res2_list$Tanimoto_adjust <- (res2_list$Tanimoto - res2_list$mean)/ res2_list$sd

  cat("Done \n")

  result <- res2_list[,-c(2,3)]

  res_ScoreResultFP <- new("ScoreResultFP",
                       ScoreResult = as.data.frame(result),
                       Fingerprint = FP,
                       Adjust = FALSE)

  return(res_ScoreResultFP)
}

