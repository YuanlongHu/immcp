#' Calculate the pathway fingerprint.similarity between disease and prescription.
#'
#'
#' @title score_immpc
#' @param FP A ScoreFP.
#' @param n The number of times random permutation sampling.
#' @return ScoreResult object
#' @importFrom igraph graph.data.frame
#' @importFrom igraph centr_degree
#' @importFrom igraph mean_distance
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#' @importFrom magrittr %>%
#' @export
#' @author Yuanlong Hu
#' @examples
#' \dontrun{
#'   data("drugSample")
#'   FP <- extrFP(disease_biomarker = drugSample$disease_biomarker,
#'                drug_target = drugSample$herb_target,
#'                geneset = "ImmGenTop150")
#'   res <- score_fp(FP, n=100)
#'   res <- as.data.frame(res)
#' }


score_fp <- function(FP, n = 100){

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


  cat("Scoring and Resampling(1/2) \n")
  adj_list <- pbapply::pblapply(f[-1], function(x){
    score_T <- Tanimoto_f(disease = f$disease, drug = x)
    set.seed(123)
    res_r <- replicate(n, Tanimoto_f(disease = sample(f$disease), drug = x))
    res_r <- c(score_T, res_r)
    res_r
  })

  cat("Summarizing(2/2) \n")
  result <- pbapply::pblapply(adj_list, function(x){
    res_r <- x[-1]
    score_T <- x[1]
    adj_score <- (score_T - mean(res_r))/sd(res_r)
    p_value <- length(res_r[res_r>score_T])/n
    p_value <- signif(p_value*2, 3)
    res_r <- c(score_T, adj_score, p_value)
    names(res_r) <- c("Score","adj_score", "p_value")
    res_r
  })

  cat("Done \n")


  result <- as.data.frame(result) %>%
    t() %>%
    as.data.frame()

  result <- result[order(result$adj_score, decreasing = T),]
  head(result)

  res_ScoreResultFP <- new("ScoreResultFP",
                           ScoreResult = as.data.frame(result),
                           Fingerprint = FP,
                           adj = adj_list)

  return(res_ScoreResultFP)
}

