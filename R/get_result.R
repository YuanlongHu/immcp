#' Extract a table of the score result
#'
#'
#' @title get_result
#' @param result an object of class ScoreResult.
#' @param pvalueCutoff p-value cutoff.
#' @return a data.frame
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   FP <- extrFP(disease_biomarker = drugSample$disease_biomarker,
#'                drug_target = drugSample$herb_target,
#'                geneset = "ImmGenTop150")
#'   res <- score_fp(FP, n=100)
#'   res <- get_result(res)

get_result <- function(result, pvalueCutoff = 0.05){
  result <- result@ScoreResult
  result <- result[result$p_value < pvalueCutoff,]
  return(result)
}
