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
#'   data("drugResult")
#'   res <- drugResult$demoScoreFP
#'   res <- get_result(res)

get_result <- function(result, pvalueCutoff = 0.05){
  result <- result@ScoreResult
  result <- result[result$p_value < pvalueCutoff,]
  return(result)
}
