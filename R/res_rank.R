##' Rank the results by rank aggregation methods
##'
##'
##' @title res_rank
##' @param ... ScoreResult Object
##' @param method rank aggregation method, by defaylt \code{'RRA'}, other options are
##' \code{'min'}, \code{'geom.mean'}, \code{'mean'}, \code{'median'} and \code{'stuart'}
##' @return  a dataframe with two column
##' @importFrom RobustRankAggreg aggregateRanks
##' @export
##' @author Yuanlong Hu
##' @references Kolde, R., Laur, S., Adler, P., & Vilo, J. (2012). Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics, 28(4), 573-580.


res_rank <- function(..., method = "RRA"){

  a <- list(...)
  a <- lapply(a, function(x){
    res <- get_result(x)
    res <- rownames(res)
  })
  res_rank <- aggregateRanks(a, method = method)
  return(res_rank)
}
