##' @method as.data.frame enrichResult
##' @export
as.data.frame.ScoreResult <- function(x, ...) {
  x <- x@ScoreResult
  as.data.frame(x, ...)
}
