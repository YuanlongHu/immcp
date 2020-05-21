##' @method as.data.frame ScoreResult
##' @export


setMethod(f="as.data.frame", signature = "ScoreResult", function(x, ...){
  x <- x@ScoreResult
  as.data.frame(x, ...)
})


##' @method head ScoreResult
##' @export

setMethod(f="head", signature = "ScoreResult", function(x, ...){
  x <- x@ScoreResult
  head(x, ...)
})
