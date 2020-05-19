##' @method as.data.frame ScoreResult
##' @export


setMethod(f="as.data.frame", signature = "ScoreResult", function(x){
  x <- x@ScoreResult
  as.data.frame(x)
})
