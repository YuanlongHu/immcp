
setMethod(f="as.data.frame", signature = "ScoreResult", function(x){
  as.data.frame(x@ScoreResult)
})


setMethod(f="head", signature = "ScoreResult", function(x, ...){
  head(x@ScoreResult, ...)
})


setMethod(f="tail", signature = "ScoreResult", function(x, ...){
  tail(x@ScoreResult, ...)
})
