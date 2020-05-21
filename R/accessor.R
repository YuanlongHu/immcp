##' as.data.frame method
##'
##' @param x a object of class ScoreResult
##' @author Yuanlong Hu



setMethod(f="as.data.frame", signature = "ScoreResult", function(x){
  as.data.frame(x@ScoreResult)
})


##' head method
##'
##' @param x a object of class ScoreResult
##' @param ... other arguments
##' @author Yuanlong Hu

setMethod(f="head", signature = "ScoreResult", function(x, ...){
  head(x@ScoreResult, ...)
})


##' tail method
##'
##' @param x a object of class ScoreResult
##' @param ... other arguments
##' @author Yuanlong Hu

setMethod(f="tail", signature = "ScoreResult", function(x, ...){
  tail(x@ScoreResult, ...)
})
