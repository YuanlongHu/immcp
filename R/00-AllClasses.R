##' Class "ScoreResult"
##' This class represents the score result.
##'
##'
##' @name ScoreResult-class
##' @docType class
##'
##'
##' @exportClass ScoreResult
##' @author Yuanlong Hu

setClass("ScoreResult",
         slots = list(
           ScoreResult = "data.frame",
           DiseaseNetwork = "data.frame",
           DiseaseBiomarker = "character",
           Target = "data.frame",
           Adjust = "logical"
         )
)
