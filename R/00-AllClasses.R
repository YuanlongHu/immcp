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
           Fingerprint = "list",
           DiseaseNetwork = "data.frame",
           DiseaseBiomarker = "character",
           Target = "list",
           Adjust = "logical"
         )
)
