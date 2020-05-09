##' Class "ScoreFP"
##' This class represents the FP.
##'
##'
##' @name ScoreFP-class
##' @docType class
##'
##'
##' @exportClass ScoreFP
##' @author Yuanlong Hu


setClass("ScoreFP",
         slots = list(
           Fingerprint = "list",
           FPType = "character"
         ))



setClass("ScoreFP1",
         contains = "ScoreFP",
         slots = list(
           DiseaseBiomarker = "character",
           DrugTarget = "list"
         )
)


setClass("ScoreFP2",
         contains = "ScoreFP",
         slots = list(
           DiseaseExpr = "data.frame",
           DiseaseDEG = "data.frame",
           DrugExpr = "list",
           DrugDEG = "list"
         )
)

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
           Fingerprint = "ScoreFP",
           DiseaseNetwork = "data.frame",
           Adjust = "logical"
         )
)
