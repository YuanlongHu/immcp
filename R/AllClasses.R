##' Class "ScoreFP"
##' This class represents the pathway fingerprint.
##'
##'
##' @name ScoreFP-class
##' @docType class
##' @slot Fingerprint pathway fingerprint
##' @slot FPType pathway fingerprint type
##' @slot Geneset Geneset name
##' @exportClass ScoreFP
##' @author Yuanlong Hu


setClass("ScoreFP",
         slots = list(
           Fingerprint = "list",
           FPType = "character",
           Geneset = "character"
         ))



setClass("ScoreFP1",
         contains = "ScoreFP",
         slots = list(
           DiseaseBiomarker = "vector",
           DrugTarget = "list",
           Relationship = "data.frame",
           CompoundAnno = "data.frame"
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
##'
##' This class represents the result of score.
##'
##'
##' @name ScoreResult-class
##' @docType class
##' @slot ScoreResult all score reslut.
##' @slot adj distribution data
##' @exportClass ScoreResult
##' @author Yuanlong Hu

setClass("ScoreResult",
         slots = list(
           ScoreResult = "data.frame",
           adj = "list"
         )
)

setClass("ScoreResultFP",
         contains = "ScoreResult",
         slots = list(
           Fingerprint = "ScoreFP"
         )
)


setClass("ScoreResultNet",
         contains = "ScoreResult",
         slots = list(
           DiseaseNetwork = "data.frame",
           Tar = "list",
           Relationship = "data.frame"
         )
)

#' Coerce a ScoreResult object into a data frame
#'
#' @name as.data.frame
#' @aliases as.data.frame,ScoreResult-method
#' @docType methods
#' @param x A ScoreResult object
#' @param row.names NULL or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional logical. If TRUE, setting row names and converting column names (to syntactic names: see make.names) is optional. Note that all of R's base package as.data.frame() methods use optional only for column names treatment, basically with the meaning of data.frame(*, check.names = !optional). See also the make.names argument of the matrix method.
#' @param ... other arguments
#' @export
#' @author Yuanlong Hu



setMethod("as.data.frame", "ScoreResult",
          function(x, row.names=NULL, optional=FALSE, ...){
            data.frame(x@ScoreResult, ...)
          })



#' Return the First Parts of a ScoreResult Object
#'
#' @name head
#' @aliases head,ScoreResult-method
#' @docType methods
#' @param x A ScoreResult object
#' @param ... other arguments
#' @export
#' @author Yuanlong Hu

setMethod("head", "ScoreResult",
          function(x, ...){
            head(x@ScoreResult, ...)
          })


#' Return the last Parts of a ScoreResult Object
#'
#' @name tail
#' @aliases tail,ScoreResult-method
#' @docType methods
#' @param x A ScoreResult object
#' @param ... other arguments
#' @export
#' @author Yuanlong Hu

setMethod("tail", "ScoreResult",
          function(x, ...){
            tail(x@ScoreResult, ...)
          })


##' Class "BasicData"
##' This class represents the basic input data.
##'
##'
##' @name BasicData-class
##' @docType class
##' @slot BasicData Alist containing basic data.
##' @slot Key Column name of basic data.
##' @slot Relationship Relationship.
##' @slot CompoundAnno Compound ID
##' @exportClass BasicData
##' @author Yuanlong Hu

setClass("BasicData",
         slots = list(
           BasicData = "list",
           Key = "character",
           Relationship = "data.frame",
           CompoundAnno = "data.frame"
         ))
