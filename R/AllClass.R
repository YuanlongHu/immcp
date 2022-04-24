# register the S3 `igraph` class for use with S4 methods.
#' @import methods
setOldClass("igraph", igraph::make_empty_graph())
selectMethod("show", "igraph")

#' Class \code{BasicData}
#' This class represents the basic input data.
#'
#'
#' @name BasicData-class
#' @docType class
#' @slot drugnet A directed graph
#' @slot vertices Vertices of drug graph.
#' @slot diseasenet Disease network.
#' @slot biomarker Disease-related gene.
#' @exportClass BasicData
#' @author Yuanlong Hu

setClass("BasicData",
         slots = list(
           drugnet = "igraph",
           vertices = "data.frame",
           diseasenet = "igraph",
           biomarker = "list"
         ))

#' Class \code{BioDescr}
#' This class represents the biological descriptor data.
#'
#'
#' @name BioDescr-class
#' @docType class
#' @slot drug_geneset from drug to geneset.
#' @slot geneset_gene from geneset to gene for each drug.
#' @slot anno Geneset ID and description.
#' @slot BasicData BasicData object.
#' @exportClass BioDescr
#' @author Yuanlong Hu

setClass("BioDescr",
         slots = list(
           drug_geneset = "igraph",
           geneset_gene = "list",
           anno = "data.frame",
           BasicData = "BasicData"
         ))
