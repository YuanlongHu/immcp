#' Drug target or pathway network visualization
#'
#'
#' @title plot_network
#' @rdname plot_network
#' @param x ScoreFP or ScoreResultNet object
#' @param Drug The name of drug.
#' @param node_color The node color.
#' @param layout Character Name of network layout function to use. Default to "layout_nicely".
#' @param ... additional parameters
#' @return visNetwork object
#' @export

setGeneric("plot_network",
           function(x, Drug, node_color = c("orange", "lightblue"), layout = "layout_nicely", ...){
             standardGeneric("plot_network")
           }

)



#' Calculate the pathway fingerprints
#'
#'
#' @title extrFP
#' @rdname extrFP
#' @param disease_biomarker A character of disease biomarkers or an order ranked geneList.
#' @param drug_target A data frame or list of drug target.
#' @param method one of "enrich" and "gsea"
#' @return ScoreFP object
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   drug_target <- PrepareData(drugSample$herb_target, col1 = "herb", col2 = "target")
#'   drug_target <- CreateBasicData(drug_target)
#'   FP <- extrFP(drug_target = drug_target,
#'                disease_biomarker = drugSample$disease_biomarker,
#'                method = "enrich")


setGeneric("extrFP",
           function(drug_target, disease_biomarker, method = "enrich"){
             standardGeneric("extrFP")
           }

)
