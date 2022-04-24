#' Extract Biological descriptor
#'
#'
#' @title extr_biodescr
#' @rdname extr_biodescr
#' @param BasicData BasicData object.
#' @param geneset Charactor vector, one of "kegg"(KEGG), "mkegg"(KEGG Module), "go"(GO-BP), and "wp"(WikiPathways); a data frame and list.
#' @param arguments A list of the arguments of `clusterProfiler`, including `minGSSize`, `maxGSSize`, `pvalue`, and `qvalue`.
#' @param ref_type Charactor vector, one of "drug", "herb", "compound" or "target", defaults to "drug".
#' @param ref Charactor vector, reference drug, herb, compound or target, defaults to `NULL`.
#' @param to_ENTREZID Logical, whether to translate to ENTREZID from SYMBOL, defaults to TRUE.
#' @return A BioDescr object.
#' @author Yuanlong Hu
#' @export

setGeneric("extr_biodescr",
           function(BasicData, geneset= c("kegg", "mkegg","go"),
                    arguments=list(minGSSize=5,maxGSSize = 500,
                                   pvalue=0.05,qvalue=0.1),
                    ref_type="drug", ref=NULL, to_ENTREZID=TRUE){
             standardGeneric("extr_biodescr")
           }

)

#' Subgraph of a Drug graph
#'
#'
#' @title subset_network
#' @rdname subset_network
#' @param BasicData A BasicData object.
#' @param from The source vertex.
#' @param to The target vertex.
#' @return A BasicData object.
#' @export
#' @author Yuanlong Hu

setGeneric("subset_network",
          function(BasicData, from, to=NULL){
              standardGeneric("subset_network")
             })


