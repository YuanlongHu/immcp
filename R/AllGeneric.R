#' @title Extract Biological descriptor
#' @rdname extr_biodescr-method
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
           function(BasicData, geneset= c("kegg", "mkegg","go","wp"),
                    arguments=list(minGSSize=5,maxGSSize = 500,
                                   pvalue=0.05,qvalue=0.1),
                    ref_type="drug", ref=NULL, to_ENTREZID=TRUE){
             standardGeneric("extr_biodescr")
           }

)




#' @title Plot Disease-Drug Network
#' @rdname plot_graph-method
#' @param graph graph.
#' @param Isolated  Whether to delect Isolated nodes.
#' @param drug drug.
#' @param disease disease.
#' @param vis one of "igraph", "visNetwork" and "shiny".
#' @param color Nodes Color
#' @param width Edges width
#' @param size Nodes size
#' @param ... Arguments
#' @return Returns NULL, invisibly.
#' @author Yuanlong Hu
#' @export

setGeneric("plot_graph",
           function(graph,
                    drug, disease,Isolated = TRUE,
                    vis = "visNetwork",
                    color = c(drug = "blue",
                            herb = "lightblue",
                            target = "orange"),
                    width = 1, size = 20, ...){
             standardGeneric("plot_graph")
           }

)

