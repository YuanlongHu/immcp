##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreResultNet"),
          function(x, Drug, node_color = c("red", "blue"),
                   layout = "layout_nicely", ...) {
            plot_network.ScoreResultNet(x, Drug, node_color = node_color, layout = layout, ...)
          })


##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreFP"),
          function(x, Drug, node_color = c("red", "blue"),
                   layout = "layout_nicely", ...) {
            plot_network.ScoreFP(x, Drug, node_color = node_color, layout = layout, ...)
          })


#' @rdname plot_network
#' @param node one of "target" or "all"
#' @param node_color The node color
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visIgraphLayout
#' @author Yuanlong Hu



plot_network.ScoreResultNet <- function(x,
                                        Drug,
                                        node_color = c("red", "blue"),
                                        layout = "layout_nicely",
                                        node = "target"
                                        ){

  overlap <- intersect(c(x@DiseaseNetwork[,1], x@DiseaseNetwork[,2]), x@Tar[[Drug]])

  TarNet <- x@DiseaseNetwork

  nodes_df <- data.frame(id = unique(c(TarNet[,1], TarNet[,2])),
                         label = unique(c(TarNet[,1], TarNet[,2]))
  )
  if (node == "target"){
    nodes_df$color <- rep(node_color[2], nrow(nodes_df))
    nodes_df <- nodes_df[nodes_df$id %in% overlap,]
    TarNet <- TarNet[TarNet[,1] %in% overlap & TarNet[,2] %in% overlap,]

  }

  if (node == "all"){
    nodes_df$color <- ifelse(nodes_df$id %in% overlap, node_color[1], node_color[2])
  }


  names(TarNet) <- c("from", "to")

  if (layout == "none") {
    visNetwork(nodes = nodes_df,
    edges = TarNet) %>%
      visOptions(highlightNearest = TRUE)
  }else{
    visNetwork(nodes = nodes_df,
    edges = TarNet) %>%
      visOptions(highlightNearest = TRUE) %>%
      visIgraphLayout(layout = layout)
  }

}

#' @rdname plot_network
#' @param highlight A character vector of gene.
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom corrr shave
#' @importFrom corrr as_cordf
#' @importFrom corrr stretch


plot_network.ScoreFP <- function(x,
                                 Drug,
                                 node_color = c("blue", "red"),
                                 layout = "layout_nicely",
                                 highlight = NULL
                                 ){

  FP1 <- as.data.frame(x@Fingerprint)

  # The intersection of disease and drug pathways
  pathway_overlap <- FP1[,c("disease", Drug)]
  pathway_overlap <- pathway_overlap[pathway_overlap[,1] == 1 & pathway_overlap[,2] == 1,]
  pathway_overlap <- rownames(pathway_overlap)

  if (x@Geneset == "KEGG") {
    geneset0 <- genesetlist$KEGGPATHID2EXTID
    geneset0 <- geneset0[,-1]
  }else{
    stop("The geneset for extracting pathway fingerprints must come from KEGG.")
  }
  tar <- x@DrugTarget[[Drug]]
  geneset0 <- geneset0[geneset0$from %in% pathway_overlap,]
  geneset1 <- geneset0[geneset0$SYMBOL %in% tar,]

  geneset1 <- to_list(geneset1)

  # Number of overlapping genes between pathways
  mat_overlap <- overlap_count(geneset1)
  mat_overlap <- shave(as_cordf(mat_overlap)) %>%
    stretch() %>%
    as("data.frame")

  colnames(mat_overlap) <- c("from", "to", "width")

  mat_overlap <- mat_overlap[!is.na(mat_overlap$width),]
  mat_overlap <- mat_overlap[mat_overlap$width != 0,]
  mat_overlap <- as.data.frame(mat_overlap)

  nodes <- genesetlist$KEGGPATHID2NAME[genesetlist$KEGGPATHID2NAME$from %in% unique(c(mat_overlap$from, mat_overlap$to)),]

  colnames(nodes) <- c("id", "label")
  if(!is.null(highlight)){

    highlight_pathway <- lapply(to_list(geneset0), function(x){
    length(intersect(x, highlight))
  })

    highlight_pathway <- unlist(highlight_pathway)
    highlight_pathway <- names(highlight_pathway[highlight_pathway>0])
    nodes$color <- ifelse(nodes$id %in% highlight_pathway, node_color[2], node_color[1])
  }else{
    nodes$color <- rep(node_color[2], nrow(nodes))
  }

  message(
    paste("------ Summary ------ \n",
          "> Pathway: \n",
          paste(nodes$label, collapse = ", "),
          "\n",
          "> Pathway Number: \n",
          length(nodes$label))
  )

  if (layout == "none"){
    visNetwork(nodes = nodes,edges = mat_overlap[,-3]) %>%
      visOptions(highlightNearest = TRUE)
  }else{
    visNetwork(nodes = nodes,edges = mat_overlap[,-3]) %>%
      visOptions(highlightNearest = TRUE) %>%
      visIgraphLayout(layout = layout)
  }
}


##' Count overlap gene
##'
##' @title overlap_count
##' @param list A list
##' @return a data frame
##' @author Yuanlong Hu
##' @noRd

overlap_count <- function(list){

  b <- lapply(list, function(x){
    a <- x
    b <- lapply(list, function(x){
      length(intersect(x,a))
    })

    b <- unlist(b)
    return(b)
  })

 b <-  as.data.frame(b)

 return(b)
}



##' Performs set intersection on pathyways fingerprints
##'
##'
##' @title overlap_pathway
##' @param FP A ScoreFP Object
##' @param Drug The drug names
##' @return a vector or data frame
##' @export
##' @author Yuanlong Hu

overlap_pathway <- function(FP, Drug){
  Drug <- c("disease", Drug)
  FP_d <- FP@Fingerprint[Drug]
  FP_d <- lapply(FP_d, function(x){
    a <- names(x)[x==1]
    return(a)
  })

  res <- Reduce(intersect,FP_d)
  if (FP@Geneset == "KEGG") {
    geneset <- genesetlist$KEGGPATHID2NAME
    res <- geneset[geneset$from %in% res,]
    colnames(res) <- c("ID", "Description")
    rownames(res) <- res$ID
  }
  return(res)
}

