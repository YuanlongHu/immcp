#' Drug-Target network visualization
#'
#'
#' @title plot_network
#' @param ResNet A ScoreResultNet object
#' @param Drug A character of drug name
#' @param visIgraphLayout use layout of igraph
#' @return visNetwork
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visIgraphLayout
#' @export
#' @author Yuanlong Hu

plot_network <- function(ResNet, Drug, visIgraphLayout="none"){

  overlap <- intersect(c(ResNet@DiseaseNetwork[,1], ResNet@DiseaseNetwork[,2]), ResNet@Tar[[Drug]])

  TarNet <- ResNet@DiseaseNetwork
  TarNet <- TarNet[TarNet[,1] %in% overlap & TarNet[,2] %in% overlap,]

  if (visIgraphLayout == "none") {
    visNetwork(nodes = data.frame(id = unique(c(TarNet[,1], TarNet[,2])),
                          label = unique(c(TarNet[,1], TarNet[,2]))
    ),
    edges = TarNet) %>%
      visOptions(highlightNearest = TRUE)
  }else{
    visNetwork(nodes = data.frame(id = unique(c(TarNet[,1], TarNet[,2])),
                          label = unique(c(TarNet[,1], TarNet[,2]))
    ),
    edges = TarNet) %>%
      visOptions(highlightNearest = TRUE) %>%
      visIgraphLayout(layout = visIgraphLayout)
  }

}

#' Drug-KEGG Pathyway network visualization
#'
#'
#' @title plot_KEGG_network
#' @param FP A ScoreFP object
#' @param Drug A character of drug name
#' @param layout use layout of igraph
#' @return visNetwork
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom corrr shave
#' @importFrom corrr as_cordf
#' @importFrom corrr stretch
#' @export
#' @author Yuanlong Hu


plot_KEGG_network <- function(FP, Drug, layout = "layout_in_circle"){

  FP1 <- as.data.frame(FP@Fingerprint)

  pathway_overlap <- FP1[,c("disease", Drug)]
  pathway_overlap <- pathway_overlap[pathway_overlap[,1]==1 & pathway_overlap[,2]==1,]
  pathway_overlap <- rownames(pathway_overlap)

  if (FP@Geneset == "KEGG") {
    geneset0 <- genesetlist$KEGGPATHID2EXTID
    geneset0 <- geneset0[,-1]
  }else{
    stop("The geneset for extracting pathway fingerprints must come from KEGG.")
  }
  tar <- FP@DrugTarget[[Drug]]
  geneset0 <- geneset0[geneset0$from %in% pathway_overlap,]
  geneset0 <- geneset0[geneset0$SYMBOL %in% tar,]

  geneset0 <- to_list(geneset0)

  mat_overlap <- overlap_count(geneset0)

  mat_overlap <- shave(as_cordf(mat_overlap)) %>%
    stretch() %>%
    as("data.frame")

  colnames(mat_overlap) <- c("from", "to", "width")

  mat_overlap <- mat_overlap[!is.na(mat_overlap$width),]
  mat_overlap <- mat_overlap[mat_overlap$width != 0,]
  mat_overlap <- as.data.frame(mat_overlap)


  nodes <- genesetlist$KEGGPATHID2NAME[genesetlist$KEGGPATHID2NAME$from %in% unique(c(mat_overlap$from, mat_overlap$to)),]

  colnames(nodes) <- c("id", "label")


  message(
    paste(nodes$label, collapse = ", ")
  )
  visNetwork(nodes = nodes,edges = mat_overlap[,-3]) %>%
    visOptions(highlightNearest = TRUE) %>%
    visIgraphLayout(layout = layout)
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
    res <- genesetlist$KEGGPATHID2NAME[genesetlist$KEGGPATHID2NAME$from %in% res,]
    colnames(res) <- c("ID", "Description")
  }
  return(res)
}

