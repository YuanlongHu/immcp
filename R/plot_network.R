##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreResultNet"),
          function(x, Drug, node_color = c("orange", "lightblue","orange", "green"), layout = "layout_nicely", manipulation = FALSE, node_type = "target", background = "drug",neighbor = FALSE) {
            plot_network.ScoreResultNet(x, Drug, node_color = node_color, layout = layout, manipulation = manipulation, node_type = node_type, background = background, neighbor = neighbor)
          })


##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreFP"),
          function(x, Drug, node_color = c("orange", "lightblue","orange", "green"), layout = "layout_nicely", manipulation = FALSE, highlight = NULL, width = FALSE) {
            plot_network.ScoreFP(x, Drug, node_color = node_color, layout = layout, manipulation = manipulation, highlight = highlight, width = width)
          })


##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreResultFP"),
          function(x, Drug, node_color = c("orange", "lightblue","orange", "green"), layout = "layout_nicely", manipulation = FALSE, highlight = NULL, width = FALSE) {
            x <- x@Fingerprint
            plot_network.ScoreFP(x = x, Drug, node_color = node_color, layout = layout, manipulation = manipulation, highlight = highlight, width = width)
          })


#' @rdname plot_network
#' @param node_type network type. one of "herb-target","herb-compound-target" and "target".
#' @param node_color The node color.
#' @param background one of "drug" or "disease"
#' @param neighbor logical.
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visIgraphLayout
#' @author Yuanlong Hu



plot_network.ScoreResultNet <- function(x,
                                        Drug,
                                        node_color = c("orange", "lightblue","orange", "lightred"),
                                        layout = "layout_nicely",
                                        manipulation = FALSE,
                                        node_type = "target",
                                        background = "drug",
                                        neighbor = FALSE
                                        ){

  Relationship <- x@Relationship
  overlap <- intersect(c(x@DiseaseNetwork[,1], x@DiseaseNetwork[,2]), x@Tar[[Drug]])

  TarNet <- x@DiseaseNetwork

  nodes_df <- data.frame(id = unique(c(TarNet[,1], TarNet[,2])),
                         label = unique(c(TarNet[,1], TarNet[,2]))
  )

  res_df <- create_network(overlap=overlap,
                           Drug = Drug,
                           Relationship=Relationship,
                           nodes=nodes_df, edges=TarNet,
                           node_color=node_color, node_type=node_type,
                           background=background, neighbor=neighbor)

  if (layout == "none") {
    visNetwork(nodes = res_df$nodes,
               edges = res_df$edges) %>%
      visOptions(highlightNearest = TRUE, manipulation = manipulation)
  }else{
    visNetwork(nodes = res_df$nodes,
               edges = res_df$edges) %>%
      visOptions(highlightNearest = TRUE, manipulation = manipulation) %>%
      visIgraphLayout(layout = layout)
  }


}

#' @rdname plot_network
#' @param highlight A character vector of gene.
#' @param width A logical. The number of overlapping genes between the two pathways is used as the width of the edges.
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom corrr shave
#' @importFrom corrr as_cordf
#' @importFrom corrr stretch


plot_network.ScoreFP <- function(x,
                                 Drug,
                                 node_color = c("orange", "lightblue"),
                                 layout = "layout_nicely",
                                 manipulation = FALSE,
                                 highlight = NULL,
                                 width = FALSE
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
    nodes$color <- node_color[2]
  }

  message(
    paste("------ Summary ------ \n",
          ">>> Pathway: \n",
          paste(nodes$label, collapse = ", "),
          "\n",
          ">>> Pathway Number: \n",
          length(nodes$label))
  )

  if (width) {
    mat_overlap$width <- mat_overlap$width/5
  }else{
    mat_overlap <- mat_overlap[,-3]
  }

  if (layout == "none") {
    visNetwork(nodes = nodes,edges = mat_overlap) %>%
      visOptions(highlightNearest = TRUE, manipulation = manipulation)
  }else{
    visNetwork(nodes = nodes,edges = mat_overlap) %>%
      visOptions(highlightNearest = TRUE, manipulation = manipulation) %>%
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

##' Create nodes and edges data frame
##'
##'
##' @title create_network
##' @param overlap overlop gene.
##' @param Drug drug name.
##' @param Relationship Relationship.
##' @param nodes A data frame of nodes.
##' @param edges A data frame of edges.
##' @param node_color nodes color.
##' @param node_type network type. one of "herb-target","herb-compound-target" and "target".
##' @param background one of drug or disease.
##' @return a list
##' @noRd
##' @author Yuanlong Hu

create_network <- function(overlap, Drug, Relationship, nodes, edges, node_color, node_type, background, neighbor){

  if ( !node_type %in% c("herb-target","herb-compound-target","target")){
    stop("The 'node_type' must be one of 'herb-target','herb-compound-target' and 'target'")
  }


  if (background == "drug"){
    nodes$color <- node_color[2]
    nodes <- nodes[nodes$id %in% overlap,]
    edges <- edges[edges[,1] %in% overlap & edges[,2] %in% overlap,]
    names(edges) <- c("from", "to")
  }

  if (background == "disease"){
    nodes$color <- ifelse(nodes$id %in% overlap, node_color[1], node_color[2])
    names(edges) <- c("from", "to")

    if (neighbor){
      edges <- edges[edges$from %in% overlap | edges$to %in% overlap,]
      nodes <- nodes[nodes$label %in% unique(c(edges$from, edges$to)),]
    }
  }


  if (node_type == "herb-target"){
    drug_data <- Relationship[Relationship$col1 == "drug",]
    drug_data <- drug_data[drug_data$from == Drug,]
    target_data <- Relationship[Relationship$col2 == "target",]
    target_data <- target_data[target_data$from %in% unique(drug_data$to),]
    target_data <- target_data[target_data$to %in% unique(nodes$label),]

    nodes_df2 <- data.frame(id = unique(target_data$from),
                            label = unique(target_data$from)
    )
    nodes_df2$color <- node_color[3]
    nodes <- rbind(nodes,nodes_df2)
    edges <- rbind(edges, target_data[,1:2])
  }

  if (node_type == "herb-compound-target"){
    drug_data <- Relationship[Relationship$col1 == "drug",]
    drug_data <- drug_data[drug_data$from == Drug,]
    head(drug_data)

    herb_compound <- Relationship[Relationship$col1 == "herb" & Relationship$col2 == "compound",]
    herb_compound <- herb_compound[herb_compound$from %in% unique(drug_data$to),1:2]
    head(herb_compound)

    compound_target <- Relationship[Relationship$col1 == "compound" & Relationship$col2 == "target",]
    head(compound_target)

    compound_target <- compound_target[compound_target$from %in% unique(herb_compound$to),]
    compound_target <- compound_target[compound_target$to %in% unique(nodes$label),1:2]



    nodes_df2 <- data.frame(id = unique(herb_compound$from),
                            label = unique(herb_compound$from)
    )
    nodes_df2$color <- node_color[3]
    nodes_df3 <- data.frame(id = unique(compound_target$from),
                            label = unique(compound_target$from)
    )
    nodes_df3$color <- node_color[4]
    nodes <- Reduce("rbind", list(nodes, nodes_df2, nodes_df3))

    herb_compound <- herb_compound[herb_compound$to %in% unique(nodes_df3$label),]
    edges <- Reduce("rbind", list(edges, herb_compound, compound_target))
  }


  res_df <- list(nodes = nodes, edges = edges)

  return(res_df)
}
