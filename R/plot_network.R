##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreResultNet"),
          function(x, Drug, node_color = c("lightblue","orange", "red", "green"), layout = "layout_nicely", manipulation = FALSE, node_type = "target", background = "drug",neighbor = FALSE) {
            plot_network.ScoreResultNet(x, Drug, node_color = node_color, layout = layout, manipulation = manipulation, node_type = node_type, background = background, neighbor = neighbor)
          })


##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreFP"),
          function(x, Drug, node_type="herb-compound-target", node_color = c("lightblue","orange", "red", "green"), layout = "layout_nicely", manipulation = FALSE, highlight = NULL, width = FALSE) {
            plot_network.ScoreFP(x, Drug, node_type=node_type, node_color = node_color, layout = layout, manipulation = manipulation, highlight = highlight, width = width)
          })


##' @rdname plot_network
##' @exportMethod plot_network

setMethod("plot_network", signature(x = "ScoreResultFP"),
          function(x, Drug, node_type="herb-compound-target", node_color = c("lightblue","orange", "red", "green"), layout = "layout_nicely", manipulation = FALSE, highlight = NULL, width = FALSE) {
            x <- x@Fingerprint
            plot_network.ScoreFP(x = x, Drug, node_type = node_type, node_color = node_color, layout = layout, manipulation = manipulation, highlight = highlight, width = width)
          })


#' @rdname plot_network
#' @param node_type network type. one of "herb-target","herb-compound-target" and "target".
#' @param node_color The node color.
#' @param background one of "drug" or "disease"
#' @param neighbor logical.
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visIgraphLayout



plot_network.ScoreResultNet <- function(x,
                                        Drug,
                                        node_color = c("lightblue","orange", "red", "green"),
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


plot_network.ScoreFP <- function(x,
                                 Drug,
                                 node_type = "herb-compound-pathway",
                                 node_color = c("lightblue","orange", "red", "green"),
                                 layout = "layout_nicely",
                                 manipulation = FALSE,
                                 highlight = NULL,
                                 width = FALSE
                                 ){

  FP1 <- as.data.frame(x@Fingerprint)

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
  val <- unique(paste0(Relationship$col1, Relationship$col2))

  CompoundAnno <- x@CompoundAnno
  Relationship <- x@Relationship

  if (node_type == "herb-compound-pathway"){


    if ("compoundtarget" %in% val & "herbcompound" %in% val & "drugherb" %in% val){

    compound_target <- Relationship[Relationship$col1 == "compound" & Relationship$col2 == "target",]
    compound_target$from <- as.character(compound_target$id)

    compound_target <- to_list(compound_target)

    compound_target <- c(geneset1, compound_target)
    mat_overlap2 <- overlap_count(compound_target)
    mat_overlap2 <- mat_overlap2[mat_overlap2$from %in% names(geneset1) | mat_overlap2$to %in% names(geneset1),]

    # unique(c(mat_overlap2$from, mat_overlap2$to))
    nodes1 <- genesetlist$KEGGPATHID2NAME[genesetlist$KEGGPATHID2NAME$from %in% pathway_overlap,]

    nodes1$color <- node_color[1]
    colnames(nodes1) <- c("id", "label","color")
    nodes2 <- data.frame(from = unique(c(mat_overlap2$from, mat_overlap2$to)),
                         to = unique(c(mat_overlap2$from, mat_overlap2$to))
    )


    nodes2 <- nodes2[!nodes2$from %in% unique(nodes1$from),]


    nodes2 <- merge(nodes2, CompoundAnno, by.x = "from", by.y = "id")
    nodes2 <- nodes2[,-2]
    nodes2$color <- node_color[2]

    colnames(nodes2) <- c("id", "label","color")


    herb_compound <- Relationship[Relationship$col1 == "herb" & Relationship$col2 == "compound",]
    herb_compound <- herb_compound[herb_compound$to %in% nodes2$label,]

    nodes3 <- herb_compound[,1:2]
    nodes3 <- merge(nodes3, CompoundAnno, by.x = "to", by.y = "compound")

    mat_overlap3 <- nodes3[,c(2,3)]
    colnames(mat_overlap3) <- c("from", "to")
    mat_overlap3$width <- 1

    nodes3 <- data.frame(id = unique(nodes3$from), label = unique(nodes3$from))
    nodes3$color <- node_color[3]

    nodes <- Reduce("rbind", list(nodes1, nodes2, nodes3))
    edges <- rbind(mat_overlap2, mat_overlap3)
    }else{
      stop("Lack of information !")
    }

  }


  if (node_type == "herb-pathway"){

    if ("compoundtarget" %in% val & "herbcompound" %in% val & "drugherb" %in% val){

      compound_target <- Relationship[Relationship$col1 == "compound" & Relationship$col2 == "target",]

      herb_compound <- Relationship[Relationship$col1 == "herb" & Relationship$col2 == "compound",]

      herb_target <- merge(herb_compound,compound_target, by.x = "to", by.y = "from")
      herb_target <- herb_target[,c(2,6,3,8,5)]

      colnames(herb_target) <- c("from","to","col1","col2","id")
      herb_target <- herb_target[!duplicated(paste0(herb_target$from, herb_target$to)),]
      # herb_target <- to_list(herb_target)
      # herb_target <- c(geneset1, herb_target)
    }

    if ("herbtarget" %in% val & "drugherb" %in% val){

      drug_herb <- Relationship[Relationship$from == Drug & Relationship$col2 == "herb",]

      herb_target <- Relationship[Relationship$from %in% unique(drug_herb$to) & Relationship$col2 == "target",]
      herb_target <- herb_target[!duplicated(paste0(herb_target$from, herb_target$to)),]

    }

      herb_target <- to_list(herb_target)
      herb_target <- c(geneset1, herb_target)
      edges <- overlap_count(herb_target)
      edges <- edges[edges$from %in% names(geneset1) | edges$to %in% names(geneset1),]

      nodes1 <- genesetlist$KEGGPATHID2NAME[genesetlist$KEGGPATHID2NAME$from %in% pathway_overlap,]

      nodes1$color <- node_color[1]
      colnames(nodes1) <- c("id", "label","color")

      nodes2 <- data.frame(from = unique(c(edges$from, edges$to)),
                         to = unique(c(edges$from, edges$to))
      )
      nodes2 <- nodes2[!nodes2$from %in% unique(nodes1$id),]
      nodes2$color <- node_color[2]

      colnames(nodes2) <- c("id", "label","color")
      nodes <- rbind(nodes1, nodes2)


  }

  if (node_type == "pathway"){

    mat_overlap <- overlap_count(geneset1)
    mat_overlap <- mat_overlap[mat_overlap$from %in% names(geneset1) | mat_overlap$to %in% names(geneset1),]

    nodes1 <- genesetlist$KEGGPATHID2NAME[genesetlist$KEGGPATHID2NAME$from %in% pathway_overlap,]

    nodes1$color <- node_color[1]
    colnames(nodes1) <- c("id", "label","color")

    nodes <- nodes1
    edges <- mat_overlap
  }


  if(!is.null(highlight)){

    highlight_pathway <- lapply(to_list(geneset0), function(x){
      length(intersect(x, highlight))
    })

    highlight_pathway <- unlist(highlight_pathway)
    highlight_pathway <- names(highlight_pathway[highlight_pathway>0])
    nodes$color <- ifelse(nodes$id %in% highlight_pathway, node_color[2], node_color[1])
  }

  message(
    paste("------ Summary ------ \n",
          ">>> Pathway: \n",
          paste(nodes1$label, collapse = ", "),
          "\n",
          ">>> Pathway Number: \n",
          length(nodes1$label))
  )

  if (width) {
    edges$width <- edges$width/5
  }else{
    edges <- edges[,-3]
  }

  if (layout == "none") {
    visNetwork(nodes = nodes,edges = edges) %>%
      visOptions(highlightNearest = TRUE, manipulation = manipulation)
  }else{
    visNetwork(nodes = nodes,edges = edges) %>%
      visOptions(highlightNearest = TRUE, manipulation = manipulation) %>%
      visIgraphLayout(layout = layout)
  }
}


##' Count overlap gene
##'
##' @title overlap_count
##' @param list A list
##' @return a data frame
##' @importFrom corrr shave
##' @importFrom corrr as_cordf
##' @importFrom corrr stretch
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

 b <- shave(as_cordf(b)) %>%
   stretch() %>%
   as("data.frame")

 colnames(b) <- c("from", "to", "width")

 b <- b[!is.na(b$width),]
 b <- b[b$width != 0,]
 b <- as.data.frame(b)

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
##' @param CompoundAnno CompoundAnno.
##' @param nodes A data frame of nodes.
##' @param edges A data frame of edges.
##' @param node_color nodes color.
##' @param node_type network type. one of "herb-target","herb-compound-target" and "target".
##' @param background one of drug or disease.
##' @return a list
##' @noRd
##' @author Yuanlong Hu

create_network <- function(overlap, Drug, Relationship, CompoundAnno, nodes, edges, node_color, node_type, background, neighbor){

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
