#' Create Disease-Drug Network
#'
#'
#' @title CreateDisDrugNet
#' @param BasicData BasicData object.
#' @param drug Character vector, the drug。
#' @param disease Character vector, the disease.
#' @return A igraph object.
#' @importFrom dplyr %>%
#' @importFrom igraph induced_subgraph
#' @importFrom igraph neighbors
#' @importFrom igraph delete.vertices
#' @importFrom igraph as.undirected
#' @importFrom dplyr filter
#' @importFrom igraph subcomponent
#' @importFrom rlang .data
#' @export
#' @author Yuanlong Hu
#' @examples
#' data(drugdemo)
#' drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
#' herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
#' compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
#' disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")
#' BasicData <- CreateBasicData(drug_herb, herb_compound, compound_target, diseasenet = disease)
#' DisDrugNet <- CreateDisDrugNet(BasicData, drug = "Drug1", disease = "disease")

CreateDisDrugNet <- function(BasicData, drug, disease){

  v_drug <- lapply(as.list(drug), function(x){
    path <- subcomponent(BasicData@drugnet, v = x, mode = "out")
    path <-  data.frame(name = path$name, type = path$type) %>%
      filter(.data$type == "target")
    path <- path$name
    return(path)
  }) %>% unlist() %>% unique()

  v_dis <- neighbors(BasicData@diseasenet, v = disease, mode = "all")$name
  v <- intersect(v_drug, v_dis)

  # There are the same nodes of target between DrugNet and DisNet.
  DrugNet <- subset_network(BasicData, from = drug, to = v)
  DisNet <- induced_subgraph(BasicData@diseasenet, v)

  if (length(drug) == 1) {
    drugnet <- delete.vertices(DrugNet@drugnet, v = drug)
    DisDrugNet <- union2(DisNet, as.undirected(drugnet)) %>%
      as.undirected()

  }else{
    DisDrugNet <- union2(DisNet, as.undirected(DrugNet@drugnet))%>%
      as.undirected()
  }

  return(DisDrugNet)
}

#' Plot Disease-Drug Network
#'
#' @rdname plot_graph-method
#' @title Plot Disease-Drug Network
#' @param graph graph.
#' @param drug drug.
#' @param disease disease.
#' @param vis one of "igraph", "visNetwork" and "shiny".
#' @param color Color
#' @param width width
#' @param ... Arguments
#' @return Returns NULL, invisibly.
#' @exportMethod plot_graph

setMethod("plot_graph", signature(graph = "BasicData"),
          function(graph,
                   drug, disease,
                   vis = "visNetwork",
                   color= c(drug="blue",
                            herb="lightblue",
                            target="orange"),
                   width = 1,
                   ...){

            DisDrugNet <- CreateDisDrugNet(BasicData = graph,
                                          drug = drug, disease = disease)
            plot_graph_internal(graph=DisDrugNet, vis = vis, color=color, width = width, ...)
          })

#' @rdname plot_graph-method
#' @exportMethod plot_graph

setMethod("plot_graph", signature(graph = "igraph"),
          function(graph,
                   vis = "visNetwork",
                   color= c(drug = "blue",
                            herb = "lightblue",
                            target = "orange"),
                   width = 1, ...){

            plot_graph_internal(graph, vis = vis, color=color, width = width, ...)
          })

#' @importFrom dplyr %>%
#' @importFrom dplyr recode
#' @importFrom rlang !!!
#' @importFrom igraph induced_subgraph
#' @importFrom igraph plot.igraph
#' @importFrom igraph as.undirected
#' @importFrom igraph E
#' @importFrom igraph E<-
#' @importFrom igraph V<-
#' @importFrom igraph V
#' @importFrom visNetwork toVisNetworkData
#' @importFrom visNetwork visNetwork
#' @importFrom visNetwork visOptions
#' @importFrom visNetwork visEdges
#' @importFrom visNetwork visNetworkEditor

plot_graph_internal <- function(graph,
                       vis = "visNetwork",
                       color=c(drug = "blue",
                               herb = "lightblue",
                               target = "orange"),
                      width = 1,
                           ...){

  # set color
  net <- graph
  V(net)$color <- recode(V(net)$type, !!!color)
  ck <- V(net)$color %in% V(net)$type
  if(sum(ck)>0) V(net)$color <- ifelse(ck, "blue", V(net)$color)

  if (vis == "igraph") {
    E(net)$width <- width
    plot.igraph(net, ...)
    }

  if (vis == "visNetwork"){
    data_visNetwork <- toVisNetworkData(net)
    data_visNetwork$edges$width <- width
    p <- visNetwork(nodes = data_visNetwork$nodes,
               edges = data_visNetwork$edges,
               height = "700px", width = "100%") %>%
      visOptions(selectedBy = "type",
                 manipulation = TRUE,
                 highlightNearest = TRUE) %>%
      visEdges(smooth = FALSE)
    return(p)
  }

  if (vis == "shiny"){
    data_visNetwork <- toVisNetworkData(net)
    data_visNetwork$edges$width <- width
    visNetwork(nodes = data_visNetwork$nodes,
               edges = data_visNetwork$edges,
               height = "700px", width = "100%") %>%
      visOptions(selectedBy = "type",
                 manipulation = TRUE,
                 highlightNearest = TRUE)%>%
      visEdges(smooth = FALSE) %>%
      visNetworkEditor()

  }
}


#' Export an xlsx file to Cytoscape.
#'
#'
#' @rdname exportCytoscape
#' @title Export an xlsx file to Cytoscape
#' @param graph igraph object.
#' @param file file
#' @return A workbook object
#' @importFrom dplyr %>%
#' @importFrom igraph as_data_frame
#' @importFrom openxlsx write.xlsx
#' @export
#' @author Yuanlong Hu

exportCytoscape <- function(graph, file){

  # net <- union2(as.undirected(BasicData@drugnet), BasicData@diseasenet)
  as_data_frame(graph, what = "both") %>%
    write.xlsx(data, file = file)
}
