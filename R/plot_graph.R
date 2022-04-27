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
