#' @rdname plot_graph-method
#' @exportMethod plot_graph

setMethod("plot_graph", signature(graph = "BasicData"),
          function(graph,
                   drug, disease,
                   vis = "visNetwork",
                   color= c(drug="blue",
                            herb="lightblue",
                            target="orange"),
                   width = 1, size = 1, ...){

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
                   width = 1, size=1, ...){

            plot_graph_internal(graph, vis = vis, color = color, width = width, ...)
          })

#' @rdname plot_graph-method
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph delete.vertices
#' @importFrom igraph degree
#' @exportMethod plot_graph

setMethod("plot_graph", signature(graph = "HerbResult"),
          function(graph, Isolated = TRUE,
                   vis = "visNetwork",
                   color= "lightblue",
                   width = 1, size = 20,...){

            graph <- graph@Herb_Herb
            if(Isolated) graph <- delete.vertices(graph, which(degree(graph)==0))

            if (size == "frequency") size <- V(graph)$frequency_relative*30

            if (width == "support") width <- E(graph)$support*10
            if (width == "confidence") width <- E(graph)$confidence*10
            if (width == "coverage") width <- E(graph)$coverage*10
            if (width == "lift") width <- E(graph)$lift*1.1

            plot_graph_internal(graph, vis = vis, color=c(herb = color), width = width, size = size ,...)
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
                                size = 1,
                                ...){

  # set color
  net <- graph
  V(net)$color <- recode(V(net)$type, !!!color)
  ck <- V(net)$color %in% V(net)$type
  if(sum(ck)>0) V(net)$color <- ifelse(ck, "blue", V(net)$color)

  if (vis == "igraph") {
    E(net)$width <- width
    V(net)$size <- size
    plot.igraph(net, ...)
  }

  if (vis == "visNetwork"){
    data_visNetwork <- toVisNetworkData(net)
    if(nrow(data_visNetwork$edges)>0) data_visNetwork$edges$width <- width
    if(nrow(data_visNetwork$nodes)>0) data_visNetwork$nodes$size <- size
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
    if(nrow(data_visNetwork$edges)>0) data_visNetwork$edges$width <- width
    if(nrow(data_visNetwork$nodes)>0) data_visNetwork$nodes$size <- size
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
