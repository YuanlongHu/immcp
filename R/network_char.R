#' Calculate the natural connectivity
#'
#'
#' @title natural_connectivity
#' @param graph A igraph object.
#' @return A numeric vector.
#' @importFrom igraph as_adjacency_matrix
#' @export
#' @author Yuanlong Hu


natural_connectivity <- function(graph) {

  adj_matrix <- as_adjacency_matrix(graph,sparse=F)
  adj_matrix[abs(adj_matrix) != 0] <- 1

  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)

  lambda_sum <- 0
  N = length(lambda)
  for (i in 1:N) lambda_sum <- lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum/N, base = exp(1))
  return(lambda_average)
}

#' Calculate the network characters
#'
#'
#' @title network_char
#' @param graph The graph.
#' @param total_network Calculate for total network or each nodes.
#' @return A number vector or data frame.
#' @importFrom igraph degree
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @importFrom igraph evcent
#' @importFrom igraph transitivity
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph vertex.connectivity
#' @importFrom igraph edge.connectivity
#' @importFrom igraph average.path.length
#' @importFrom igraph diameter
#' @importFrom igraph graph.density
#' @importFrom igraph vcount
#' @importFrom igraph centralization.betweenness
#' @importFrom igraph centralization.degree
#' @export
#' @author Yuanlong Hu

network_char <- function(graph, total_network=FALSE){
  if(vcount(graph)==0) stop("The number of vertices of this graph is 0!")
  transitivity <- transitivity(graph, 'local', vids = V(graph))
  transitivity <- ifelse(is.na(transitivity), 0, transitivity)

  if(total_network){

    net_df <- c(
      nodes_num = length(V(graph)),    #number of nodes
      edges_num = length(E(graph)),    #number of edges
      average_degree = mean(degree(graph)),    #average degree
      average_path_length = average.path.length(graph, directed = FALSE),    #average path length
      average_closeness = mean(closeness(graph)),
      average_betweenness = mean(betweenness(graph)),
      average_eigenvector = mean(evcent(graph)$vector),
      average_transitivity = mean(transitivity),
      natural_connectivity = natural_connectivity(graph), #natural connectivity
      graph_diameter = diameter(graph, directed = FALSE),    #diameter
      graph_density = graph.density(graph),    #density
      clustering_coefficient = transitivity(graph),    #clustering coefficient
      betweenness_centralization = centralization.betweenness(graph)$centralization,    #betweenness centralization
      degree_centralization = centralization.degree(graph)$centralization    #degree centralization.
    )
    return(net_df)
  }else{
    node_df <- data.frame(
      degree = degree(graph),
      closeness_centrality = closeness(graph),
      betweenness_centrality = betweenness(graph),
      eigenvector_centrality = evcent(graph)$vector,
      transitivity = transitivity
    )
    rownames(node_df) <- V(graph)$name
    return(node_df)
  }
}


