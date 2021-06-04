#' Calculate the natural connectivity
#'
#'
#' @title natural_connectivity
#' @param graph The graph.
#' @return a numeric vector
#' @importFrom igraph as_adjacency_matrix
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#' data("drugSample")
#' graph <- igraph::graph.data.frame(drugSample$disease_network)
#' natural_connectivity(graph)


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
#' @return a data frame
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
#' @importFrom igraph centralization.betweenness
#' @importFrom igraph centralization.degree
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#' data("drugSample")
#' graph <- igraph::graph.data.frame(drugSample$disease_network)
#' network_char(graph)

network_char <- function(graph, total_network=FALSE){

  transitivity <- transitivity(graph, 'local', vids = V(graph))
  transitivity <- ifelse(is.na(transitivity), 0, transitivity)

  if(total_network){

    net_df <- data.frame(
      nodes_num = length(V(graph)),    #number of nodes
      edges_num = length(E(graph)),    #number of edges
      average_degree = mean(degree(graph)),    #average degree
      average_path_length = average.path.length(graph, directed = FALSE),    #average path length
      average_closeness = mean(closeness(graph)),
      average_betweenness = mean(betweenness(graph)),
      average_eigenvector = mean(evcent(graph)$vector),
      average_transitivity = mean(transitivity),
      #nodes_connectivity = vertex.connectivity(graph),    #vertex connectivity
      #edges_connectivity = edge.connectivity(graph),    #dges connectivity
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
    #node_df$transitivity <- ifelse(is.na(node_df$transitivity), 0, node_df$transitivity)
    rownames(node_df) <- V(graph)$name
    return(node_df)
  }
}

#' Kolmogorov-Smirnov tests for node characters between networks
#'
#'
#' @title network_node_ks
#' @param graph list; network node characters
#' @param target drug target
#' @param replicate the number of conduct bootstrapping sampling replications
#' @return a data frame
#' @importFrom pbapply pblapply
#' @importFrom stats ks.test
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#' data("drugSample")
#' graph <- igraph::graph.data.frame(drugSample$disease_network)
#' network_char <- network_char(graph)
#' network_node_ks(network_char)

network_node_ks <- function(graph, target=NULL, replicate=1000){

  #if(is.list(graph) & !is.null(target)) stop("The target must be null!")

  if(is.list(graph) & is.null(target)){
    network_char <- lapply(graph, function(x){
      network_char <- suppressWarnings(network_char(x))
      return(network_char)
    })
  }else{
    target <- intersect(target, V(graph)$name)
    graph2 <- delete.vertices(graph, target)
    graph <- list(graph, graph2)
    network_char <- lapply(graph, function(x){
      network_char <- suppressWarnings(network_char(x))
      return(network_char)
    })
  }


  # bootstrap
  message("----- Bootstrapping Sampling -----")
  network_char_boot <- pblapply(network_char,
                              function(x){
                                set.seed(12345)
                                node_df <- lapply(x, function(y) replicate(replicate, sample(y, 1, replace = TRUE)))
                                node_df <- as.data.frame(node_df)
                                names(node_df) <- c("boot_degree", "boot_closeness",
                                                    "boot_betweenness", "boot_eigenvector",
                                                    "boot_transitivity")
                                return(node_df)
                              })

  # Kolmogorov-Smirnov tests
  message("----- Kolmogorov-Smirnov Tests -----")
  result <- list()
  n <- length(network_char_boot) # a list

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {

      ks <- lapply(as.list(c("boot_degree", "boot_closeness",
                             "boot_betweenness", "boot_eigenvector",
                             "boot_transitivity")),
             function(x){
               node_i <- network_char_boot[[i]][ ,x]
               node_j <- network_char_boot[[j]][ ,x]
               ks <- suppressWarnings(ks.test(node_i, node_j, alternative = 'two.sided'))
               return(ks$p.value)
             })

      # Summary p-value
      result <- c(result,
                  list(
                        c("Network1" = names(network_char_boot)[i],
                          "Network2" = names(network_char_boot)[j],
                          "ks_degree" = ks[[1]],
                          "ks_closeness" = ks[[2]],
                          "ks_betweenness" = ks[[3]],
                          "ks_eigenvector" = ks[[4]],
                          "ks_transitivity" = ks[[5]])
                        )
                  )
    }
  }

  message("----- Summary Result -----")

  if(length(result)==1){
      result <- result[[1]]
      result <- as.data.frame(result)
      result <- as.data.frame(t(result))
    }else{
      result <- Reduce(rbind, result)
      rownames(result) <- 1:nrow(1)
    }

  #result <- as.data.frame(result)
  return(result)    # p<0.05
}
