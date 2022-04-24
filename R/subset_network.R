##' @rdname subset_network
##' @exportMethod subset_network

setMethod("subset_network", signature(BasicData = "BasicData"),
          function(BasicData, from="drug", to = NULL) {
            subset_network.BasicData(BasicData=BasicData, from=from, to = to)
          })


#' Subgraph of a Drug graph
#'
#'
#' @title subset_network
#' @param BasicData A BasicData object.
#' @param from The source vertex.
#' @param to The target vertex.
#' @return  a BasicData object.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom igraph all_simple_paths
#' @importFrom igraph induced_subgraph
#' @importFrom rlang .data
#' @export
#' @author Yuanlong Hu

subset_network.BasicData <- function(BasicData, from, to=NULL){

  if(is.null(to)) to <- V(BasicData@drugnet)$name
  path <- all_simple_paths(BasicData@drugnet, from = from, to = to, mode = "out") %>%
    lapply(function(x) x$name) %>%
    unlist() %>% unique()

  BasicData@drugnet <- induced_subgraph(BasicData@drugnet, path)
  BasicData@vertices <- BasicData@vertices %>%
                          filter(.data$name %in% V(BasicData@drugnet)$name)
  return(BasicData)
}


#' Union of graphs created by `PrepareData`
#'
#'
#' @title unionnet
#' @param ... Drug graph from `PrepareData`.
#' @return A new igraph object.
#' @importFrom dplyr %>%
#' @importFrom igraph union
#' @importFrom igraph V
#' @importFrom igraph V<-
#' @export
#' @author Yuanlong Hu


unionnet <- function(...){

  v <-  list(...) %>% lapply(function(x) as_data_frame(x, "vertices"))
  v <- Reduce(rbind, v) %>% distinct()
  diseasenet <- union(...)
  v <- v[V(diseasenet)$name,]
  V(diseasenet)$type <- v$type
  return(diseasenet)
}
