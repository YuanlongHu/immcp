#' Subgraph of a Drug graph
#'
#' @rdname subset_network
#' @title subset_network
#' @param BasicData A BasicData object.
#' @param from The source vertex or vertices.
#' @param to The target vertex or vertices.
#' @return  A BasicData object.
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom igraph all_simple_paths
#' @importFrom igraph induced_subgraph
#' @importFrom rlang .data
#' @noRd
#' @author Yuanlong Hu

subset_network <- function(BasicData, from, to=NULL){

  if(is.null(to)) to <- V(BasicData@drugnet)$name

  v <- lapply(as.list(from), function(x){
    all_simple_paths(BasicData@drugnet, from=x,
                     to = to,
                   mode = "out")

  })

  v <- Reduce(c, v) %>%
     lapply(function(x) x$name)%>%
    unlist() %>% unique()

  BasicData@drugnet <- induced_subgraph(BasicData@drugnet, v)
  BasicData@vertices <- BasicData@vertices %>%
                          filter(.data$name %in% V(BasicData@drugnet)$name)
  return(BasicData)
}
