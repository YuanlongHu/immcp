#' Prepare input format.
#'
#'
#' @title PrepareData
#' @param ... data frame, containing interaction information.
#' @param from A charactor vector, containing "drug", "herb", "compound", or "target".
#' @param to A character vector, containing "drug", "herb", "compound", or "target".
#' @param diseaseID Charactor vector, diseaseID
#' @param format one of "single" or "basket".
#' @param sep Separator.
#' @return A igraph object.
#' @importFrom dplyr %>%
#' @importFrom dplyr distinct
#' @importFrom igraph graph.data.frame
#' @importFrom rlang !!!
#' @export
#' @author Yuanlong Hu
#' @examples
#' data(drugdemo)
#' drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
#' herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
#' compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
#' disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")

PrepareData <- function(..., from, to, diseaseID, format = "single", sep){

  data <- list(...)

  if (! from %in% c("drug","herb","compound","target")) stop("The 'to' must be one of 'drug', 'herb', 'compound', or 'target'")
  if (! to %in% c("drug","herb","compound","target")) stop("The 'from' must be one of 'drug', 'herb', 'compound', or 'target'")

  g <-  lapply(data, function(x){
    if (! is.data.frame(x)) stop("Input data must be a data.frame")

    if (format == "basket") data <- to_list(x[,1:2], input = format, sep = sep) %>% to_df()
    if (format == "single") data <- x[,1:2]
    names(data) <- c("from", "to")
  # create a data frame of vertices
    vertices_from <- data.frame(name= unique(data$from),type= from)
    vertices_to <- data.frame(name= unique(data$to),type= to)
    vertices <- rbind(vertices_from, vertices_to) %>% distinct()

  # create a undirected network
  # from disease to target, from target to target

  if (from == "target" & to == "target") {
    vertices <- rbind(vertices, data.frame(name=diseaseID, type="disease"))
    data <- rbind(data, data.frame(from=diseaseID, to=unique(unlist(data))))
    net <- graph.data.frame(data, directed = FALSE, vertices = vertices)
  }else{
    net <- graph.data.frame(data, directed = TRUE, vertices = vertices)
  }
  return(net)
  })

  if(length(g)>1) {
    g <- union2(!!!g)
  }else{
    g <- g[[1]]
  }
  return(g)
}
