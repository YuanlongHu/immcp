#' Create BasicData Object
#'
#'
#' @title CreateBasicData
#' @param ... Drug graph from `PrepareData`.
#' @param diseasenet A graph of Disease-related gene from `PrepareData`.
#' @param biomarker Character vector, the vector of Disease-related gene.
#' @return A BasicData object.
#' @importFrom dplyr %>%
#' @importFrom igraph as_data_frame
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph V
#' @importFrom igraph make_empty_graph
#' @importFrom igraph neighbors
#' @export
#' @author Yuanlong Hu
#' @examples
#' data(drugdemo)
#' drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
#' herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
#' compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
#' disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")
#' BasicData <- CreateBasicData(drug_herb, herb_compound, compound_target, diseasenet = disease)

CreateBasicData <- function(..., diseasenet=NULL, biomarker=NULL){

  net <- union2(...)

  if (is.null(diseasenet)) diseasenet <- make_empty_graph()
  if (is.null(biomarker)){
    # get disease ID
    vertices_disease <- V(diseasenet)$name[V(diseasenet)$type=="disease"]
    disease_list <- as.list(vertices_disease)
    names(disease_list) <- vertices_disease
    biomarker <- lapply(disease_list, function(x){
      biomarker <- neighbors(diseasenet,v=x, mode = "all")$name
      return(biomarker)
    })
  }
  #v <- as_data_frame(net, "vertices")
  res <- new("BasicData",
             drugnet = net,
             #vertices = v,
             diseasenet = diseasenet,
             biomarker = biomarker # a list of biomarker, the disease name
  )
  return(res)
}


#' New union function
#'
#'
#' @title union2
#' @param ... igraph object, the vertices of graph must be the same.
#' @return A igraph object.
#' @importFrom dplyr %>%
#' @importFrom dplyr distinct
#' @importFrom magrittr %<>%
#' @importFrom igraph union
#' @importFrom igraph as_data_frame
#' @importFrom igraph graph_from_data_frame
#' @noRd
#' @author Yuanlong Hu

union2 <- function(...){
  net <- union(...)
  v <- list(...) %>%
    lapply(function(x) as_data_frame(x, "vertices"))
  v <- Reduce(rbind, v) %>% distinct()

  net %<>%
    as_data_frame() %>%
    graph_from_data_frame(directed = TRUE, vertices = v)
  return(net)
}
