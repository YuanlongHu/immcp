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
