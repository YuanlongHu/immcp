#' Create Disease-Drug Network
#'
#'
#' @title CreateDisDrugNet
#' @param BasicData BasicData object.
#' @param drug Character vector, the drug。
#' @param disease Character vector, the disease.
#' @return A BasicData object.
#' @importFrom dplyr %>%
#' @importFrom igraph induced_subgraph
#' @importFrom igraph all_simple_paths
#' @importFrom igraph neighbors
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom rlang .data
#' @export
#' @author Yuanlong Hu


CreateDisDrugNet <- function(BasicData, drug, disease){

  path <- all_simple_paths(BasicData@drugnet, from = drug, mode = "out") %>%
    lapply(function(x) data.frame(name=x$name, type=x$type))
  path <- Reduce(rbind, path) %>%
    distinct() %>%
    filter(.data$type=="target")

  v_drug <- path$name
  v_dis <- neighbors(BasicData@diseasenet,v = disease, mode = "all")$name
  v <- intersect(v_drug, v_dis)

  # There are the same nodes of target between DrugNet and DisNet.
  DrugNet <- subset_network(BasicData, from = drug, to=v)
  DisNet <- induced_subgraph(BasicData@diseasenet, v)

  v <- list(v)
  names(v) <- disease

  DrugNet@diseasenet <- DisNet
  DrugNet@biomarker <- v

  return(DrugNet)
}


#' Plot Disease-Drug Network
#'
#'
#' @title plot_BasicData
#' @param BasicData BasicData object.
#' @param drug drug
#' @param disease disease
#' @param color Color
#' @param ... Arguments
#' @return Returns NULL, invisibly.
#' @importFrom dplyr %>%
#' @importFrom dplyr recode
#' @importFrom rlang !!!
#' @importFrom igraph induced_subgraph
#' @importFrom igraph plot.igraph
#' @importFrom igraph as.undirected
#' @export
#' @author Yuanlong Hu

plot_BasicData <- function(BasicData,
                           drug=NULL, disease=NULL,
                           color=c(drug="blue",
                                   herb="lightblue",
                                   target="orange"),
                           ...){

  if (!is.null(drug) & !is.null(disease)){
    BasicData <- CreateDisDrugNet(BasicData=BasicData,
                                drug=drug,
                                disease=disease)
  }

  net <- union2(as.undirected(BasicData@drugnet), BasicData@diseasenet)

  # set color
  V(net)$color <- recode(V(net)$type, !!!color)
  ck <- V(net)$color %in% V(net)$type
  if(sum(ck)>0) V(net)$color <- ifelse(ck, "blue", V(net)$color)

  plot.igraph(net, ...)
}


#' Export an xlsx file to Cytoscape.
#'
#'
#' @title exportCytoscape
#' @param BasicData BasicData object.
#' @param file file
#' @return A workbook object
#' @importFrom dplyr %>%
#' @importFrom igraph as_data_frame
#' @importFrom openxlsx write.xlsx
#' @export
#' @author Yuanlong Hu

exportCytoscape <- function(BasicData, file){

  net <- unionnet(BasicData@drugnet, BasicData@diseasenet)
  as_data_frame(net, what = "both") %>%
    write.xlsx(data, file = file)
}
