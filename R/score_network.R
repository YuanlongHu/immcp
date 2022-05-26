#' Calculating differences in disease network characteristics before and after removal of drug targets
#'
#'
#' @title score_network
#' @param BasicData A BasicData object.
#' @param n Number vector, the number of times random permutation sampling, default to 1000.
#' @return A list.
#' @importFrom pbapply pblapply
#' @importFrom igraph neighbors
#' @importFrom igraph induced_subgraph
#' @importFrom igraph delete.vertices
#' @importFrom igraph vcount
#' @importFrom igraph as_data_frame
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
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
#' res <- score_network(BasicData, n = 100)

score_network <- function(BasicData, n = 1000){

  #drug_list <- BasicData@vertices %>% filter(.data$type=="drug")

  drug_list <- as_data_frame(BasicData@drugnet, "vertices") %>%
    filter(.data$type=="drug")
  drug_list <- drug_list$name
  g_drug_list <- lapply(as.list(drug_list), function(x){
      BasicData <- subset_network(BasicData = BasicData, from = x)
      #drugtarget <- BasicData@vertices %>% filter(.data$type == "target")
      drugtarget <- as_data_frame(BasicData@drugnet, "vertices") %>%
        filter(.data$type == "target")
      drugtarget$name
  })
  names(g_drug_list) <- drug_list

  dis_list <- V(BasicData@diseasenet)$name[V(BasicData@diseasenet)$type=="disease"]
  g_dis_list <- lapply(as.list(dis_list), function(x){
    v_dis <- neighbors(BasicData@diseasenet, v = x, mode = "all")$name
    induced_subgraph(BasicData@diseasenet, v_dis)
  })
  names(g_dis_list) <- dis_list

  message(">>>>> Calculating Network Characters and Tests <<<<<")
  res_list <- lapply(g_dis_list, function(x){
    res <- pblapply(g_drug_list, function(y){
      del_v <- intersect(V(x)$name, y)
      x2 <- delete.vertices(x, v = del_v)
      if(vcount(x2) == 0) stop("The number of vertices of this graph is 0!")
      res1 <- network_node_ks(graph1 = x, graph2 = x2, replicate = n)
      res2 <- diff_network_char(graph1 = x, graph2 = x2, output_all = TRUE)
      res <- c(res1, res2)
      return(res)
    })
    return(res)
  })

  message(">>>>> Summarizing all results <<<<<")

  res_list <- lapply(res_list, function(x){
    res <- Reduce(rbind, x)
    rownames(res) <- drug_list
    return(res)
  })
  message(">>>>> Done <<<<<")
  return(res_list)
}



#' Calculate the difference of network characters in two network
#'
#'
#' @title diff_network_char
#' @param graph1 A igraph object.
#' @param graph2 A igraph object.
#' @param output_all FALSE
#' @return A number vector.
#' @importFrom igraph graph.data.frame
#' @importFrom igraph delete.vertices
#' @export
#' @author Yuanlong Hu


diff_network_char <- function(graph1, graph2, output_all = FALSE){

  netchar_g1 <- network_char(graph1, T)
  netchar_g2 <- network_char(graph2, T)
  change <- (netchar_g2 - netchar_g1)/netchar_g1

  # Summary
  names(netchar_g1) <- paste0("G1_", names(netchar_g1))
  names(netchar_g2) <- paste0("G2_", names(netchar_g2))
  names(change) <- paste0("Ratio_", names(change))
  change[is.na(change)] <- 0


  if (output_all) {
    res_network <- c(netchar_g1, netchar_g2, change)
  }else{
    res_network <- change
  }
  return(res_network)
}
