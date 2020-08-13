#' Calculate the network score
#'
#'
#' @title score_network
#' @param Tar A BasicData object containing drug target.
#' @param DNet A data frame of disease network containing two columns.
#' @param n The number of times random permutation sampling.
#' @param two_tailed a logical: select a two-tailed p-value.
#' @return ScoreResultNet object
#' @importFrom pbapply pblapply
#' @importFrom igraph graph.data.frame
#' @importFrom igraph centr_degree
#' @importFrom igraph mean_distance
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   drug_herb <- PrepareData(drugSample$drug_herb, col1 = "drug", col2 = "herb")
#'   herb_target <- PrepareData(drugSample$herb_target,
#'                              col1 = "herb", col2 = "target",
#'                              format = "basket", sep = ", ")
#'   drug_target <- CreateBasicData(drug_herb, herb_target)
#'   res <- score_network(Tar = drug_target, DNet = drugSample$disease_network)
#'   res <- get_result(res)

score_network <- function(Tar, DNet, n = 100, two_tailed = TRUE){
  Relationship <- Tar@Relationship
  Tar <- Tar@BasicData

  score_network_s <- function(DNet, target, method = "all"){

    DNet <- as.data.frame(DNet[,1:2])
    colnames(DNet)<- c("node1","node2")
    target <- intersect(target, unique(c(DNet$node1, DNet$node2)))

    DNet2 <- DNet[!DNet$node1 %in% target,]
    DNet2 <- DNet2[!DNet2$node2 %in% target,]

    g1 <- graph.data.frame(DNet, directed = F)
    g2 <- graph.data.frame(DNet2, directed = F)

    degree <- (mean(centr_degree(g2)$res) - mean(centr_degree(g1)$res))/mean(centr_degree(g1)$res)
    mean_distance <- (mean_distance(g2, directed = F, unconnected = TRUE) - mean_distance(g1, directed = F, unconnected = TRUE))/mean_distance(g1, directed = F, unconnected = TRUE)
    Total_disturbance_rate <- mean_distance - degree


    if(method == "all"){
      res_network <- c(degree, mean_distance, Total_disturbance_rate)
      names(res_network) <- c("Mean_degree_disturbance_rate", "Mean_distance_disturbance_rate", "Total_disturbance_rate")
      return(res_network)
    }

    if (method == "degree") {
      return(degree)
    }

    if (method == "distance") {
      return(mean_distance)
    }

    if (method == "total") {
      return(Total_disturbance_rate)
    }
  }


  message("Calculating score... \n")
  net1 <- pbapply::pblapply(Tar, function(x){

    score <- score_network_s(DNet = DNet, target = x, method = "all")
    marker <- intersect(x,unique(c(DNet[,1], DNet[,2])))
    marker <- paste0(marker, collapse=", ")
    set.seed(1234)
    adj_total <- replicate(n, score_network_s(DNet = data.frame(node1 = sample(DNet[,1]), node2 = sample(DNet[,2]), stringsAsFactors = F), target = x, method = "total"))
    res <- c(marker, score, adj_total)
    res
  })

  message("Summarizing all results... \n")
  result <- pbapply::pblapply(net1, function(x){
    x1 <- as.numeric(x[4])
    x2 <- as.numeric(x[-c(1:4)])

    z_score <- (x1 - mean(x2))/sd(x2)
    if(two_tailed){
      p_value <- (length(x2[abs(x2) > abs(x1)])+1)/(n+1)
    }else{
      p_value <- (length(x2[x2 > x1])+1)/(n+1)
    }
    p_value <- signif(p_value, 3)
    res <- c(as.numeric(x[2:4]), z_score, p_value, x[1])
    names(res) <- c("ChangeDegree", "ChangeDistance", "TotalScore", "adj_TotalScore", "p_value", "Target")
    res
  })

  result <- as.data.frame(result) %>%
    t() %>%
    as.data.frame()
  result <- result[order(result$adj_TotalScore, decreasing = T),]
  adj <- lapply(net1, function(x) x[-c(1:2)])

  message("Done... \n")
  res_ScoreResult <- new("ScoreResultNet",
                         ScoreResult = as.data.frame(result),
                         DiseaseNetwork = DNet,
                         Tar = Tar,
                         Relationship = Relationship,
                         adj = adj)

  return(res_ScoreResult)
}



##' @rdname imm_centr
##' @exportMethod imm_centr

setMethod("imm_centr", signature(x = "data.frame"),
          function(x) {
            imm_centr.data.frame(x)
          })


##' @rdname imm_centr
##' @exportMethod imm_centr

setMethod("imm_centr", signature(x = "ScoreResultNet"),
          function(x, drug, node = "target", net = "disease") {
            imm_centr.ScoreResultNet(x, drug, node = node, net = net)
          })


#' @rdname imm_centr
#' @importFrom igraph degree
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @importFrom igraph eigen_centrality
#' @importFrom igraph V
#' @author Yuanlong Hu


imm_centr.data.frame <- function(x, node, net){
  x <- as.data.frame(x)[,1:2]
  names(x) <- c("from", "to")
  g <- graph.data.frame(x, directed = F)

  res <- list(
    degree = degree(g, v = V(g), mode = "all"),
    closeness = closeness(g, vids = V(g), mode = "all"),
    betweenness = betweenness(g, v = V(g)),
    eigen = eigen_centrality(g)$vector
  )
  res <- as.data.frame(res)
  return(res)
}

#' @rdname imm_centr
#' @param drug drug name
#' @param node Nodes that need to be evaluated. one of "disease" and "target.
#' @param net Network. one of "disease" and "target.
#' @importFrom igraph degree
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @importFrom igraph eigen_centrality
#' @importFrom igraph V
#' @author Yuanlong Hu

imm_centr.ScoreResultNet <- function(x, drug, node, net){

  g <- x@DiseaseNetwork
  names(g) <- c("from", "to")

  if (net == "target"){
    Target <- x@ScoreResult[drug,"Target"]
    Target <- strsplit(Target, split = ", ")[[1]]
    g <- g[g$from %in% Target & g$to %in% Target,]
  }

  g <- graph.data.frame(g, directed = F)

  res <- list(
    degree = degree(g, v = V(g), mode = "all"),
    closeness = closeness(g, vids = V(g), mode = "all"),
    betweenness = betweenness(g, v = V(g)),
    eigen = eigen_centrality(g)$vector
  )

  res <- as.data.frame(res)
  if (node == "target"){
    Target <- x@ScoreResult[drug,"Target"]
    Target <- strsplit(Target, split = ", ")[[1]]
    res <- res[Target,]

  }


  res <- res[order(res$degree, decreasing = T),]
  return(res)
}
