#' Calculate the network score
#'
#'
#' @title score_network
#' @param Tar A list containing drug target and disease biomarker.
#' @param DNet A data frame of disease network containing two columns.
#' @param n The number of times random permutation sampling.
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
#' \dontrun{
#'   data("drugSample")
#'   Tar <- drugSample$herb_target
#'   disease_network <- drugSample$disease_network
#'   res <- score_network(Tar, DNet)
#'   res <- as.data.frame(res)
#' }

score_network <- function(Tar, DNet, n = 100){

  if(class(Tar) == "list") Tar
  if(class(Tar) == "data.frame") Tar <- to_list(Tar)

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


  cat("Scoring \n")
  net1 <- pbapply::pblapply(Tar, function(x){

    score <- score_network_s(DNet = DNet, target = x, method = "all")
    set.seed(1234)
    adj_total <- replicate(n, score_network_s(DNet = data.frame(node1 = sample(DNet[,1]), node2 = sample(DNet[,2]), stringsAsFactors = F), target = x, method = "total"))
    res <- c(score, adj_total)

    res
  })

  cat("Summarizing \n")
  result <- pbapply::pblapply(net1, function(x){
    x2 <- x[-c(1:3)]
    z_score <- (x[3] - mean(x2))/sd(x2)
    p_value <- (length(x2[x2 > x[3]])+1)/(n+1)
    res <- c(x[1:3],z_score,p_value)
    names(res) <- c("ChangeDegree", "ChangeDistance", "TotalScore", "adj_TotalScore", "p_value")
    res
  })

  result <- as.data.frame(result) %>%
    t() %>%
    as.data.frame()
  result <- result[order(result$adj_TotalScore, decreasing = T),]
  adj <- lapply(net1, function(x) x[-c(1:2)])
  res_ScoreResult <- new("ScoreResultNet",
                         ScoreResult = as.data.frame(result),
                         DiseaseNetwork = DNet,
                         Tar = Tar,
                         adj = adj)

  return(res_ScoreResult)

}
