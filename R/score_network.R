#' Calculate the network score
#'
#'
#' @title score_network
#' @param targetlist A list containing drug target and disease biomarker.
#' @param disease_network A data frame of disease network containing two columns.
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


score_network <- function(targetlist, disease_network, n = 100){

  score_network_s <- function(disease_network, target, method = "all"){

    disease_network <- as.data.frame(disease_network[,1:2])
    colnames(disease_network)<- c("node1","node2")
    target <- intersect(target, unique(c(disease_network$node1, disease_network$node2)))

    disease_network2 <- disease_network[!disease_network$node1 %in% target,]
    disease_network2 <- disease_network2[!disease_network2$node2 %in% target,]

    g1 <- graph.data.frame(disease_network, directed = F)
    g2 <- graph.data.frame(disease_network2, directed = F)

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


  cat("Calculating \n")
  net1 <- pbapply::pblapply(targetlist, function(x){
    score_network_s(disease_network = disease_network, target = x, method = "all")
  })

  net2 <- as.data.frame(net1) %>%
    t() %>%
    as.data.frame()

  net2 <- data.frame(Drug = names(targetlist),
                     Mean_degree_disturbance_rate = net2$Mean_degree_disturbance_rate,
                     Mean_distance_disturbance_rate = net2$Mean_distance_disturbance_rate,
                     Total_disturbance_rate = net2$Total_disturbance_rate)
  cat("Done \n")

  cat("Adjusting \n")
  res3_list <- pbapply::pblapply(targetlist, function(x){
    degree <- replicate(n, score_network_s(disease_network = data.frame(node1 = sample(disease_network[,1]), node2 = sample(disease_network[,2]), stringsAsFactors = F), target = x, method = "degree"))
    distance <- replicate(n, score_network_s(disease_network = data.frame(node1 = sample(disease_network[,1]), node2 = sample(disease_network[,2]), stringsAsFactors = F),target = x,method = "distance"))
    total <- replicate(n, score_network_s(disease_network = data.frame(node1 = sample(disease_network[,1]), node2 = sample(disease_network[,2]), stringsAsFactors = F), target = x, method = "total"))
    res_r2 <- c(mean(degree), sd(degree),
                mean(distance), sd(distance),
                mean(total), sd(total))
    names(res_r2) <- c("mean_degree", "sd_degree",
                       "mean_distance", "sd_distance",
                       "mean_total", "sd_total")
    return(res_r2)
  })
  res3_list <- as.data.frame(res3_list) %>%
    t() %>%
    as.data.frame()
  res3_list$Drug <- names(targetlist)

  result <- merge(res3_list, net2, by = "Drug")

  result$Mean_degree_adjust <- (result$Mean_degree_disturbance_rate - result$mean_degree)/result$sd_degree
  result$Mean_distance_adjust <- (result$Mean_distance_disturbance_rate - result$mean_distance)/result$sd_distance
  result$Total_network_score_adjust <- (result$Total_disturbance_rate - result$mean_total)/result$sd_total
  result <- result[,-c(2:7)]
  cat("Done \n")


  res_ScoreResult <- new("ScoreResultNet",
                         ScoreResult = as.data.frame(result),
                         DiseaseNetwork = disease_network,
                         Targetlist = targetlist,
                         Adjust = FALSE)

  return(res_ScoreResult)

}
