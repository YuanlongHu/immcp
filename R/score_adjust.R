#' Adjust Score for Permutation Test
#'
#' @title score_adjust
#' @param result ScoreResult object
#' @param n integer.
#'

#setMethod("score_adjust",
#          signature(result = "ScoreResult", n = "integer"),
score_adjust <- function(result, n = 100){

  #disease <- result@Fingerprint$disease
  #drug <- result@Fingerprint[-1]

  if(result@Adjust == TRUE){
    stop("This result has been adjusted !")
  }

  Tanimoto_f <- function(disease, drug){
    f <- rbind(disease, drug)
    res_f <- as.matrix(f) %>%
      Matrix::Matrix() %>%
      proxyC::simil(method = "jaccard") %>%
      as.matrix() %>%
      as.data.frame()
    res_f <- res_f[1,2]
    return(res_f)
  }

  score_network <- function(disease_network = disease_network, target = target, method){

    disease_network <- as.data.frame(disease_network[,1:2])
    colnames(disease_network)<- c("node1","node2")
    target <- unique(target)
    target <- intersect(target, unique(c(disease_network$node1, disease_network$node2)))

    disease_network2 <- disease_network[!disease_network$node1 %in% target,]
    disease_network2 <- disease_network2[!disease_network2$node2 %in% target,]
    g1 <- igraph::graph.data.frame(disease_network, directed = F)
    g2 <- igraph::graph.data.frame(disease_network2, directed = F)

    if (method == "degree") {
      degree <- (mean(igraph::centr_degree(g2)$res) - mean(igraph::centr_degree(g1)$res))/mean(igraph::centr_degree(g1)$res)
      return(degree)
    }

    if (method == "distance") {
      distance <- (igraph::mean_distance(g2, directed = F, unconnected = TRUE) - igraph::mean_distance(g1, directed = F, unconnected = TRUE))/igraph::mean_distance(g1, directed = F, unconnected = TRUE)
      return(distance)
    }

    if (method == "total") {
      degree <- (mean(igraph::centr_degree(g2)$res) - mean(igraph::centr_degree(g1)$res))/mean(igraph::centr_degree(g1)$res)
      distance <- (igraph::mean_distance(g2, directed = F, unconnected = TRUE) - igraph::mean_distance(g1, directed = F, unconnected = TRUE))/igraph::mean_distance(g1, directed = F, unconnected = TRUE)
      total_disturbance_rate <- distance - degree
      return(total_disturbance_rate)
    }
  }

  cat("Adjust Tanimoto \n")
  res2_list <- pbapply::pblapply(result@Fingerprint[-1], function(x){
    res_r <- replicate(n, Tanimoto_f(disease = sample(result@Fingerprint$disease), drug = x))
    res_r <- c(mean(res_r), sd(res_r))
    names(res_r) <- c("mean", "sd")
    res_r
  })

  res2_list <- as.data.frame(res2_list) %>%
    t() %>%
    as.data.frame()

  res2_list <- data.frame(Drug = names(result@Target),
                          mean= res2_list$mean,
                          sd = res2_list$sd)

  res2_list <- merge(res2_list, result@ScoreResult, by="Drug")
  res2_list$Tanimoto_adjust <- (res2_list$Tanimoto - res2_list$mean)/ res2_list$sd

  cat("Adjust Network \n")
  res3_list <- pbapply::pblapply(result@Target, function(x){
    degree <- replicate(n, score_network(disease_network = data.frame(node1 = sample(result@DiseaseNetwork[,1]),
                                                                      node2 = sample(result@DiseaseNetwork[,2]),
                                                                      stringsAsFactors = F),
                                         target = x,
                                         method = "degree"))

    distance <- replicate(n, score_network(disease_network = data.frame(node1 = sample(result@DiseaseNetwork[,1]),
                                                                        node2 = sample(result@DiseaseNetwork[,2]),
                                                                        stringsAsFactors = F
    ),
    target = x,
    method = "distance"))

    total <- replicate(n, score_network(disease_network = data.frame(node1 = sample(result@DiseaseNetwork[,1]),
                                                                     node2 = sample(result@DiseaseNetwork[,2]),
                                                                     stringsAsFactors = F),
                                        target = x,
                                        method = "total"))

    res_r2 <- c(mean(degree), sd(degree),
                mean(distance), sd(distance),
                mean(total), sd(total))
    names(res_r2) <- c("mean_degree", "sd_degree",
                       "mean_distance", "sd_distance",
                       "mean_total", "sd_total")
    res_r2
  })

  res3_list <- as.data.frame(res3_list) %>%
    t() %>%
    as.data.frame()

  res3_list$Drug <- names(result@Target)
  res3_list <- merge(res3_list, res2_list, by="Drug",sort = F)

  res3_list$Mean_degree_adjust <- (res3_list$Mean_degree_disturbance_rate - res3_list$mean_degree)/res3_list$sd_degree
  res3_list$Mean_distance_adjust <- (res3_list$Mean_distance_disturbance_rate - res3_list$mean_distance)/res3_list$sd_distance
  res3_list$Total_adjust <- (res3_list$Total_disturbance_rate - res3_list$mean_total)/res3_list$sd_total

  result@ScoreResult <- res3_list[,-c(2:9)]
  result@Adjust <- TRUE
  return(result)
}
