#' Adjust Score for Permutation Test
#'
#' @title score_adjust
#' @param result ScoreResult object
#' @param n integer.
#'

#setMethod("score_adjust",
#          signature(result = "ScoreResult", n = "integer"),
score_adjust <- function(result, n = 100){

  disease <- result@Fingerprint$disease
  drug <- result@Fingerprint[,-1]

  cat("Adjust Tanimoto \n")
  Tanimoto_f <- function(disease, drug){
    f <- rbind(disease, drug)
    f <- Matrix::Matrix(f)
    res_Tanimoto <- proxyC::simil(f, method = "jaccard")
    res_Tanimoto <- res_Tanimoto[1, 2]
    return(res_Tanimoto)
  }
  res_rep0 <- NULL
  pb <- tkProgressBar(title="Progress",label="Completed %", min=0, max=100, initial = 0, width = 300)
  u <- c(1:ncol(drug))
  for (i in u) {
    res_rep <- replicate(n,Tanimoto_f(disease=sample(disease), drug=drug[,i]))
    res_rep <- (result@ScoreResult$Tanimoto[i] - mean(res_rep))/sd(res_rep)
    names(res_rep) <- colnames(drug)[i]
    res_rep0 <- c(res_rep0, res_rep)

    info <- sprintf("Completed %d%%", round(i*100/length(u)))
    setTkProgressBar(pb, value = i*100/length(u), title = sprintf("Progress (%s)",info),label = info)
  }
  close(pb)
  result@ScoreResult$Tanimoto_adj <- res_rep0

  # target is character
  score_network <- function(disease_network = disease_network, target = target, method = "degree"){

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
      total_disturbance_rate <- distance - degree
      return(total_disturbance_rate)
    }
  }

  cat("Adjust Network")
  degree0 <- NULL
  distance0 <- NULL
  total_disturbance_rate0 <- NULL

  Target <- result@Target
  pb <- tkProgressBar(title="Progress",label="Completed %", min=0, max=100, initial = 0, width = 300)
  u <- c(1:length(unique(Target[,1])))
  for (i in u) {
    j <- as.character(unique(Target[,1])[i])
    degree <- replicate(n, score_network(disease_network = data.frame(node1 = sample(result@DiseaseNetwork[,1]),
                                                                      node2 = sample(result@DiseaseNetwork[,2]),
                                                                      stringsAsFactors = F
    ),
    target = as.character(Target$c2[Target$c1 == j]),
    method = "degree")
    )
    degree
    distance <- replicate(n, score_network(disease_network = data.frame(node1 = sample(result@DiseaseNetwork[,1]),
                                                                        node2 = sample(result@DiseaseNetwork[,2]),
                                                                        stringsAsFactors = F
    ),
    target = Target$c2[Target$c1 == j],
    method = "distance")
    )
    total <- replicate(n, score_network(disease_network = data.frame(node1 = sample(result@DiseaseNetwork[,1]),
                                                                     node2 = sample(result@DiseaseNetwork[,2]),
                                                                     stringsAsFactors = F
    ),
    target = Target$c2[Target$c1 == j],
    method = "total")
    )

    degree <- (result@ScoreResult$Degree_disturbance_rate[result@ScoreResult$Drug == j] - mean(degree))/sd(degree)
    distance <- (result@ScoreResult$Mean_distance_disturbance_rate[result@ScoreResult$Drug == j] - mean(distance))/sd(distance)
    total <- (result@ScoreResult$Total_disturbance_rate[result@ScoreResult$Drug == j] - mean(total))/sd(total)

    names(degree) <- j
    names(distance) <- j
    names(total) <- j

    degree0 <- c(degree0, degree)
    distance0 <- c(distance0, distance)
    total_disturbance_rate0 <- c(total_disturbance_rate0, total)

    info <- sprintf("Completed %d%%", round(i*100/length(u)))
    setTkProgressBar(pb, value = i*100/length(u), title = sprintf("Progress (%s)",info),label = info)

  }
   close(pb)
   res_network <- data.frame(Drug = names(degree0),
                             degree_adj = degree0,
                             distance_adj = distance0,
                             total_disturbance_rate_adj = total_disturbance_rate0)
   res_network
   result@ScoreResult <- merge(result@ScoreResult, res_network, by="Drug")

   return(result)
}
