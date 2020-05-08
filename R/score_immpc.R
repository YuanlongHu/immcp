#' Calculate the similarity of immune fingerprint and the characteristics of network topology
#'
#'
#' @title score_immpc
#' @param disease_biomarker A character.
#' @param disease_network A data frame of disease network.This data frame containing a symbolic edge list in the first two columns.
#' @param target A data frame of drug target.
#' @param geneset A list.
#' @return ScoreResult object
#' @importFrom igraph graph.data.frame
#' @importFrom igraph centr_degree
#' @importFrom igraph mean_distance
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#' @export
#' @author Yuanlong Hu


score_immpc <- function(FP,
                         disease_network,
                         target){

  f <- FP@Fingerprint
  target <- FP@Target
  # geneset list to data.frame

  #if(is.null(geneset)){
  #  geneset0 <- genesetlist
  #}

  #if(class(geneset) == "list"){
  #  for(i in names(geneset)){
  #    geneset1 <- data.frame(feature = rep(i, length(geneset[i])), genesymbol = geneset[i])
  #    names(geneset1) <- c("feature", "genesymbol")
  #    geneset0 <- rbind(geneset0, geneset1)
  # }
  #}

  #  output_list <- function(disease_biomarker, target){
  #    target <- target[,c(1,2)]
  #    names(target) <- c("c1", "c2")
  #    target0 <- list(disease_biomarker)
  #   for (i in unique(target[,1])) {
  #    target1 <- target$c2[target$c1 == i]
  #    target1 <- list(target1)
  #    target0 <- c(target0, target1)
  #  }
  # names(target0) <- c("disease", unique(target[,1]))
  #  return(target0)
  #}

#  enrich_f <- function(target_character, geneset = geneset0){
#    names(geneset) <- c("c1", "c2")
  #    enrich_drug <- clusterProfiler::enricher(target_character,
  #                                             TERM2GENE = geneset,
                                             #                                             minGSSize = 2,maxGSSize = Inf,
  #                                             pvalueCutoff = 0.05,
  #                                            qvalueCutoff = 0.1)
  #   enrich_drug <- enrich_drug@result
  #   enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalue<0.1,]

  #   fingerprint_drug <- ifelse(unique(geneset$c1) %in% enrich_drug$ID, 1, 0)
  #  names(fingerprint_drug) <- unique(geneset$c1)
  #  return(fingerprint_drug)
  # }

  score_network <- function(disease_network = disease_network, target = target){

    disease_network <- as.data.frame(disease_network[,1:2])
    colnames(disease_network)<- c("node1","node2")

    target <- intersect(target, unique(c(disease_network$node1, disease_network$node2)))

    disease_network2 <- disease_network[!disease_network$node1 %in% target,]
    disease_network2 <- disease_network2[!disease_network2$node2 %in% target,]

    g1 <- igraph::graph.data.frame(disease_network, directed = F)
    g2 <- igraph::graph.data.frame(disease_network2, directed = F)

    degree <- (mean(igraph::centr_degree(g2)$res) - mean(igraph::centr_degree(g1)$res))/mean(igraph::centr_degree(g1)$res)
    mean_distance <- (igraph::mean_distance(g2, directed = F, unconnected = TRUE) - igraph::mean_distance(g1, directed = F, unconnected = TRUE))/igraph::mean_distance(g1, directed = F, unconnected = TRUE)
    Total_disturbance_rate <- mean_distance - degree

    res_network <- c(degree, mean_distance, Total_disturbance_rate)

    names(res_network) <- c("Mean_degree_disturbance_rate",
                            "Mean_distance_disturbance_rate",
                            "Total_disturbance_rate")
    return(res_network)
  }

  #  cat("Extract immune fingerprint \n")

  #  target <- output_list(disease_biomarker, target)
  #  f <- pbapply::pblapply(target, function(x){
    #    enrich_f(x, geneset = geneset0)
#  })

  cat("Calculate the Tanimoto coefficient \n")
  f1 <- as.data.frame(f) %>%
    t() %>%
    Matrix::Matrix() %>%
    proxyC::simil(method = "jaccard") %>%
    as.matrix() %>%
    as.data.frame()
  f1 <- data.frame(Drug = names(f[-1]), Tanimoto = f1[-1,1], stringsAsFactors = F)

  cat("Calculate the characteristics of network topology \n")
  net1 <- pbapply::pblapply(target, function(x){
    score_network(disease_network = disease_network, target = x)
  })

  net2 <- as.data.frame(net1) %>%
    t() %>%
    as.data.frame()

  net2 <- data.frame(Drug = names(target),
                     Mean_degree_disturbance_rate = net2$Mean_degree_disturbance_rate,
                     Mean_distance_disturbance_rate = net2$Mean_distance_disturbance_rate,
                     Total_disturbance_rate = net2$Total_disturbance_rate)

  result <- merge(f1, net2, by = "Drug")

#  res_ScoreResult <- new("ScoreResult",
 #                        ScoreResult = as.data.frame(result),
 #                        Fingerprint = f,
 #                        DiseaseNetwork = as.data.frame(disease_network),
 #                        DiseaseBiomarker = as.character(target[[1]]),
  #                       Target = target[-1],
  #                       Adjust = FALSE)

  FP@ScoreResult <- as.data.frame(result)
  FP@Adjust <- FALSE
  return(FP)
}

