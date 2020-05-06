#' Calculate the similarity of immune fingerprint and the characteristics of network topology
#'
#'
#' @title score_immpc
#' @param disease_biomarker
#' @param disease_network A data frame of disease network.This data frame containing a symbolic edge list in the first two columns.
#' @param target A data frame of drug target.
#' @return ScoreResult object
#' @export
#' @author Yuanlong Hu


score_immpc <- function(disease_biomarker,
                        disease_network,
                        target,
                        geneset = NULL){
  # geneset list to data.frame

  if(is.null(geneset)){
    geneset0 <- genesetlist
  }

  if(class(geneset) == "list"){
    for(i in names(geneset)){
      geneset1 <- data.frame(feature = rep(i, length(geneset[i])), genesymbol = geneset[i])
      names(geneset1) <- c("feature", "genesymbol")
      geneset0 <- rbind(geneset0, geneset1)
    }
  }

  target <- target[,c(1,2)]
  names(target) <- c("c1", "c2")
  disease <- data.frame(c1 = rep("disease", length(disease_biomarker)), c2 = disease_biomarker)

  enrich_f <- function(target_character, geneset = geneset0){
    names(geneset) <- c("c1", "c2")
    enrich_drug <- clusterProfiler::enricher(target_character,
                                             TERM2GENE = geneset,
                                             minGSSize = 2,maxGSSize = Inf,
                                             pvalueCutoff = 0.05,
                                             qvalueCutoff = 0.1)
    enrich_drug <- enrich_drug@result
    enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalue<0.1,]
    fingerprint_drug <- data.frame(feature = unique(geneset$c1), var = ifelse(unique(geneset$c1) %in% enrich_drug$ID, 1, 0))
    names(fingerprint_drug) <- c("feature","drug")
    return(fingerprint_drug)
  }

  cat("Extract immune fingerprint \n")
  res_enrich0 <- enrich_f(disease$c2)
  names(res_enrich0) <- c("feature", "disease")
  for (i in unique(target$c1)) {
    res_enrich <- enrich_f(target$c2[target$c1 == i])
    names(res_enrich) <- c("feature", i)
    res_enrich0 <- merge(res_enrich0, res_enrich)
  }
  rownames(res_enrich0) <- res_enrich0$feature
  res_enrich0 <- res_enrich0[,-1]

  # target is a data frame; disease is a data frame
  Tanimoto_f <- function(f = res_enrich0){
    f <- t(f)
    res_Tanimoto <- 1 - philentropy::distance(as.matrix(f), method="jaccard", use.row.names = TRUE)
    res_Tanimoto <- reshape2::melt(res_Tanimoto,value.name = "Tanimoto")
    res_Tanimoto <- res_Tanimoto[res_Tanimoto$Var1 == "disease" & res_Tanimoto$Var2 != "disease",]
    res_Tanimoto <- res_Tanimoto[,c(2,3)]
    names(res_Tanimoto) <- c("Drug", "Tanimoto")
    return(res_Tanimoto)
  }

  score_change <- function(disease_network = disease_network, target = target){

    disease_network <- as.data.frame(disease_network[,1:2])

      colnames(disease_network)<- c("node1","node2")
      target <- as.data.frame(target[,1:2])
      colnames(target) <- c("herb","target")

      degree0 <- NULL
      mean_distance0 <- NULL
      Total_disturbance_rate0 <- NULL
      for (i in unique(target[,1])) {
        target_drug <- unique(target$target[target$herb == i])
        target_drug <- intersect(target_drug, unique(c(disease_network$node1, disease_network$node2)))
        disease_network2 <- disease_network[!disease_network$node1 %in% target_drug,]
        disease_network2 <- disease_network2[!disease_network2$node2 %in% target_drug,]
        g1 <- igraph::graph.data.frame(disease_network, directed = F)
        g2 <- igraph::graph.data.frame(disease_network2, directed = F)
        degree <- (mean(igraph::centr_degree(g2)$res) - mean(igraph::centr_degree(g1)$res))/mean(igraph::centr_degree(g1)$res)
        mean_distance <- (igraph::mean_distance(g2, directed = F, unconnected = TRUE) - igraph::mean_distance(g1, directed = F, unconnected = TRUE))/igraph::mean_distance(g1, directed = F, unconnected = TRUE)
        Total_disturbance_rate <- mean_distance - degree

        names(degree) <- i
        names(mean_distance) <- i
        names(Total_disturbance_rate) <- i

        degree0 <- c(degree0, degree)
        mean_distance0 <- c(mean_distance0, mean_distance)
        Total_disturbance_rate0 <- c(Total_disturbance_rate0, Total_disturbance_rate)
    }


    res_network <- data.frame(Drug = names(degree0),
                              Degree_disturbance_rate = degree0,
                              Mean_distance_disturbance_rate = mean_distance0,
                              Total_disturbance_rate = Total_disturbance_rate0)

    return(res_network)
  }

  cat("Calculate the Tanimoto coefficient \n")
  res_Tanimoto <- Tanimoto_f(f = res_enrich0)
  res_Tanimoto2 <- Tanimoto_f(f = res_enrich0[res_enrich0$disease == 1,])
  names(res_Tanimoto2) <- c("Drug", "Tanimoto2")
  if(!is.null(disease_network)){
    cat("Calculate the characteristics of network topology \n")
    res_network <- score_change(disease_network = disease_network, target = target)
    result <- merge(res_Tanimoto, res_network, by = "Drug")
    result <- merge(result, res_Tanimoto2, by = "Drug")
  }else{
    result <- res_Tanimoto
  }

  res_ScoreResult <- new("ScoreResult",
                            ScoreResult = as.data.frame(result),
                            Fingerprint = as.data.frame(res_enrich0),
                            DiseaseNetwork = as.data.frame(disease_network),
                            DiseaseBiomarker = as.character(disease$c2),
                            Target = as.data.frame(target),
                            Adjust = FALSE)
  return(res_ScoreResult)
}

