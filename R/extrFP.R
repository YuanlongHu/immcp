#' Calculate the Drug Fingerprints
#'
#'
#' @title extrFP
#' @param disease_biomarker A character.
#' @param drug_target A data frame of drug target.
#' @param geneset A list.
#' @return ScoreResult object
#' @importFrom clusterProfiler enricher
#' @importFrom pbapply pblapply
#' @export
#' @author Yuanlong Hu


extrFP <- function(disease_biomarker, drug_target, geneset){

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


  output_list <- function(disease_biomarker, drug_target){
    drug_target <- drug_target[,c(1,2)]
    names(drug_target) <- c("c1", "c2")
    target0 <- list(disease_biomarker)
    for (i in unique(drug_target[,1])) {
      target1 <- drug_target$c2[drug_target$c1 == i]
      target1 <- list(target1)
      target0 <- c(target0, target1)
    }
    names(target0) <- c("disease", unique(target[,1]))
    return(target0)
  }

  enrich_f <- function(target_character, geneset = geneset0){
    names(geneset) <- c("c1", "c2")
    enrich_drug <- enricher(target_character,
                                             TERM2GENE = geneset,
                                             minGSSize = 2,maxGSSize = Inf,
                                             pvalueCutoff = 0.05,
                                             qvalueCutoff = 0.1)
    enrich_drug <- enrich_drug@result
    enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalue<0.1,]

    fingerprint_drug <- ifelse(unique(geneset$c1) %in% enrich_drug$ID, 1, 0)
    names(fingerprint_drug) <- unique(geneset$c1)
    return(fingerprint_drug)
  }

  target <- output_list(disease_biomarker, drug_target)
  f <- pblapply(target, function(x){
    enrich_f(x, geneset = geneset0)
  })

  res_ScoreResult <- new("ScoreResult",
                         ScoreResult = NULL,
                         Fingerprint = f,
                         DiseaseNetwork = as.data.frame(disease_network),
                         DiseaseBiomarker = as.character(target[[1]]),
                         Target = target[-1],
                         Adjust = NULL)
  return(res_ScoreResult)
}
