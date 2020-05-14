#' Calculate the Drug Fingerprints
#'
#'
#' @title extrFP
#' @param disease_biomarker A character.
#' @param drug_target A data frame of drug target.
#' @param geneset A list.
#' @return ScoreResult object
#' @importFrom pbapply pblapply
#' @importFrom clusterProfiler enricher
#' @export
#' @author Yuanlong Hu


extrFP <- function(disease_biomarker, drug_target, geneset = NULL){

  if(is.null(geneset)){
    geneset0 <- genesetlist$ImmGenTop150
  }

  if(class(geneset) == "list"){
    for(i in names(geneset)){
      geneset1 <- data.frame(feature = rep(i, length(geneset[i])), genesymbol = geneset[i])
      names(geneset1) <- c("feature", "genesymbol")
      geneset0 <- rbind(geneset0, geneset1)
    }
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

  target <- c(disease=list(disease_biomarker), drug_target)
  f <- pblapply(target, function(x){
    enrich_f(x, geneset = geneset0)
  })

  res_ScoreFP1 <- new("ScoreFP1",
                      Fingerprint = f,
                      DiseaseBiomarker = as.character(target[[1]]),
                      DrugTarget = target[-1],
                      FPType = "enrich"
                        )
  return(res_ScoreFP1)
}
