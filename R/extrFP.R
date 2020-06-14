#' Calculate the pathway fingerprints
#'
#'
#' @title extrFP
#' @param disease_biomarker A character of disease biomarkers or an order ranked geneList.
#' @param drug_target A data frame or list of drug target.
#' @param method one of "enrich" and "gsea"
#' @param geneset one of "ImmGenTop150" and "KEGG"
#' @return ScoreFP object
#' @importFrom pbapply pblapply
#' @importFrom clusterProfiler enricher
#' @importFrom clusterProfiler GSEA
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   FP <- extrFP(disease_biomarker = drugSample$disease_biomarker,
#'                drug_target = drugSample$herb_target,
#'                method = "enrich",
#'                geneset = "ImmGenTop150")


extrFP <- function(disease_biomarker, drug_target, method = "enrich", geneset = "ImmGenTop150"){

  if (geneset == "ImmGenTop150"){
    geneset0 <- genesetlist$ImmGenTop150
  }

  if (geneset == "KEGG"){
    geneset0 <- genesetlist$KEGGPATHID2EXTID
    geneset0 <- geneset0[!geneset0$from %in% genesetlist$KEGGPATHID2NAME_out,]
    geneset0 <- geneset0[,-1]
  }

  if (class(geneset) == "list"){
    for(i in names(geneset)){
      geneset1 <- data.frame(feature = rep(i, length(geneset[i])), genesymbol = geneset[i])
      names(geneset1) <- c("feature", "genesymbol")
      geneset0 <- rbind(geneset0, geneset1)
    }
  }

  if (class(geneset)== "data.frame") geneset0 <- geneset

  if (class(drug_target)=="data.frame") drug_target <- to_list(drug_target)

  enrich_f <- function(target_character, geneset = geneset0, method){
    names(geneset) <- c("c1", "c2")
    if (method=="enrich"){
       enrich_drug <- enricher(target_character,
                               TERM2GENE = geneset,
                               minGSSize = 2, maxGSSize = Inf,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.1)
       enrich_drug <- enrich_drug@result
       enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalue<0.1,]
    }
    if (method=="gsea"){
      target_character <- sort(target_character,decreasing = TRUE)
       enrich_drug <- GSEA(target_character,
                           exponent = 1,
                           minGSSize = 1, maxGSSize = Inf,
                           pvalueCutoff = 0.05, pAdjustMethod = "BH",
                           TERM2GENE = geneset,
                           verbose = FALSE, seed = FALSE,
                           by = "fgsea")
       enrich_drug <- enrich_drug@result
       enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalues<0.1,]
    }


    fingerprint_drug <- ifelse(unique(geneset$c1) %in% enrich_drug$ID, 1, 0)
    names(fingerprint_drug) <- unique(geneset$c1)
    return(fingerprint_drug)
  }

  message("Calculating pathway fingerprints of drug...")
  f <- pblapply(drug_target, function(x){
    enrich_f(x, geneset = geneset0, method = "enrich")
  })

  message("Calculating pathway fingerprints of disease...")
  f_disease <- suppressWarnings(enrich_f(disease_biomarker, geneset = geneset0, method))
  f_disease <- list(f_disease)
  names(f_disease) <- "disease"
  f <- c(f_disease, f)

  message("Done...")

  ifelse(class(geneset)=="character",geneset, "CusDef")
  res_ScoreFP1 <- new("ScoreFP1",
                      Fingerprint = f,
                      DiseaseBiomarker = disease_biomarker,
                      DrugTarget = drug_target,
                      FPType = "enrich",
                      Geneset = geneset
                        )
  return(res_ScoreFP1)
}
