##' @rdname extrFP
##' @exportMethod extrFP

setMethod("extrFP", signature(drug_target = "BasicData"),
          function(drug_target, disease_biomarker, method = "enrich") {
            extrFP.BasicData(drug_target = drug_target, disease_biomarker=disease_biomarker, method = method)
          })



#' @rdname extrFP
#' @param disease_biomarker A character of disease biomarkers or an order ranked geneList.
#' @param drug_target A data frame or list of drug target.
#' @param method one of "enrich" and "gsea"
#' @return ScoreFP object
#' @importFrom pbapply pblapply
#' @importFrom clusterProfiler enricher
#' @importFrom clusterProfiler GSEA


extrFP.BasicData <- function(drug_target, disease_biomarker, method = "enrich"){

  Relationship <- drug_target@Relationship
  CompoundAnno <- drug_target@CompoundAnno
  drug_target <- drug_target@BasicData

  geneset <- "KEGG"
  if (geneset == "KEGG"){
    geneset0 <- genesetlist$KEGGPATHID2EXTID
    geneset0 <- geneset0[!geneset0$from %in% genesetlist$KEGGPATHID2NAME_out,]
    geneset0 <- geneset0[,-1]
  }

  if (class(geneset) == "list"){
    geneset0 <- NULL
    for(i in names(geneset)){
      geneset1 <- data.frame(feature = rep(i, length(geneset[i])), genesymbol = geneset[i])
      names(geneset1) <- c("feature", "genesymbol")
      geneset0 <- rbind(geneset0, geneset1)
    }
  }

  if (class(geneset)== "data.frame") geneset0 <- geneset

  #if (class(drug_target)=="data.frame") drug_target <- to_list(drug_target)

  enrich_f <- function(target_character, geneset = geneset0, method){
    names(geneset) <- c("c1", "c2")
    if (method=="enrich"){
       enrich_drug <- enricher(target_character,
                               TERM2GENE = geneset,
                               minGSSize = 5, maxGSSize = 500,
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.1)
       enrich_drug <- enrich_drug@result
       enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalue<0.1,]
    }
    if (method=="gsea"){
      target_character <- sort(target_character,decreasing = TRUE)
       enrich_drug <- GSEA(target_character,
                           exponent = 1,
                           minGSSize = 5, maxGSSize = 500,
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
                      Geneset = geneset,
                      Relationship = Relationship,
                      CompoundAnno = CompoundAnno
                        )
  return(res_ScoreFP1)
}


##' @param data A list of order ranked geneList.
##' @param search A vecter of order ranked geneList.
##' @param geneset A data.frame of 2 column with term and gene
##' @importFrom clusterProfiler GSEA
##' @importFrom pbapply pblapply
##' @return
##' @noRd


MoASim <- function(data, search=NULL, geneset){
  if (is.null(search)) {
    data
  }else{
    data <- c(list(search), data)
  }
  colnames(geneset) <- c("name","gene")
  geneset_temp <- unique(geneset[,1])

  geneset_temp <- data.frame(index=c(1:length(geneset_temp)), name=geneset_temp)
  geneset <- merge(geneset, geneset_temp, by = "name")
  TERM2GENE <- data.frame(term=geneset$index, gene=geneset$gene)
  TERM2NAME <- data.frame(term=geneset$index, name=geneset$name)
  extrNES <- function(genelist, TERM2GENE = TERM2GENE, TERM2NAME =TERM2NAME){
    genelist <- sort(genelist,decreasing = TRUE)
    res <- suppressMessages(GSEA(genelist,
                                 exponent = 1,
                                 minGSSize = 5, maxGSSize = 500,
                                 pvalueCutoff = 1, pAdjustMethod = "BH",
                                 TERM2GENE = TERM2GENE,
                                 TERM2NAME = TERM2NAME,
                                 verbose = TRUE, seed = FALSE,
                                 by = "fgsea"))
    res <- res@result
    res <- res[order(res$ID),]
    #enrich_drug <- enrich_drug[enrich_drug$pvalue<0.05 & enrich_drug$qvalues<0.1,]
    res2 <- res$NES
    names(res2) <- res$ID
    return(res2)
  }
  res_NES <- pblapply(data, function(x){
    suppressWarnings(extrNES(x, TERM2GENE = TERM2GENE, TERM2NAME =TERM2NAME))
  })

  return(res_NES)
}

##' Select features related to phenotype using Boruta
##'
##'
##' @title getF
##' @param expr A matrix of expression values where rows correspond to genes and columns correspond to samples.
##' @param pdata A character of phenotype.
##' @param level one of the gene or pathway
##' @param geneset A data frame of geneset containing two columns.
##' @param withTentative If set to TRUE, Tentative attributes will be also returned.
##' @return A character of features.
##' @importFrom GSVA gsva
##' @importFrom Boruta Boruta
##' @importFrom Boruta getSelectedAttributes
##' @export
##' @author Yuanlong Hu

getF <- function(expr, pdata, level = "gene", withTentative = TRUE, geneset){
  expr <- as.matrix(expr)



  if (level == "pathway"){
    if (class(geneset)=="list") geneset
    if (class(geneset) == "data.frame") geneset <- to_list(geneset)
    message("Run ssgsea ...\n")
    res1 <- gsva(expr = expr,
                 gset.idx.list = geneset,
                 method = "ssgsea",
                 abs.ranking = FALSE,
                 min.sz = 1,
                 max.sz = Inf,
                 mx.diff = TRUE,
                 ssgsea.norm = TRUE,
                 verbose = TRUE,
    )
  }

  if (level == "gene"){
    res1 <- expr
  }

  res2 <- t(res1)
  res2 <- as.data.frame(res2)
  res2$group <- pdata
  res2$group <- factor(res2$group)

  message("Run Boruta ...\n")
  set.seed(1234)
  res_B <- Boruta(group~.,data = res2,doTrace = 0)
  features <- getSelectedAttributes(res_B, withTentative = withTentative)

  message("Done...")
  return(features)
}
