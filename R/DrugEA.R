#' Enrichment analysis for drug target
#'
#'
#' @title DrugEA
#' @param gene 	a vector of gene symbol
#' @param internal_geneset one of "ImmGenTop150", "L1000_DN", "L1000_UP", "CMAP_DN", "CMAP_UP"
#' @param input_geneset input_geneset is a data frame or list, when internal_geneset is NULL.
#' @param output_df whether to return to the data frame
#' @param pvalueCutoff 	pvalue cutoff on enrichment tests to report
#' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param qvalueCutoff qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @return a data frame or enrichResult object
#' @importFrom clusterProfiler enricher
#' @export
#' @author Yuanlong Hu



DrugEA <- function(gene, internal_geneset = "L1000_UP",input_geneset, output_df = TRUE,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2){

  if (is.null(internal_geneset) & is.data.frame(input_geneset)){
    geneset <- input_geneset
  }

  if (is.null(internal_geneset) & is.list(input_geneset)){
    geneset <- list_to_df(input_geneset)
  }

  if (!is.null(internal_geneset)){
    geneset <- genesetlist[[internal_geneset]]
  }

  res_enrich <- enricher(
                         gene,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                        # universe,
                        minGSSize = 10, maxGSSize = Inf,
                        qvalueCutoff = qvalueCutoff,
                        TERM2GENE = geneset
  )

  if (output_df = T){
    res_enrich <- res_enrich@result
    return(res_enrich)
  }else{
    return(res_enrich)
  }

}
