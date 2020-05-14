#' Select Top gene as biomarker
#'
#'
#' @title get_DEG
#' @param data data frame of expression values where rows correspond to genes and columns correspond to samples.
#' @param pdata character vector of phenotype
#' @param contrasts character vector specifying contrasts
#' @return list object
#' @importFrom magrittr %>%
#' @importFrom stats model.matrix
#' @importFrom limma makeContrasts
#' @importFrom limma lmFit
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @export
#' @author Yuanlong Hu


get_DEG <- function(data, pdata, contrasts){
  Group <- factor(pdata)
  design <- model.matrix(~0 + Group)
  colnames(design) <- unique(pdata)
  contrast_matrix <- makeContrasts(contrasts,
                                   levels = design)
  fit <- lmFit(data, design) %>%
    contrasts.fit(contrast_matrix) %>%
    eBayes(fit)

  res_DEG0 <- NULL
  for (i in 1:length(contrasts)) {
    res_DEG <-  topTable(fit, adjust.method = "fdr", number = Inf, coef = i) %>%
      list()
    names(res_DEG) <- contrasts[i]
    res_DEG0 <- c(res_DEG0, res_DEG)
  }
  return(res_DEG0)
}
