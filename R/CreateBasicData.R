#' Prepare input data.
#'
#'
#' @title CreateBasicData
#' @param drug_herb A data frame or listcontaining prescription and herb
#' @param herb_target A data.frame or list containing herb and target
#' @return a list
#' @export
#' @author Yuanlong Hu

CreateBasicData <- function(drug_herb, herb_target){
  if (class(drug_herb) == "list") {
    drug_herb <- to_df(drug_herb)
  }
  if (class(herb_target) == "list") {
    herb_target <- to_df(herb_target)
  }

  drug_herb <- drug_herb[,1:2]
  herb_target <- herb_target[,1:2]

  colnames(drug_herb) <- c("drug", "herb")
  colnames(herb_target) <- c("herb", "target")

  # not matching
  w <- unique(drug_herb$herb)
  w <- w[!w %in% unique(herb_target$herb)]
  if (length(w)>0) {
    w <- paste0(w, collapse = ", ")
    warning(paste0("The following herb do not find the target information: \n", w))
  }

  drug_herb <- to_list(drug_herb)
  data <- lapply(drug_herb, function(x){
    d <- herb_target$target[herb_target$herb %in% unique(x)]
    d <- unique(d)
  })

  v <- lapply(data, function(x) length(x) != 0)

  v <- unlist(v)
  data <- data[v]

  return(data)
}


