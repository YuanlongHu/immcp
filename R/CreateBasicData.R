#' Prepare input data.
#'
#'
#' @title CreateBasicData
#' @param drug_herb A data frame containing prescription and herb
#' @param herb_target A data.frame containing herb and target
#' @return a list
#' @export
#' @author Yuanlong Hu



CreateBasicData <- function(drug_herb, herb_target){

  colnames(drug_herb) <- c("drug", "herb")
  colnames(herb_target) <- c("herb", "target")
  drug_herb <- to_list(drug_herb)
  data <- lapply(drug_herb, function(x){
    d <- herb_target$target[herb_target$herb %in% unique(x)]
    d <- unique(d)
  })

  v <- lapply(data, function(x){
    length(x) != 0
  })

  v <- unlist(v)

  data <- data[v]

  return(data)
}


