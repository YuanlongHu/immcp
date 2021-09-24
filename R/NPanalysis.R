#' NP analysis.
#'
#'
#' @title createNP
#' @param data BasicData object.
#' @param DisNet data frame
#' @return a data frame
#' @export
#' @author Yuanlong Hu


createNP <- function(data, DisNet){

  Disgene <- unique(c(DisNet[,1], DisNet[,2]))
  data <- data@Relationship

  x_target <- data[data$col2 == "target",]
  x_target <- x_target[x_target$to %in% Disgene,]
  Disgene <- unique(x_target$to)

  DisNet <- DisNet[,1:2]
  colnames(DisNet) <- c("from","to")
  DisNet$col1 <- "target"
  DisNet$col2 <- "target"
  DisNet <- DisNet[DisNet$from %in% Disgene,]
  DisNet <- DisNet[DisNet$to %in% Disgene,]

  if(unique( x_target$col1)=="compound"){
    compound <- unique( x_target$from)
    herb_compound <- data[data$col2 == "compound",]
    herb_compound <- herb_compound[data$to %in% compound,]
    herb <- unique(herb_compound$from)
    drug_herb <- data[data$col1 == "drug" & data$to %in% herb,]
    data <- Reduce(rbind, list(x_target[,1:4], herb_compound[,1:4], drug_herb[,1:4], DisNet))
  }

  if(unique( x_target$col1)=="herb"){
    herb <- unique(x_target$from)
    drug_herb <- data[data$col1 == "drug" & data$to %in% herb,]
    data <- Reduce(rbind, list(x_target[,1:4], drug_herb[,1:4], DisNet))
  }
  return(data)
}
