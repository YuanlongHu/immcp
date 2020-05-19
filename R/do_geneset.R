#' Create a new list from a data.frame of drug target and disease biomarker as input
#'
#'
#' @title to_list
#' @param dataframe a data frame of 2 column with term/drug and gene
#' @return list
#' @export
#' @author Yuanlong Hu


to_list <- function(dataframe){
  dataframe <- dataframe[,c(1,2)]
  names(dataframe) <- c("terms", "gene")
  target0 <- list()
  for (i in unique(dataframe[,1])) {
    target1 <- dataframe$gene[dataframe$terms == i]
    target1 <- list(target1)
    target0 <- c(target0, target1)
  }
  names(target0) <- unique(dataframe[,1])
  return(target0)
}


#' Convert list to data.frame
#'
#'
#' @title to_df
#' @param list a list containing gene sets
#' @return data frame
#' @export
#' @author Yuanlong Hu



to_df <- function(list){
  list0 <- NULL
  for (i in 1:length(list)) {
    list1 <- data.frame(term = rep(names(list)[i]),
                        gene = list[[i]])

    list0 <- rbind(list0, list1)
  }
  return(list0)
}


#' prints data frame to a gmt file
#'
#'
#' @title write_gmt
#' @param geneset A data.frame of 2 column with term/drug and gene
#' @param gmt_file A character of gmt file name.
#' @return gmt file
#' @export
#' @author Yuanlong Hu


write_gmt <- function(geneset, gmt_file){

  geneset <<- tapply(geneset[,1],as.factor(geneset[,2]),function(x) x)
  sink(gmt_file)
  for (i in 1:length(geneset)){
    cat(names(geneset)[i])
    cat('\tNA\t')
    cat(paste(geneset[[i]], collapse = '\t'))
    cat('\n')
  }
  sink()
}


#' parse gmt file to a data.frame
#'
#'
#' @title write_gmt
#' @param gmtfile A GMT file name or URL containing gene sets.
#' @param out_type A character vector of object name. one of "data.frame", "list", "GeneSetCollection"
#' @return data.frame, list or GeneSetCollection
#' @importFrom utils stack
#' @importFrom GSEABase getGmt
#' @importFrom clusterProfiler read.gmt
#' @export
#' @author Yuanlong Hu


read_gmt <- function(gmtfile, out_type = "data.frame"){

 if (out_type == "list") geneset <- to_list(read.gmt(gmtfile))
 if (out_type == "data.frame") geneset <- read.gmt(gmtfile)
 if (out_type == "GeneSetCollection") geneset <- getGmt(gmtfile)
 return(geneset)
}
