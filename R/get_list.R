#' Create a new list from a data.frame of drug target and disease biomarker as input
#'
#'
#' @title get_list
#' @param dataframe A data frame of drug target. The first column is drug name, and the second column is drug target.
#' @return list
#' @export
#' @author Yuanlong Hu



get_list <- function(dataframe){
  dataframe <- dataframe[,c(1,2)]
  names(dataframe) <- c("terms", "gene")
  target0 <- list()
  for (i in unique(dataframe[,1])) {
    target1 <- dataframe$c2[dataframe$terms == i]
    target1 <- list(target1)
    target0 <- c(target0, target1)
  }
  names(target0) <- unique(dataframe[,1])
  return(target0)
}


#' Convert list to data.frame
#'
#'
#' @title list_to_df
#' @param list A list
#' @return data frame
#' @export
#' @author Yuanlong Hu



list_to_df <- function(list){
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
#' @param geneset A data frame
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
#' @param gmtfile gmt file
#' @return data.frame
#' @importFrom utils stack
#' @export
#' @author Yuanlong Hu


read_gmt <- function(gmtfile, input_list = FALSE){
  x <- readLines(gmtfile)
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))

  geneset <- stack(res)
  geneset <- geneset[, c("ind", "values")]
  colnames(geneset) <- c("terms", "gene")

  if (input_list) {
    geneset <- get_list(geneset)
  }else(
    return(geneset)
  )


}
