#' Create a new list from a data.frame of drug target and disease biomarker as input
#'
#'
#' @title to_list
#' @param dataframe a data frame of 2 column with term/drug and gene
#' @param input one of the single or basket
#' @param sep When 'input' is 'basket'.
#' @return list
#' @author Yuanlong Hu
#' @export
#' @examples
#' \dontrun{
#'   to_list(dataframe)
#' }

to_list <- function(dataframe, input="single", sep = ", "){
  dataframe <- dataframe[,c(1,2)]
  names(dataframe) <- c("terms", "gene")

  target0 <- list()
  if (input == "single"){
    for (i in unique(dataframe[,1])) {
      target1 <- dataframe$gene[dataframe$terms == i]
      target1 <- list(target1)
      target0 <- c(target0, target1)
    }
    names(target0) <- unique(dataframe[,1])
  }

  if (input == "basket"){
    for(i in 1:nrow(dataframe)){
      A <- dataframe[i,2]
      A <- strsplit(A,split = sep)[[1]]
      A <- list(A)
      names(A) <- dataframe[i,1]
      target0 <- c(target0, A)
    }
  }
  return(target0)
}


#' Convert list to data.frame
#'
#'
#' @title to_df
#' @param list A list containing gene sets.
#' @return A data frame.
#' @author Yuanlong Hu
#' @export
#' @examples
#' \dontrun{
#'   to_df(list)
#' }

to_df <- function(list){
  list0 <- NULL
  for (i in 1:length(list)) {
    list1 <- data.frame(term = rep(names(list)[i]),
                        gene = list[[i]])

    list0 <- rbind(list0, list1)
  }
  return(list0)
}

#' Convert BioDescr object to a list of adjacency matrix
#'
#'
#' @title to_biodescr
#' @param BioDescr A BioDescr object.
#' @return A list.
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph as_data_frame
#' @author Yuanlong Hu
#' @export
#' @examples
#' \dontrun{
#'   to_biodescr(BioDescr)
#' }

to_biodescr <- function(BioDescr){

  bd <- as_adjacency_matrix(BioDescr@drug_geneset,
                            type = "both",
                            sparse = FALSE)
  v <- as_data_frame(BioDescr@drug_geneset, "vertices")
  # bd_disease <- bd[v$name[v$type=="disease"], v$name[v$type=="pathway"]]
  # bd_drug <- bd[v$name[v$type=="drug"], v$name[v$type=="pathway"]]

  bd <- bd[v$name[v$type %in% c("disease", "drug")], v$name[v$type=="pathway"]]
  bd <- list(bd = bd,
             bd_type = v[v$type %in% c("disease", "drug"),])
  return(bd)
}

#' prints data frame to a gmt file
#'
#'
#' @title write_gmt
#' @param geneset A data.frame of 2 column with term/drug and gene.
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
#' @param out_dataframe TRUE or FALSE
#' @return data.frame, list
#' @importFrom clusterProfiler read.gmt
#' @export
#' @author Yuanlong Hu


read_gmt <- function(gmtfile, out_dataframe=TRUE){

 if (out_dataframe){
   geneset <- read.gmt(gmtfile)
 }else{
   geneset <- to_list(read.gmt(gmtfile))
 }
 return(geneset)
}
