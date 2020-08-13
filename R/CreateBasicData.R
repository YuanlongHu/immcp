#' Prepare input format.
#'
#'
#' @title PrepareData
#' @param data A data frame containing interaction information.
#' @param col1 A charactor containing "drug", "herb", "compound", or "target".
#' @param col2 A character containing "drug", "herb", "compound", or "target".
#' @param format one of "single" or "basket".
#' @param sep Separator.
#' @return a list
#' @importFrom magrittr %>%
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   drug_herb <- PrepareData(drugSample$drug_herb, col1 = "drug", col2 = "herb")
#'   herb_target <- PrepareData(drugSample$herb_target,
#'                              col1 = "herb", col2 = "target",
#'                              format = "basket", sep = ", ")

PrepareData <- function(data, col1, col2, format = "single", sep){

  if (! col1 %in% c("drug","herb","compound","target")){
     stop("The 'col1' must be one of 'drug', 'herb', 'compound', or 'target'")
   }

  if (! col2 %in% c("drug","herb","compound","target")){
    stop("The 'col2' must be one of 'drug', 'herb', 'compound', or 'target'")
  }

  if (! class(data) %in% c("list","data.frame")) {
    stop("Input data must be a list or data.frame")
  }

  if (class(data) == "list") {
    data <- to_df(data)
  }

  if (ncol(data)>2) {

    anno <- data.frame(from = data[,c(col1, col2, "compound_id")== "compound"],
                       id = paste0("c",data[,3]))
    data <- as.data.frame(data[,1:2])
  }else{
    anno <- data.frame(from = unique(data[,1]),
                       id = paste0("c",1:length(unique(data[,1])))
                       )
    data <- as.data.frame(data[,1:2])
  }


  if (format == "single"){
    data <- data
  }

  if (format == "basket"){
    data <- to_list(data, input = format, sep = sep) %>%
      to_df()
  }
  colnames(data) <- c("from", "to")
  data$col1 <- col1
  data$col2 <- col2

  data <- merge(data, anno, by= "from")


  return(data)
}


#' Prepare input data.
#'
#'
#' @title CreateBasicData
#' @param ... A data frame from `PrepareData`.
#' @return a list
#' @importFrom magrittr %>%
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   drug_herb <- PrepareData(drugSample$drug_herb, col1 = "drug", col2 = "herb")
#'   herb_target <- PrepareData(drugSample$herb_target,
#'                              col1 = "herb", col2 = "target",
#'                              format = "basket", sep = ", ")
#'   drug_target <- CreateBasicData(drug_herb, herb_target)


CreateBasicData <- function(...){

  da <- list(...)
  rp <- length(da)-1

  index <- lapply(da, function(x){

    index <- unique(x$col1) == "compound"
    return(index)
  })

  index <- unlist(index)

  if(length(da[index]) > 0){
    compound <- da[index][[1]]
    compound <- data.frame(compound = compound$from,
                           id = compound$id)
    compound <- compound[!duplicated(compound$compound),]
  }else{
    compound <- data.frame()
  }


  if(rp == 0){

    da <- da[[1]]
    target <- da
    relation <- unique(c(target$col1, target$col2))

  }else{
    da <- Reduce("rbind",da) %>%
      as.data.frame()

    relation <- unique(c(da$col1, da$col2))

    target <- da[da$col1 == "target"|da$col2 == "target",]
    re1 <- unique(c(target$col1, target$col2))

    re1 <- re1[re1 != "target"]
    re2 <- "target"

    for (i in 1:rp) {
      m_data <- da[da$col1 != re2 & da$col2 != re2,]
      m_data <- m_data[m_data$col1 == re1|m_data$col2 == re1,]
      m_data

      re2 <- re1
      re <- unique(c(m_data$col1, m_data$col2))
      re1 <- re[re != re1]

      target <- M_df(target, m_data) %>%
        to_df()
    }

  }

  target <- to_list(target)

  res <- new("BasicData",
             BasicData = target,
             Key = relation,
             Relationship = da,
             CompoundAnno = compound
             )
  return(res)
}

#' Prepare input data.
#'
#'
#' @title M_df
#' @param drug_herb A data frame or listcontaining prescription and herb
#' @param herb_target A data.frame or list containing herb and target
#' @return a list
#' @noRd
#' @author Yuanlong Hu

M_df <- function(herb_target, drug_herb){
  herb_target <- herb_target[,1:2]
  drug_herb <- drug_herb[,1:2]
  colnames(drug_herb) <- c("from","to")
  colnames(herb_target) <- c("from","to")

  # not matching
  w <- unique(drug_herb$to)
  w <- w[!w %in% unique(herb_target$from)]
  if (length(w)>0) {
    w <- paste0(w, collapse = ", ")
    warning(paste0("The following herb do not find the target information: \n", w))
  }

  drug_herb <- to_list(drug_herb)
  data <- lapply(drug_herb, function(x){
    d <- herb_target$to[herb_target$from %in% unique(x)]
    d <- unique(d)
    return(d)
  })

  v <- lapply(data, function(x) length(x) != 0)

  v <- unlist(v)
  data <- data[v]
  return(data)
}

