#' Calculating drug-disease similarity based on biological descriptors
#'
#'
#' @title score_sim
#' @param BioDescr BioDescr object.
#' @param method method, default "jaccard".
#' @param n number
#' @return a list
#' @importFrom dplyr %>%
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#' @importFrom dplyr select
#' @importFrom purrr map2
#' @importFrom igraph as_adjacency_matrix
#' @importFrom pbapply pblapply
#' @importFrom stats sd
#' @export
#' @author Yuanlong Hu
#' @examples
#' data(drugdemo)
#' drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
#' herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
#' compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
#' disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")
#' BasicData <- CreateBasicData(drug_herb, herb_compound, compound_target, diseasenet = disease)
#' biodescr <- extr_biodescr(BasicData, geneset= "kegg")
#' res <- score_sim(biodescr, method="jaccard", n=100)
score_sim <- function(BioDescr, method="jaccard", n=1000){

  bd <- as_adjacency_matrix(BioDescr@drug_geneset,
                            type = "both",
                            sparse = FALSE)
  v <- as_data_frame(BioDescr@drug_geneset, "vertices")
  bd_disease <- bd[v$name[v$type=="disease"], v$name[v$type=="pathway"]]
  bd_drug <- bd[v$name[v$type=="drug"], v$name[v$type=="pathway"]]

  if(!is.matrix(bd_disease)){
    bd_disease <- as.matrix(bd_disease) %>% t()
    rownames(bd_disease) <- v$name[v$type=="disease"]
  }
  if(!is.matrix(bd_drug)){
    bd_drug <- as.matrix(bd_drug) %>% t()
    rownames(bd_drug) <- v$name[v$type=="drug"]
  }

  res1 <- simil(x=Matrix(bd_drug, sparse = TRUE),
                        y=Matrix(bd_disease, sparse = TRUE),
                  method = method) %>%
    as.matrix() %>%
    as.data.frame()

  res2 <- lapply(t(bd_disease) %>% as.data.frame(),function(x){
             test_simil(x,bd_drug,n=n)
  })

  res1_list <- lapply(res1, function(x){
             names(x) <- rownames(res1)
             as.list(x)
  })

  diseaseID_list <- v$name[v$type=="disease"] %>% as.list()
  names(diseaseID_list) <- v$name[v$type=="disease"]

  res_summary <- lapply(diseaseID_list, function(x){
    z_score <- map2(res1_list[[x]],res2[[x]],
                           ~((.x-mean(.y))/sd(.y))
                           )
    p_value <- map2(res1_list[[x]],res2[[x]],
                ~((length(.y[abs(.y)>abs(.x)])+1)/(n+1))
                )

    res <- data.frame(score = res1_list[[x]][[1]],
                      z_score=unlist(z_score),
                      p_value=unlist(p_value))
    return(res)

  })

  return(res_summary)


}


#' sample bd
#'
#'
#' @title sample_bd
#' @param a a
#' @param b b
#' @return a list
#' @importFrom dplyr %>%
#' @importFrom Matrix Matrix
#' @importFrom proxyC simil
#' @noRd
#' @author Yuanlong Hu


sample_bd <- function(a, b, method="jaccard"){

  a_b <- c(a, b)
  n <- length(a_b)
  n1 <- length(a)+1

  s_ab <- sample(a_b, n)

  a <- s_ab[1:length(a)]
  b <- s_ab[n1: n]

  res <- data.frame(a=a, b=b) %>%
    t() %>%
    Matrix::Matrix(sparse = TRUE) %>%
    proxyC::simil(method = method) %>%
    as.matrix()
  res <- res["a","b"]
  return(res)
}

#' test_simil
#'
#'
#' @title test_simil
#' @param bd_disease bd_disease
#' @param bd_drug bd_drug
#' @return a list
#' @importFrom dplyr %>%
#' @importFrom pbapply pblapply
#' @noRd
#' @author Yuanlong Hu
test_simil <- function(bd_disease,bd_drug, n=1000){

  bd_drug <- t(bd_drug) %>% as.data.frame()
  res <- pblapply(bd_drug, function(x){
    replicate(n, sample_bd(bd_disease, x))
  })
  return(res)
}
