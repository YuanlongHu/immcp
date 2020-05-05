#'
#'
#'
#'
#'
#'


score_adj <- function(result){
  disease <- result@Fingerprint$disease
  drug <- result@Fingerprint[,-1]

  Tanimoto_f <- function(disease, drug){
    f <- rbind(disease, drug)
    res_Tanimoto <- 1 - philentropy::distance(as.matrix(f), method="jaccard")
    return(res_Tanimoto)
  }
  res_rep0 <- NULL
  for (i in 1:ncol(drug)) {
    res_rep <- replicate(1000,Tanimoto_f(disease=sample(disease), drug=drug[,i]))
    res_rep <- (result@ScoreResult$Tanimoto[i] - mean(res_rep))/sd(res_rep)
    names(res_rep) <- colnames(drug)[i]
    res_rep0 <- c(res_rep0, res_rep)
  }
  result@ScoreResult$Tanimoto_z_score <- res_rep0


}

