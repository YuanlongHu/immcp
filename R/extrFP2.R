



extrFP2 <- function(disease_biomarker, drug_target, geneset){

  if(is.null(geneset)){
    geneset0 <- genesetlist
  }

  if(class(geneset) == "list"){
    for(i in names(geneset)){
      geneset1 <- data.frame(feature = rep(i, length(geneset[i])), genesymbol = geneset[i])
      names(geneset1) <- c("feature", "genesymbol")
      geneset0 <- rbind(geneset0, geneset1)
    }
  }


}
