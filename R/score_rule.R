#' Mine herb-herb association rules of prescription using the Apriori algorithm.
#'
#'
#' @title Mining herb-herb associations with Apriori
#' @rdname score_rule
#' @param BasicData BasicData object.
#' @param drug Charactor vector of drug names to analyze, default to `NULL`.
#' @param support A numeric value for the minimal support of an item set, default to 0.1.
#' @param confidence A numeric value for the minimal confidence of an item set, default to 0.8.
#' @return A HerbResult object.
#' @importFrom dplyr %>%
#' @importFrom arules apriori
#' @importFrom arules inspect
#' @importFrom arules itemFrequency
#' @importFrom methods as
#' @importFrom igraph graph.data.frame
#' @export
#' @author Yuanlong Hu
#' \dontrun{
#' data(drugdemo)
#' drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
#' herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
#' compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
#' disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")
#' BasicData <- CreateBasicData(drug_herb, herb_compound, compound_target, diseasenet = disease)
#' res <- score_rule(BasicData, support = 0.1,confidence = 0.8)
#' }

score_rule <- function(BasicData, drug = NULL, support = 0.1,confidence = 0.8){

  vertices <- BasicData@vertices

  druglist <- lapply(as.list(vertices$name[vertices$type == "drug"]), function(x){
    x <- subset_network(BasicData, from = x, to = vertices$name[vertices$type == "herb"])
    as_data_frame(x@drugnet)
  })
  druglist <- Reduce(rbind, druglist)
  if(!is.null(drug)) druglist <- druglist[druglist$from %in% drug,]

  parameter <- list(minlen = 2,
                    maxlen = 2,
                    support = support,
                    confidence = confidence)
  data <- as(split(druglist$to, druglist$from),"transactions")
  rules <-  suppressMessages(apriori(data,parameter=parameter) %>%
    inspect())
  rules <- rules[rules$lift > 1,-2]
  rules$lhs <- gsub("[{*}]", "", rules$lhs)
  rules$rhs <- gsub("[{*}]", "", rules$rhs)
  names(rules)[1] <- "from"
  names(rules)[2] <- "to"

  herb <- data.frame(name=names(itemFrequency(data, type = "relative")),
                     frequency_relative = itemFrequency(data, type = "relative"),
                     frequency_absolute = itemFrequency(data, type = "absolute"),
                     type = "herb")

  g <- graph.data.frame(rules, directed = TRUE, vertices = herb)

  HerbResult <- new("HerbResult",
                    Drug_Herb = druglist,
                    Herb_Herb = g)
  return(HerbResult)
}



