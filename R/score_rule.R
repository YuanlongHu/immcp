#' Plot Disease-Drug Network
#'
#'
#' @title score_rule
#' @param BasicData BasicData object.
#' @param support cutoff of support
#' @param confidence cutoff of confidence
#' @param maxlen maxlen
#' @return Returns NULL, invisibly.
#' @importFrom dplyr %>%
#' @importFrom arules apriori
#' @importFrom arules inspect
#' @importFrom arules itemFrequency
#' @importFrom methods as
#' @importFrom igraph induced_subgraph
#' @importFrom igraph plot.igraph
#' @importFrom igraph as.undirected
#' @export
#' @author Yuanlong Hu

score_rule <- function(BasicData, support = 0.1,confidence = 0.8, maxlen = 2){

  vertices <- BasicData@vertices

  druglist <- lapply(as.list(vertices$name[vertices$type == "drug"]), function(x){
    x <- subset_network(BasicData, from = x, to = vertices$name[vertices$type == "herb"])
    as_data_frame(x@drugnet)
  })
  druglist <- Reduce(rbind, druglist)
  parameter <- list(minlen = 2,
                    support = support,
                    confidence = confidence,
                    maxlen= maxlen)
  data <- as(split(druglist$to, druglist$from),"transactions")
  rules <- data %>%
   apriori(parameter=parameter) %>%
    inspect()
  rules <- rules[rules$lift > 1,]

  herb <- data.frame(herb=names(itemFrequency(data, type="relative")),
                     frequency_relative = itemFrequency(data, type="relative"),
                     frequency_absolute = itemFrequency(data, type="absolute"))

  HerbResult <- new("HerbResult",
                    Drug_Herb = druglist,
                    Rule = rules,
                    Herb = herb)
  return(rules)
}
