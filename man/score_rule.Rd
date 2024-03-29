% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_rule.R
\name{score_rule}
\alias{score_rule}
\title{Mining herb-herb associations with Apriori}
\usage{
score_rule(BasicData, drug = NULL, support = 0.1, confidence = 0.8, lift = 1)
}
\arguments{
\item{BasicData}{BasicData object.}

\item{drug}{Charactor vector of drug names to analyze, default to \code{NULL}.}

\item{support}{A numeric value for the minimal support of an item set, default to 0.1.}

\item{confidence}{A numeric value for the minimal confidence of an item set, default to 0.8.}

\item{lift}{A numeric value for the minimal lift of an item set, default to 1.}
}
\value{
A HerbResult object.
}
\description{
Mine herb-herb association rules of prescription using the Apriori algorithm.
}
\examples{
\dontrun{
data(drugdemo)
drug_herb <- PrepareData(drugdemo$drug_herb, from = "drug", to="herb")
herb_compound <- PrepareData(drugdemo$herb_compound, from = "herb", to="compound")
compound_target <- PrepareData(drugdemo$compound_target, from = "compound", to="target")
disease <- PrepareData(drugdemo$disease, diseaseID = "disease",from = "target", to="target")
BasicData <- CreateBasicData(drug_herb, herb_compound, compound_target, diseasenet = disease)
res <- score_rule(BasicData, support = 0.1,confidence = 0.8, lift = 1)
}
}
\author{
Yuanlong Hu
}
