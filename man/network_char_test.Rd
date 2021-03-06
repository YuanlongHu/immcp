% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_network.R
\name{network_char_test}
\alias{network_char_test}
\title{network_char_test}
\usage{
network_char_test(DNet, target, n = 100, method = "PT", samplingdata = FALSE)
}
\arguments{
\item{DNet}{A data frame of disease network containing two columns.}

\item{target}{character; drug target.}

\item{n}{The number of times random sampling.}

\item{method}{the method of test; "PT":Permutation test.}

\item{samplingdata}{Whether to output the random data.}
}
\value{
a data frame
}
\description{
Test the change of network characters
}
\author{
Yuanlong Hu
}
