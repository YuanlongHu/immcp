#' Kolmogorov-Smirnov tests for node characters between networks
#'
#'
#' @title network_node_ks
#' @param graph1 A igraph object.
#' @param graph2 A igraph object.
#' @param replicate Number vector, the number of conduct bootstrapping sampling replications.
#' @return A data frame
#' @importFrom pbapply pblapply
#' @importFrom stats ks.test
#' @export
#' @author Yuanlong Hu

network_node_ks <- function(graph1, graph2, replicate=1000){


  graph <- list(graph1, graph2)
  network_char <- lapply(graph, function(x){
    network_char <- suppressWarnings(network_char(x))
    return(network_char)
  })


  # bootstrap
  message(">>>>> Bootstrapping Sampling <<<<<")
  network_char_boot <- lapply(network_char,
                              function(x){
                                set.seed(12345)
                                node_df <- lapply(x, function(y) replicate(replicate, sample(y, 1, replace = TRUE)))
                                node_df <- as.data.frame(node_df)
                                names(node_df) <- c("boot_degree", "boot_closeness",
                                                    "boot_betweenness", "boot_eigenvector",
                                                    "boot_transitivity")
                                return(node_df)
                              })

  # Kolmogorov-Smirnov tests
  message(">>>>> Kolmogorov-Smirnov Tests <<<<<")
  result <- list()
  n <- length(network_char_boot) # a list

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {

      ks <- lapply(as.list(c("boot_degree", "boot_closeness",
                             "boot_betweenness", "boot_eigenvector",
                             "boot_transitivity")),
                   function(x){
                     node_i <- network_char_boot[[i]][ ,x]
                     node_j <- network_char_boot[[j]][ ,x]
                     ks <- suppressWarnings(ks.test(node_i, node_j, alternative = 'two.sided'))
                     return(ks$p.value)
                   })

      # Summary p-value
      result <- c(result,
                  list(
                    c("Network1" = names(network_char_boot)[i],
                      "Network2" = names(network_char_boot)[j],
                      "ks_degree" = ks[[1]],
                      "ks_closeness" = ks[[2]],
                      "ks_betweenness" = ks[[3]],
                      "ks_eigenvector" = ks[[4]],
                      "ks_transitivity" = ks[[5]])
                  )
      )
    }
  }

  message(">>>>> Summary Result <<<<<")

  result <- result[[1]]
  return(result)
}
