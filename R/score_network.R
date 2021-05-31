#' Calculate the network score
#'
#'
#' @title score_network
#' @param Tar A BasicData object containing drug target.
#' @param DNet A data frame of disease network containing two columns.
#' @param n The number of times random permutation sampling.
#' @return ScoreResultNet object
#' @importFrom pbapply pblapply
#' @importFrom igraph graph.data.frame
#' @importFrom igraph centr_degree
#' @importFrom igraph mean_distance
#' @importFrom magrittr %>%
#' @importFrom stats sd
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
#'   res <- score_network(Tar = drug_target, DNet = drugSample$disease_network)

score_network <- function(Tar, DNet, n = 100){
  Relationship <- Tar@Relationship
  Tar <- Tar@BasicData

  message("----- Calculating Network Characters and Tests -----")
  res <- pbapply::pblapply(Tar, function(x){
    res <- network_char_test(DNet=DNet, target=x, n=n, samplingdata=TRUE)
    return(res)
  })

  message("----- Summarizing all results -----")

  random <- lapply(res, function(x){
    random <- x$random
    return(random)
  })
  res <- lapply(res, function(x){
    res <- x$res
    return(res)
  })
  res <- Reduce(rbind, res)
  rownames(res) <- names(Tar)
  names(random) <- names(Tar)

  message("----- Done -----")
  res_ScoreResult <- new("ScoreResultNet",
                         ScoreResult = res,
                         DiseaseNetwork = DNet,
                         Tar = Tar,
                         Relationship = Relationship,
                         adj = random) #list
  return(res_ScoreResult)
}

#' Test the change of network characters
#'
#'
#' @title network_char_test
#' @param DNet A data frame of disease network containing two columns.
#' @param target character; drug target.
#' @param n The number of times random sampling.
#' @param method the method of test; "PT":Permutation test.
#' @param samplingdata Whether to output the random data.
#' @return a data frame
#' @importFrom igraph graph.data.frame
#' @importFrom igraph delete.vertices
#' @importFrom igraph erdos.renyi.game
#' @importFrom igraph vcount
#' @importFrom igraph ecount
#' @importFrom igraph V
#' @importFrom igraph V<-
#' @importFrom purrr map2
#' @importFrom stats sd
#' @importFrom pbapply pblapply
#' @export
#' @author Yuanlong Hu

network_char_test <- function(DNet, target, n=100, method="PT", samplingdata=FALSE){

  g_DNet <- graph.data.frame(DNet[,1:2])
  g_DNet_char <- network_char_change(DNet=g_DNet,target=target)

  res_rand <- pblapply(as.list(1:n), function(x){
    g_rand <- erdos.renyi.game(n = vcount(g_DNet),
                               p.or.m = ecount(g_DNet),
                               type = 'gnm', weight = FALSE, mode = 'undirected'
                               )
    V(g_rand)$name <- V(g_DNet)$name
    g_rand_char <- network_char_change(DNet=g_rand, target=target, output_all=FALSE)
    rownames(g_rand_char) <- x
    return(g_rand_char)
  })

  res_rand <- Reduce(rbind, res_rand)
  g_DNet_char_change <- g_DNet_char[,-c(1:ncol(g_DNet_char)/3*2)]

  if(method=="PT"){
    z0 <- list()
    p0 <- list()
    for (i in 3:ncol(g_DNet_char_change)){
      x <- g_DNet_char_change[,i]
      y <- res_rand[,i]
      z <- (x-mean(y))/sd(y)
      p <- (length(y[abs(y) > abs(x)])+1)/(n+1)
      z0 <- c(z0, list(z))
      p0 <- c(p0, list(p))
    }

     names(z0) <- paste0("z_", names(g_DNet_char_change)[-c(1,2)])
     names(p0) <- paste0("p_", names(g_DNet_char_change)[-c(1,2)])
     res2 <- data.frame(g_DNet_char,
                        as.data.frame(z0),
                        as.data.frame(p0))
   }
  if(samplingdata){
    res2 <- list(res=res2, random=res_rand)
  }
  return(res2)
}

#' Calculate the change of network characters
#'
#'
#' @title network_char_change
#' @param DNet A data frame of disease network containing two columns.
#' @param target character; drug target.
#' @param output_all output all result
#' @return a data frame
#' @importFrom igraph graph.data.frame
#' @importFrom igraph delete.vertices
#' @export
#' @author Yuanlong Hu
#' @examples
#'
#'   data("drugSample")
#'   network_char_change(DNet=drugSample$disease_network,
#'                       target=unique(drugSample$herb_target[,2]))


network_char_change <- function(DNet, target, output_all=TRUE){

  if(class(DNet)=="igraph"){
    g1 <- DNet
    target <- intersect(target, V(g1)$name)
  }else{

    DNet <- as.data.frame(DNet[,1:2])
   # colnames(DNet)<- c("node1","node2")
    # overlapping gene between drug targets and disease network
    target <- intersect(target, unique(c(DNet[,1], DNet[,2])))
    g1 <- graph.data.frame(DNet, directed = F)
  }

  # Change
  g2 <- delete.vertices(g1, target)

  netchar_g1 <- network_char(g1, T)
  netchar_g2 <- network_char(g2, T)
  change <- (netchar_g2 - netchar_g1)/netchar_g1

  # Summary
  names(netchar_g1) <- paste0("G1_", names(netchar_g1))
  names(netchar_g2) <- paste0("G2_", names(netchar_g2))
  names(change) <- paste0("Change_", names(change))
  change[is.na(change)] <- 0

  if(output_all){
    res_network <- data.frame(netchar_g1, netchar_g2, change)
  }else{
    res_network <- change
  }
  return(res_network)
}

##' @rdname imm_centr
##' @exportMethod imm_centr

setMethod("imm_centr", signature(x = "data.frame"),
          function(x) {
            imm_centr.data.frame(x)
          })


##' @rdname imm_centr
##' @exportMethod imm_centr

setMethod("imm_centr", signature(x = "ScoreResultNet"),
          function(x, drug, node = "target", net = "disease") {
            imm_centr.ScoreResultNet(x, drug, node = node, net = net)
          })


#' @rdname imm_centr
#' @importFrom igraph degree
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @importFrom igraph eigen_centrality
#' @importFrom igraph V
#' @author Yuanlong Hu


imm_centr.data.frame <- function(x, node, net){
  x <- as.data.frame(x)[,1:2]
  names(x) <- c("from", "to")
  g <- graph.data.frame(x, directed = F)

  res <- list(
    degree = degree(g, v = V(g), mode = "all"),
    closeness = closeness(g, vids = V(g), mode = "all"),
    betweenness = betweenness(g, v = V(g)),
    eigen = eigen_centrality(g)$vector
  )
  res <- as.data.frame(res)
  return(res)
}

#' @rdname imm_centr
#' @param drug drug name
#' @param node Nodes that need to be evaluated. one of "disease" and "target.
#' @param net Network. one of "disease" and "target.
#' @importFrom igraph degree
#' @importFrom igraph closeness
#' @importFrom igraph betweenness
#' @importFrom igraph eigen_centrality
#' @importFrom igraph V
#' @author Yuanlong Hu

imm_centr.ScoreResultNet <- function(x, drug, node, net){

  g <- x@DiseaseNetwork
  names(g) <- c("from", "to")

  if (net == "target"){
    Target <- x@ScoreResult[drug,"Target"]
    Target <- strsplit(Target, split = ", ")[[1]]
    g <- g[g$from %in% Target & g$to %in% Target,]
  }

  g <- graph.data.frame(g, directed = F)

  res <- list(
    degree = degree(g, v = V(g), mode = "all"),
    closeness = closeness(g, vids = V(g), mode = "all"),
    betweenness = betweenness(g, v = V(g)),
    eigen = eigen_centrality(g)$vector
  )

  res <- as.data.frame(res)
  if (node == "target"){
    Target <- x@ScoreResult[drug,"Target"]
    Target <- strsplit(Target, split = ", ")[[1]]
    res <- res[Target,]

  }
  res <- res[order(res$degree, decreasing = T),]
  return(res)
}
