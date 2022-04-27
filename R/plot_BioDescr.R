#' Plot Biological descriptor
#'
#' @rdname plot_BioDescr
#' @title Plot Biological descriptor
#' @param BioDescr BioDescr object.
#' @param type one of "heatmap" and "clusterplot".
#' @param cluster_k Number vector, number of cluster.
#' @param colors vector of colors.
#' @return Returns NULL, invisibly.
#' @importFrom factoextra fviz_dend
#' @importFrom dplyr %>%
#' @importFrom ggheatmap ggheatmap
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom stats hclust
#' @importFrom stats dist
#' @export

plot_BioDescr <- function(BioDescr, type="heatmap", cluster_k=2, colors=c("#2E9FDF", "#E7B800")){

  BioDescr2 <- to_biodescr(BioDescr)

  if (type == "heatmap") {
    ggheatmap(BioDescr2$bd,
            color = c("white", colors[1]),
            border = "grey")+
    theme(axis.text.x = element_text(angle = 90,face = "bold"),
          axis.text.y = element_text(colour = "red",face = "bold"),
          legend.position = "none")
  }

  if (type == "clusterplot") {
    BioDescr2$bd %>%
      dist(method = "euclidean") %>%
      hclust(method = "ward.D2") %>%
    fviz_dend(k = cluster_k,
              cex = 0.5,
              k_colors = colors,
              color_labels_by_k = TRUE,
              rect = TRUE )
  }

}
