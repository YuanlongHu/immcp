#' Plot smoothed density estimates for adjusted score
#'
#'
#' @title plot_density
#' @param result an object of class ScoreResult.
#' @param drug a character of drug name.
#' @param view the index
#' @param fill fill color.
#' @return a ggplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 labs
#' @export
#' @author Yuanlong Hu


plot_density <- function(result, drug, view="Change_average_degree", fill="#6495ED"){
  dat <- data.frame(drug = rep(drug, nrow(result@adj[[drug]])),
                    var = result@adj[[drug]][,view])
  plot <- ggplot(dat, aes_(x=~var)) +
    geom_density(alpha=0.2, fill=fill)+
    geom_vline(xintercept=result@res[drug,view], color="red", linetype="dashed", size=1)+
    theme_light()+
    labs(x="")
  return(plot)
}
