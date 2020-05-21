#' plot smoothed density estimates for adjusted score
#'
#'
#' @title plot_density
#' @param result an object of class ScoreResult.
#' @param drug a character of drug name.
#' @param fill fill color.
#' @return a ggplot
#' @export
#' @author Yuanlong Hu
#' @examples
#' \dontrun{
#'   data("drugSample")
#'   FP <- extrFP(disease_biomarker = drugSample$disease_biomarker,
#'                drug_target = drugSample$herb_target,
#'                geneset = "ImmGenTop150")
#'   res <- score_fp(FP, n=100)
#'   plot_density(res, drug="BAN_XIA_XIE_XIN_TANG")
#' }


plot_density <- function(result, drug, fill="#6495ED" ){
  dat <- data.frame(drug = rep(drug, length(res@adj[[drug]][-1])),
                    var = res@adj[[drug]][-1])
  plot <- ggplot(dat, aes(x=var)) +
    geom_density(alpha=0.2, fill=fill)+
    geom_vline(xintercept=res@adj[[drug]][1], color="red", linetype="dashed", size=1)+
    theme_light()
  return(plot)
}
