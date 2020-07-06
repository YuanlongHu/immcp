##' Drug target or pathway network visualization
##'
##'
##' @title plot_network
##' @rdname plot_network
##' @param x enrichment result
##' @param Drug number of enriched terms to display
##' @param node_color The node color.
##' @param layout Character Name of network layout function to use. Default to "layout_nicely".
##' @param ... additional parameters
##' @return visNetwork object
##' @export

setGeneric("plot_network",
           function(x, Drug, node_color, layout = "layout_nicely", ...)
             standardGeneric("plot_network")
)
