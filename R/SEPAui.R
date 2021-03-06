#' SEPAui
#' 
#' Launch the SEPA user interface in local machine
#'
#' This function will automatically launch the SEPA user interface in a web browser. 
#' It provides a much easier and more convenient way to analysis gene expression patterns and perform GO analysis for single-cell RNA-seq data.
#' The user interface can also be accessed by http://zhiji.shinyapps.io/SEPA. Neither R nor any packages are required in this online version.
#' However, it is highly recommended that the user interface be launched locally for faster running speed.
#' 
#' @export
#' @import shiny ggplot2 topGO reshape2 segmented
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' \dontrun{
#'    SEPAui()
#' }


SEPAui <- function() {
      shiny::runApp(system.file("shiny",package="SEPA"))
}
