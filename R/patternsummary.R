#' patternsummary
#' 
#' Count table for pattern
#'
#' This function generates a count table of number of genes with each pattern.
#' 
#' @param pattern The direct output of truetimepattern or pseudotimepattern function.
#' @return A data.frame object. First column: pattern; Second column: number of genes
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(HSMMdata)
#' pattern <- truetimepattern(HSMMdata,truetime,removeconstant=TRUE)
#' patternsummary(pattern)

patternsummary <- function(pattern) {                
      if (is.list(pattern)) {
            tmp <- sapply(names(pattern$pattern),function(i) {
                  if (grepl("_",i)) {
                        nrow(pattern$pattern[[i]])
                  } else {
                        length(pattern$pattern[[i]])
                  }
            })      
      } else {
            tmp <- table(pattern)            
      }
      data.frame(Pattern=names(tmp),Number=as.vector(tmp))
}

