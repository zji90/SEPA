#' truetimepattern
#' 
#' Identify Pattern for True Experiment Time
#'
#' Identify the gene expression patterns for true experiment time. For the expressions of each gene, the function performs t-tests for cells from neighboring time points.
#' The expression pattern for cells from neiboring time points could be increasing, decreasing or constant. All patterns are concatenated using "-" to form the final pattern.
#' 
#' @param expr The matrix of gene expression profile.
#' @param truetime A character data.frame or matrix of true experimental time. First column: Cell name; Second column: experiment time.
#' @param simplify Whether to simplify pattern so that same neiboring patterns will be reduced to one. For example "up_up_constant" will be simplied to "up_constant".
#' @param removeconstant Whether to remove all constant patterns. For example "up_up_constant" will be simplied to "up_up". This step will be performed before simplify.
#' @return A named vector or patterns. The names corresponds to gene names.
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(HSMMdata)
#' truetimepattern(HSMMdata,truetime)

truetimepattern <- function(expr,truetime,simplify=T,removeconstant=F) {
      tmp <- truetime[,2]
      names(tmp) <- truetime[,1]
      truetime <- factor(tmp)      
      apply(expr, 1, function(e) {
            if (sum(expr) == 0) {
                  "zero"
            } else {
                  ttestpval <- sapply(1:(length(levels(truetime))-1), function(i) {
                        sign(mean(e[truetime==levels(truetime)[i]])-mean(e[truetime==levels(truetime)[i+1]]))*t.test(e[truetime==levels(truetime)[i]],e[truetime==levels(truetime)[i+1]])$p.value            
                  })      
                  ttestpval <- sign(ttestpval) * p.adjust(abs(ttestpval),method="fdr")
                  pattern <- rep("constant",length(levels(truetime))-1)
                  pattern[abs(ttestpval) < 0.05 & ttestpval < 0] <- "up"
                  pattern[abs(ttestpval) < 0.05 & ttestpval > 0] <- "down"
                  if (removeconstant) 
                        pattern <- pattern[pattern!="constant"]
                  if (simplify)
                        pattern <- rle(pattern)$values
                  if (length(pattern) == 0) {
                        "NA"
                  } else {
                        paste(pattern,collapse = "_")      
                  }      
            }
      })   
}

