#' truetimevisualize
#' 
#' Visualize Gene Expression Pattern for True Experiment Time
#'
#' Identify the gene expression patterns for true experiment time. For the expressions of each gene, the function performs t-tests for cells from neighboring time points.
#' The expression pattern for cells from neiboring time points could be increasing, decreasing or constant. All patterns are concatenated using "-" to form the final pattern.
#' 
#' @param expr The matrix of gene expression profile.
#' @param truetime A character data.frame or matrix of true experimental time. First column: Cell name; Second column: experiment time.
#' @param gene A vector of gene names to be plotted.
#' @param mode A character value specifying mean or median to be displayed
#' @return A ggplot2 object.
#' @import ggplot2
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(HSMMdata)
#' truetimevisualize(HSMMdata,truetime,c("ENSG00000122180.4","ENSG00000125968.7"))

truetimevisualize <- function(expr,truetime,gene,mode=c("mean","median")) {
      tmp <- truetime[,2]
      names(tmp) <- truetime[,1]
      truetime <- factor(tmp)      
      tmp <- data.frame(time=rep(levels(truetime),length(gene)),expmean=0,Gene=rep(gene,each=length(levels(truetime))))
      tmp$expmean <- as.vector(sapply(gene, function(singlegene) {
            e <- expr[singlegene,names(truetime)]
            sapply(1:length(levels(truetime)), function(i) {
                  if (mode == "mean") {
                        mean(e[truetime==levels(truetime)[i]])      
                  } else {
                        median(e[truetime==levels(truetime)[i]])
                  }                  
            })
      }))           
      ggplot(data = tmp, aes(x=time, y=expmean, colour=Gene)) +
            geom_line(aes(group=Gene)) +
            geom_point(size=4) +             
            xlab("Experiment Time") +
            ylab("Mean Expression Values") +
            theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.text.x = element_text(size=17,color='black'),
                  axis.text.y = element_text(size=17,color='black'),
                  axis.title.x = element_text(size=20,vjust=-1),
                  axis.title.y = element_text(size=20,vjust=1),
                  strip.text.y = element_text(size=17,color='black'),                  
                  legend.text = element_text(size=15),
                  legend.title = element_text(size=15),
                  legend.position = "right"                  
            )
}

