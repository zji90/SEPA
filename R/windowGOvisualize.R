#' windowGOvisualize
#' 
#' Visualize results of GO analysis with moving window for transition patterns
#'
#' This function is specifically designed to visualize the results obtained from windowGOanalysis function. Users can choose to visualize specific GO terms or 
#' top GO terms in each time window
#' 
#' @param GOres The direct output of windowGOanalysis function.
#' @param GOTerm The name of GO term to be displayed. If NULL, top GO terms will be displayed instead
#' @param topterm The number of top GO terms to be displayed. This argument only works when GOTerm is NULL.
#' @param mode To plot in heatmap or line graph. Either "Heatmap" or "Line".
#' @return A ggplot2 object.
#' @import ggplot2
#' @importFrom grid unit
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' library(topGO)
#' data(HSMMdata)
#' pattern <- pseudotimepattern(HSMMdata,pseudotime)
#' windowGOvisualize(windowGOanalysis(pattern,type="constant_up"))

windowGOvisualize <- function(GOres,GOTerm=NULL,topterm=2,mode="Heatmap") {      
      termnum <- nrow(GOres[[1]])
      if (is.null(GOTerm)) {                  
            GOTerm <- unique(as.vector(sapply(GOres,function(i) i[1:topterm,1])))
      }
      GOdes <- NULL
      for (i in GOres) {
            GOdes <- rbind(GOdes,as.matrix(i[,1:2]))
      }
      GOdes <- unique(GOdes)
      rankres <- sapply(GOTerm,function(term) {
            sapply(GOres,function(i) {
                  tmp <- which(i[,1]==term)
                  if (length(tmp) == 0) 
                        tmp <- termnum + 1
                  tmp
            })
      })                  
      rankres <- melt(rankres,id="V1")
      maxrank <- max(rankres$value)
      if (maxrank < termnum + 1) {
            tmpyset <- scale_y_reverse() 
      } else {
            tickpos <- round(seq(0,maxrank,maxrank/5))
            tickpos <- tickpos[tickpos < termnum + 1]
            tmpyset <- scale_y_reverse(lim=c(termnum + 1,1),breaks=c(tickpos,(termnum+1)),labels=c(tickpos,paste(">",termnum)))
      }
      colnames(rankres)[2] <- "GOTerm"
      rankres[,2] <- paste0(rankres[,2],"\n",sapply(rankres[,2], function(i) GOdes[GOdes[,1]==i,2]))
      if (mode=="Heatmap") {
            p <- ggplot(data=rankres, aes(x=Var1, y=GOTerm)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "blue",high = "red",mid="white")
      } else {
            p <- ggplot(data = rankres, aes(x=Var1, y=value, colour=GOTerm)) +
                  geom_line(aes(group=GOTerm)) + tmpyset +
                  geom_point(size=4)
      }
      p + xlab("Interval") +
            ylab("Rank") +
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
                  legend.position = "right",
                  legend.key.width=unit(3,"line"),legend.key.height=unit(3,"line")
            )
}
