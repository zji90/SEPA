#' pseudotimevisualize
#' 
#' Visualize Gene Expression Pattern for Pseudo Temporal Cell Ordering
#'
#' Visualize gene expression pattern of one or multiple genes for pseudo temporal cell ordering. For one gene, a scatterplot with fitted lines will be generated.
#' For multiple genes, a heatmap with fitted or true expression values will be generated.
#' 
#' @param pattern The exact output of the pseudotimepattern function.
#' @param gene A character value or vector of gene names. Should be included in the gene expression matrix.
#' @param showtrue For the heatmap of multiple gene, whether to display true gene expression values or fitted gene expression values.
#' @return A ggplot2 object
#' @import segmented reshape2
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(HSMMdata)
#' pattern <- pseudotimepattern(HSMMdata,pseudotime)
#' pseudotimevisualize(pattern,"ENSG00000108821.9")
#' pseudotimevisualize(pattern,c("ENSG00000108821.9","ENSG00000187193.8"))

pseudotimevisualize <- function(pattern,gene,showtrue=F) {      
      findgene <- function(gene) {
            findres <- sapply(names(pattern$pattern), function(pat) {
                  if (grepl("_",pat)) {
                        gene %in% row.names(pattern$pattern[[pat]])
                  } else {
                        gene %in% names(pattern$pattern[[pat]])
                  }
            })
            names(which(findres))
      }
      genename <- gene
      x <- sort(pattern$pseudotime)
      gene <- sapply(gene, findgene)      
      if (length(gene) == 0) {
            stop("Gene is not found")
      } else if (length(gene) == 1) {
            if (grepl("_",gene)) {
                  transinfo <- pattern$pattern[[gene]][names(gene),]
                  transmean <- transinfo[1,1]                  
            } else {
                  transinfo <- pattern$pattern[[gene]][names(gene)]                  
            }            
            y <- pattern$expr[names(gene),names(x)]              
            fity <- pattern$fitexpr[names(gene),names(x)]
            p <- qplot(x,y,data=data.frame(x=x,y=y),size=3)
            gene <- gsub("constant","black",gene)
            gene <- gsub("up","green",gene)
            gene <- gsub("down","red",gene)
            gene <- strsplit(gene,"_")[[1]]
            if (length(gene)==1) {
                  p <- p + geom_line(aes(x=x,y=y),data=data.frame(x=x,y=fity),col=gene,size=2)
            } else {
                  transpointy <- (fity[2]-fity[1])/(x[2]-x[1])*(transmean-x[1])+fity[1]
                  p <- p + geom_line(aes(x=x,y=y),data=data.frame(x=c(x[x<transmean],transmean),y=c(fity[x<transmean],transpointy)),col=gene[1],size=2)
                  p <- p + geom_line(aes(x=x,y=y),data=data.frame(x=c(transmean,x[x>transmean]),y=c(transpointy,fity[x>transmean])),col=gene[2],size=2)
                  ypos <- ifelse(max(y) > 0, 1.1 * max(y), 0.9 * max(y))
                  p <- p + geom_point(aes(x=transmean,y=ypos),col="blue",size=5) + geom_segment(aes(x=transinfo[1,2],xend=transinfo[1,3],y=ypos,yend=ypos),col="blue",size=2)
            }      
            p + xlab("Pseudo-time") +
                  ylab("Gene Expression") +
                  ggtitle(genename) +
                  theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(size=17,color='black'),
                        axis.text.y = element_text(size=17,color='black'),
                        axis.title.x = element_text(size=20),
                        axis.title.y = element_text(size=20),
                        strip.text.y = element_text(size=17,color='black'),        
                        plot.title = element_text(size = 25),
                        legend.position = "none")
      } else {            
            if (!showtrue) {
                  geneexpr <- pattern$fitexpr[names(gene),names(x)]
            } else {
                  geneexpr <- pattern$expr[names(gene),names(pseudotime)]
            }            
            colnames(geneexpr) <- NULL
            geneexpr <- geneexpr[rev(row.names(geneexpr)),]                  
            linesegmat <- NULL
            for (g in 1:length(gene)) {
                  if (grepl("_",gene[g])) {
                        tmp <- unlist(pattern$pattern[[gene[g]]][names(gene[g]),])
                        tmp <- sapply(tmp,function(j) which.min(abs(j-x)))
                        linesegmat <- rbind(linesegmat,c(tmp,g))     
                  }                  
            }
            linesegmat[,4] <- length(gene) + 1 - linesegmat[,4]
            geneexpr <- t(apply(geneexpr,1,scale))
            geneexpr <- melt(geneexpr)                                          
            p <- ggplot(data=geneexpr, aes(x=Var2, y=Var1)) + geom_tile(aes(fill = value)) + scale_fill_gradient2(low = "blue",high = "red",mid="white") +
                  geom_segment(aes(x=x,xend=xend,y=y,yend=yend),data=data.frame(x=linesegmat[,2],xend=linesegmat[,3],y=linesegmat[,4],yend=linesegmat[,4])) +
                  geom_point(aes(x=x,y=y),data=data.frame(x=linesegmat[,1],y=linesegmat[,4]))                        
            p + xlab("Cell") +
                  ylab("Predicted Expression") +
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
                        legend.position = "right" )                       
            
      }      
}