#' pseudotimepattern
#' 
#' Identify Pattern for Pseudo Temporal Cell Ordering
#'
#' Identify the gene expression patterns for true experiment time. For the expressions of each gene, the function performs t-tests for cells from neighboring time points.
#' The expression pattern for cells from neiboring time points could be increasing, decreasing or constant. All patterns are concatenated using "-" to form the final pattern.
#' 
#' @param expr The matrix of gene expression profile.
#' @param pseudotime A character data.frame or matrix of pseudo-time. First column: Cell name; Second column: pseudo-time.
#' @param simplify Whether to simplify pattern so that same neiboring patterns will be reduced to one. For example "up_up_constant" will be simplied to "up_constant".
#' @param removeconstant Whether to remove all constant patterns. For example "up_up_constant" will be simplied to "up_up". This step will be performed before simplify.
#' @param plot Whether to generate plot for genes with transition points.
#' @param gap Number of first and last gap cells that will be excluded when considering transition points.
#' @return A list. expr: original expression matrix; pseudotime: original pseudotime; pattern: a list containing results of different patterns. For single patterns, it is a named
#' vector where values are the p-values of the t-test of the simple linear regression slope coefficient. The vector is ordered according to the p-values. For transition patterns, 
#' a data.frame containing the mean and confidence interval of the transition point. It is ordered according to the transition points; fitexpr: the fitted expression matrix
#' @import segmented
#' @export
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' data(HSMMdata)
#' pseudotimepattern(HSMMdata,pseudotime)

pseudotimepattern <- function(expr,pseudotime,simplify=T,removeconstant=F,plot=F,gap=10) {
      # presettings
      fastRSS <- function(X, y) {
            XtX <- crossprod(X)
            ch <- chol(XtX)
            V <- chol2inv(ch)
            sum(drop(X %*% V %*% crossprod(X, y)-y)^2)
      }
      x <- pseudotime[,2]
      names(x) <- pseudotime[,1]
      x <- sort(x)
      breakpoint <- rowMeans(cbind(x[-1],x[-length(x)]))
      breakpoint <- breakpoint[gap:(length(breakpoint)-gap+1)]
      spmat <- t(sapply(breakpoint,function(i) {
            ifelse(x > i,x - i,0)
      })) 
      # genes that are constantly zero
      expr <- expr[,names(x)]
      zerogene <- row.names(expr)[rowSums(expr) == 0]
      fitexpr <- expr <- expr[rowSums(expr) > 0,]
      # use davies test to identify whether there is a transition point      
      print("Running davies test")
      set.seed(1234)
      daviespval <- apply(expr,1,function(y) {
            tmpy <- y[(gap+1):(length(y)-gap)]
            tmpx <-x[(gap+1):(length(y)-gap)]
            if (sum(tmpy) > 0) {
                  fit <- lm(tmpy~tmpx)
                  davies.test(fit,~tmpx)$p.value      
            } else {
                  1
            }            
      })
      daviespval[is.na(daviespval)] <- 1
      daviespval <- p.adjust(daviespval,method="fdr")      
      notransgene <- names(daviespval)[daviespval >= 0.05]
      transgene <- names(daviespval)[daviespval < 0.05]
      # calculate potential transition point      
      print("Determining potential transition point positions")
      transpos <- sapply(transgene,function(gene) {
            y <- expr[gene,]
            RSS1 <- sapply(1:length(breakpoint),function(i) {
                  fastRSS(cbind(1,x,spmat[i,]),y)
            })            
            which.min(RSS1)
      }) 
      pattern <- data.frame(pattern=rep("constant",length(transgene)),transpoint=rep(0,length(transgene)),LCI=rep(0,length(transgene)),UCI=rep(0,length(transgene)),stringsAsFactors = F)      
      row.names(pattern) <- transgene
      print("Fitting segmented regression models")
      for (i in 1:length(transgene)) {
            set.seed(1234)
            y <- expr[transgene[i],]
            fit <- lm(y~x)
            o.seg1 <- tryCatch(segmented(fit,seg.Z=~x,psi=breakpoint[transpos[i]]),error=function(e) {}) 
            if (is.null(o.seg1)) {
                  suppressWarnings(o.seg1 <- segmented(fit,seg.Z=~x,psi=breakpoint[transpos[i]],control = seg.control(n.boot=0,it.max=1)))            
            }      
            if (o.seg1$psi[2] <= x[gap] || o.seg1$psi[2] >= x[length(x)-gap+1]) {
                  notransgene <- c(notransgene,transgene[i])
            } else {
                  slopecol <- rep("black",2)
                  for (j in 1:2) {
                        if (slope(o.seg1)$x[j,4] * slope(o.seg1)$x[j,5] > 0) {
                              if (slope(o.seg1)$x[j,4] > 0) {
                                    slopecol[j] <- "green"
                              } else {
                                    slopecol[j] <- "red"
                              }
                        }
                  }
                  if (removeconstant) {
                        slopecol <- slopecol[slopecol!="black"]
                  }
                  if (length(slopecol) < 2 | (simplify & slopecol[1] == slopecol[2])) {
                        notransgene <- c(notransgene,transgene[i])
                  } else {
                        if (plot) {
                              plot(x,y,pch=20,main=transgene[i])
                              plot(o.seg1,add=T,link=F,lwd=3,col=slopecol,rug=F)
                              lines(o.seg1,col="blue",pch=19,bottom=FALSE,lwd=2)
                        }                  
                        pattern[i,1] <- paste0(slopecol,collapse = "_")
                        pattern[i,2] <- confint(o.seg1)$x[1]
                        pattern[i,3] <- confint(o.seg1)$x[2]
                        pattern[i,4] <- confint(o.seg1)$x[3]      
                        fitexpr[transgene[i],] <- fitted(o.seg1)
                  }      
            }                        
      }
      pattern <- pattern[pattern[,1]!="constant",]
      pattern[,1] <- gsub("black","constant",pattern[,1])
      pattern[,1] <- gsub("green","up",pattern[,1])
      pattern[,1] <- gsub("red","down",pattern[,1])
      pattern <- pattern[order(pattern$pattern,pattern$transpoint),]
      tmp <- pattern[,1]
      pattern <- pattern[,-1]
      pattern <- split(pattern,tmp)
      if (length(zerogene) > 0)
            pattern$zero <- zerogene      
      notranspval <- sapply(notransgene,function(i) {            
            y <- expr[i,]
            fit <- lm(y~x)
            fitexpr[i,] <<- fitted(fit)
            sign(summary(fit)[[4]][2,1])*summary(fit)[[4]][2,4]            
      })
      
      notranspval <- sign(notranspval) * p.adjust(abs(notranspval),method = "fdr")
      upgene <- notranspval[notranspval > 0 & abs(notranspval) < 0.05]
      downgene <- notranspval[notranspval < 0 & abs(notranspval) < 0.05]
      constantgene <- notranspval[abs(notranspval) >= 0.05]
      pattern$up <- sort(upgene)
      pattern$down <- sort(abs(downgene))
      pattern$constant <- sort(abs(constantgene))
      list(expr=expr,pseudotime=x,pattern=pattern,fitexpr=fitexpr)
}
