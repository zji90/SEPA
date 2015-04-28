#' windowGOanalysis
#' 
#' Performs GO analysis with moving window for transition patterns
#'
#' This function is specifically designed for GO analysis of genes with a specific transition pattern. GO analyses are performed iteratively on a group of genes 
#' with similar transition points. Users can define the windowsize (the number of genes in each group) and the movesize (how many genes to move forward each time).
#' 
#' @param pattern The direct output of pseudotimepattern function.
#' @param type The type of transition pattern.
#' @param windowsize The number of genes in each group.
#' @param movesize How many genes to move forward each time.
#' @param termnum Number of top GO terms to be displayed.
#' @param identifier The identifier of the genes. It should be one of the following: "Entrez", "GenBank", "Alias", "Ensembl", "Gene", "Symbol", "GeneName" and "UniGene"
#' @param species The species of the genes. Currently only "Human" and "Mouse" are supported
#' @return A list where each element is a data.frame containing the results of GO analysis. The name of the list specifies the group of genes.
#' @export
#' @import topGO org.Hs.eg.db org.Mm.eg.db
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' library(topGO)
#' data(HSMMdata)
#' pattern <- pseudotimepattern(HSMMdata,pseudotime)
#' windowGOanalysis(pattern,type="constant_up")

windowGOanalysis <- function(pattern, type = "constant_up", windowsize = NULL, movesize = NULL, termnum = 20, identifier="ENSEMBL", species="Human") {
      if (species == "Human") {
            mapdb <- "org.Hs.eg.db"      
      } else if (species == "Mouse") {
            mapdb <- "org.Mm.eg.db"      
      }         
      if (is.null(windowsize)) {
            windowsize <- round(nrow(pattern$pattern[[type]]) / 2)
      } 
      if (is.null(movesize)) {
            movesize <- round(windowsize / 2)
      }
      maxwinnum <- ceiling(max(0,(nrow(pattern$pattern[[type]])-windowsize))/movesize + 1)
      res <- list()
      for (i in 1:maxwinnum) {
            if (i==maxwinnum) {
                  startid <- max(1,nrow(pattern$pattern[[type]])-windowsize+1)
                  endid <- nrow(pattern$pattern[[type]])                
            } else {
                  startid <- 1+movesize*(i-1)
                  endid <- movesize*(i-1)+windowsize
            }
            inputgene <- row.names(pattern$pattern[[type]])[startid:endid]
            allgene <- row.names(pattern$expr)
            
            if (identifier == "ENSEMBL") {            
                  inputgene <- sapply(inputgene,function(i) strsplit(i,"\\.")[[1]][1])                              
                  allgene <- sapply(allgene,function(i) strsplit(i,"\\.")[[1]][1])                              
            }
            geneList <- factor(as.integer(allgene %in% inputgene))
            names(geneList) <- allgene
            GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = mapdb, ID = identifier)
            resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
            res[[paste0(startid,"-",endid)]] <- GenTable(GOdata, classicFisher = resultFisher, topNodes = termnum,orderBy="classicFisher")                  
      }
      res      
}

