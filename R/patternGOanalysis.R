#' patternGOanalysis
#' 
#' Performs GO analysis on genes for each pattern
#'
#' This function is basically a wrap up of functions in the topGO package. It takes as input the output of truetimepattern or pseudotimepattern function.
#' For each pattern, the GO enrichment analysis is performed where input genes are genes with specific patterns and background genes are all genes in the 
#' expression profile. Users should correctly select identifier and species, otherwise the function may breakdown.
#' 
#' @param pattern The direct output of truetimepattern or pseudotimepattern function.
#' @param type Character value of specific pattern to perform the GO analysis. If NULL GO analysis will be performed for all patterns.
#' @param termnum Number of top GO terms to be displayed.
#' @param identifier The identifier of the genes. It should be one of the following: "Entrez", "GenBank", "Alias", "Ensembl", "Gene", "Symbol", "GeneName" and "UniGene"
#' @param species The species of the genes. Currently only "Human" and "Mouse" are supported
#' @return A list where each element is a data.frame containing the results of GO analysis.
#' @export
#' @import topGO org.Hs.eg.db org.Mm.eg.db
#' @author Zhicheng Ji, Hongkai Ji <zji4@@zji4.edu>
#' @examples
#' library(topGO)
#' data(HSMMdata)
#' pattern <- truetimepattern(HSMMdata,truetime,removeconstant=TRUE)
#' patternGOanalysis(pattern,termnum=20,identifier="ENSEMBL",species="Human")

patternGOanalysis <- function(pattern,type=NULL,termnum=20,identifier="ENSEMBL",species="Human") {
      if (species == "Human") {
            mapdb <- "org.Hs.eg.db"      
      } else if (species == "Mouse") {
            mapdb <- "org.Mm.eg.db"      
      }                   
      if (is.list(pattern)) {
            allgene <- NULL
            for (i in names(pattern$pattern)) {
                  if (grepl("_",i)) {
                        gene <- row.names(pattern$pattern[[i]])                        
                        tmp <- rep(i,length(gene))
                        names(tmp) <- gene
                        allgene <- c(allgene,tmp)
                  } else {
                        gene <- names(pattern$pattern[[i]])                        
                        tmp <- rep(i,length(gene))
                        names(tmp) <- gene
                        allgene <- c(allgene,tmp)
                  }
            }
            pattern <- allgene
      }
      if (identifier == "ENSEMBL") {            
            names(pattern) <- sapply(names(pattern),function(i) strsplit(i,"\\.")[[1]][1])                              
      }
      allgene <- names(pattern)
      res <- list()
      if (is.null(type)) {
            type <- unique(pattern)
      } else {
            type <- intersect(type,unique(pattern))
      }
      if (length(type) > 0) {
            for (i in type) {
                  inputgene <- names(pattern)[pattern==i]
                  geneList <- factor(as.integer(allgene %in% inputgene))
                  names(geneList) <- allgene
                  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = mapdb, ID = identifier)
                  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
                  res[[i]] <- GenTable(GOdata, classicFisher = resultFisher, topNodes = termnum,orderBy="classicFisher")      
            }
            res      
      } else {
            stop("No corresponding pattern is found!")     
      }      
}


