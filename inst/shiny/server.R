######################################################
##                       SEPA                       ##
##             Interactive User Interface           ##
##                     Server File                  ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################
library(shiny)
library(grid)
library(ggplot2)
library(topGO)
library(reshape2)
library(segmented)

options(shiny.maxRequestSize=30000*1024^2)

shinyServer(function(input, output,session) {
      
      output$showbusybar <- renderUI({
            tagList(
                  tags$head(
                        tags$link(rel="stylesheet", type="text/css",href="style.css"),
                        tags$script(type="text/javascript", src = "busy.js")
                  ),                  
                  div(class = "busy",  
                      p("Calculation in progress.."), 
                      img(src="ajaxloaderq.gif")
                  )
            )        
      })
      
      Maindata <- reactiveValues()
      
      ### Readexp ###
      observe({
            if (input$Readexpreadin > 0)
                  isolate({
                        FileHandle <- input$ReadexpFile
                        if (!is.null(FileHandle)) {
                              tmpdata <- read.table(FileHandle$datapath,header=T,sep=input$Readexpsep,as.is=T,blank.lines.skip=TRUE)
                              Maindata$rawexpr <- as.matrix(tmpdata[,-1])
                              row.names(Maindata$rawexpr) <- make.names(tmpdata[,1])
                        }
                  })
      })
      
      output$Readexpraddsymbolui <- renderUI({
            if (input$Readexprchooseanno!="symbol") 
                  checkboxInput("Readexpraddsymbol","Add gene symbol to row names",value=F)
      })
      
      observe({
            if (!is.null(Maindata$rawexpr)) {
            tmp <- Maindata$rawexpr
            if (input$Readexpraddsymbol) {                  
                  if (input$Readexprchoosespecies == "Human") {
                        x <- org.Hs.egSYMBOL
                        mapped_genes <- mappedkeys(x)
                        eg2symbol <- unlist(as.list(x[mapped_genes]))
                        if (input$Readexprchooseanno == "entrez") {
                              row.names(tmp) <- paste0(eg2symbol[row.names(tmp)],"_",row.names(tmp))
                        } else if (input$Readexprchooseanno == "ensembl") {                              
                              xx <- sapply(as.list(org.Hs.egENSEMBL2EG),function(i) i[1])
                              ensembl2symbol <- eg2symbol[xx]
                              names(ensembl2symbol) <- names(xx)
                              tmprowname <- sapply(row.names(tmp),function(i) strsplit(i,"\\.")[[1]][1])                              
                              row.names(tmp) <- paste0(ensembl2symbol[tmprowname],"_",row.names(tmp))                              
                        } else if (input$Readexprchooseanno == "unigene") {                              
                              x <- org.Hs.egUNIGENE2EG
                              mapped_genes <- mappedkeys(x)
                              unigene2eg <- sapply(as.list(x[mapped_genes]),function(i) i[1])
                              unigene2symbol <- eg2symbol[unigene2eg]
                              names(unigene2symbol) <- names(unigene2eg)
                              row.names(tmp) <- paste0(unigene2symbol[row.names(tmp)],"_",row.names(tmp))
                        }
                  } else if (input$Readexprchoosespecies == "Mouse") {
                        x <- org.Mm.egSYMBOL
                        mapped_genes <- mappedkeys(x)
                        eg2symbol <- unlist(as.list(x[mapped_genes]))
                        if (input$Readexprchooseanno == "entrez") {
                              row.names(tmp) <- paste0(eg2symbol[row.names(tmp)],"_",row.names(tmp))
                        } else if (input$Readexprchooseanno == "ensembl") {                              
                              xx <- sapply(as.list(org.Mm.egENSEMBL2EG),function(i) i[1])
                              ensembl2symbol <- eg2symbol[xx]
                              names(ensembl2symbol) <- names(xx)
                              tmprowname <- sapply(row.names(tmp),function(i) strsplit(i,"\\.")[[1]][1])
                              row.names(tmp) <- paste0(ensembl2symbol[tmprowname],"_",row.names(tmp))
                        } else if (input$Readexprchooseanno == "unigene") {                              
                              x <- org.Mm.egUNIGENE2EG
                              mapped_genes <- mappedkeys(x)
                              unigene2eg <- sapply(as.list(x[mapped_genes]),function(i) i[1])
                              unigene2symbol <- eg2symbol[unigene2eg]
                              names(unigene2symbol) <- names(unigene2eg)
                              row.names(tmp) <- paste0(unigene2symbol[row.names(tmp)],"_",row.names(tmp))
                        }                        
                  }                  
            }
            Maindata$nametransexpr <- tmp
            }
      })
      
      observe({
            if (!is.null(Maindata$nametransexpr)) {
                  tmp <- Maindata$nametransexpr[,!colnames(Maindata$nametransexpr) %in% input$Readexpcell]
                  if (input$Readexplogtf) {
                        if (input$Readexplogbase == "2") {
                              tmp <- log2(tmp+as.numeric(input$Readexplogpseudocount))
                        } else if (input$Readexplogbase == "10") {
                              tmp <- log10(tmp+as.numeric(input$Readexplogpseudocount))
                        } else if (input$Readexplogbase == "e") {
                              tmp <- log(tmp+as.numeric(input$Readexplogpseudocount))
                        }
                  }                  
                  tmp <- tmp[rowSums(tmp) > 0,]
                  tmp <- tmp[rowSums(tmp > as.numeric(input$Readexpgeneexpcutoff)) > ncol(tmp)*as.numeric(input$Readexpgenepercentcutoff)/100,]
                  tmp <- tmp[apply(tmp,1,sd)/rowMeans(tmp) > as.numeric(input$Readexpgenecvcutoff),]
                  Maindata$expr <- tmp      
            }            
      })
      
      observe({
            if (input$Readexprchoosespecies=="Human") {
                  library(org.Hs.eg.db)
                  mapdb <- "org.Hs.eg.db"
            } else {
                  library(org.Mm.eg.db)
                  mapdb <- "org.Mm.eg.db"
            }
            xx <- unique(unlist(annFUN.org("BP", mapping = mapdb, ID = input$Readexprchooseanno)))
            genename <- row.names(Maindata$rawexpr)
            if (input$Readexprchooseanno == "ensembl") {
                  genename <- sapply(genename,function(i) strsplit(i,"\\.")[[1]][1])
            }
            Maindata$genenumwithid <- length(intersect(xx,genename))
      })
      
      output$Readexpcellui <- renderUI(
            selectizeInput("Readexpcell","Select cells to be removed",colnames(Maindata$rawexpr),multiple = TRUE),            
      )
      
      output$Readexpshowrawtable <- renderTable(head(Maindata$expr))
      
      output$Readexpshowsummaryui <- renderUI({
            if (!is.null(Maindata$expr)) {
                  tagList(                        
                        h4("Summary of the gene expression profile:"),                        
                        h5(paste("The filtered dataset contains",nrow(Maindata$expr),"genes and",ncol(Maindata$expr),"cells")),                        
                        h5(paste(Maindata$genenumwithid,"genes have corresponding GO terms.")),
                        hr(),
                        h4("Head of the input file:"),
                        tableOutput("Readexpshowrawtable")                                                
                  )
            }
      })
      
      ### Truetime ###
      
      observe({
            if (input$Truetimereadin > 0)
                  isolate({
                        FileHandle <- input$TruetimeFile
                        if (!is.null(FileHandle)) {
                              tmp <- read.table(FileHandle$datapath,header=T,sep=input$Truetimesep,as.is=T,blank.lines.skip=TRUE)                              
                              commoncell <- intersect(tmp[,1],colnames(Maindata$expr))
                              Maindata$Truetimedata <- tmp[tmp[,1] %in% commoncell,]
                              Maindata$Truetimeexpr <- Maindata$expr[,commoncell]
                        }
                  })
      })
      
      output$Truetimeshowtime <- renderDataTable({
            Maindata$Truetimedata
      })
      
      observe({
            if (!is.null(Maindata$Truetimeexpr) && !is.null(Maindata$Truetimedata) && input$Truetimerunanalysis > 0)
                  isolate({
                        tmptime <- Maindata$Truetimedata[,2]
                        names(tmptime) <- Maindata$Truetimedata[,1]                        
                        truetime <- factor(tmptime)
                        
                        tmppattern <- rep("tmp",nrow(Maindata$Truetimeexpr))
                        withProgress(message = 'Calculation in Progress...', {                                  
                              for (i in 1:nrow(Maindata$Truetimeexpr)) {
                                    e <- Maindata$Truetimeexpr[i,]
                                    ttestpval <- sapply(1:(length(levels(truetime))-1), function(i) {
                                          sign(mean(e[truetime==levels(truetime)[i]])-mean(e[truetime==levels(truetime)[i+1]]))*t.test(e[truetime==levels(truetime)[i]],e[truetime==levels(truetime)[i+1]])$p.value            
                                    })      
                                    ttestpval <- sign(ttestpval) * p.adjust(abs(ttestpval),method="fdr")
                                    pattern <- rep("constant",length(levels(truetime))-1)
                                    pattern[abs(ttestpval) < as.numeric(input$Truetimepvalcutoff) & ttestpval < 0] <- "up"
                                    pattern[abs(ttestpval) < as.numeric(input$Truetimepvalcutoff) & ttestpval > 0] <- "down"                                    
                                    pattern <- paste(pattern,collapse = "_")
                                    tmppattern[i] <- pattern                                    
                                    incProgress(1/nrow(Maindata$Truetimeexpr), detail = paste("Calculating Pattern for Gene", i))                                    
                              }
                        })
                        names(tmppattern) <- row.names(Maindata$Truetimeexpr)
                        Maindata$Truetimeoripattern <- tmppattern
                  })
      })
      
      observe({
            if (!is.null(Maindata$Truetimeoripattern)) {
                  if (!input$Truetimeignoreconst && !input$Truetimesimplify) {
                        Maindata$Truetimepattern <- Maindata$Truetimeoripattern
                  } else {
                        Maindata$Truetimepattern <- sapply(Maindata$Truetimeoripattern, function(pattern) {
                              pattern <- strsplit(pattern,"_")[[1]]
                              if (input$Truetimeignoreconst) 
                                    pattern <- pattern[pattern!="constant"]
                              if (input$Truetimesimplify)
                                    pattern <- rle(pattern)$values    
                              paste(pattern,collapse = "_")
                        })                        
                  }
                  Maindata$Truetimenoconstantname <- names(Maindata$Truetimepattern)[grepl("up",Maindata$Truetimepattern) | grepl("down",Maindata$Truetimepattern)]
            }            
      })
      
      output$Truetimepatternsummary <- renderDataTable({
            if (!is.null(Maindata$Truetimepattern)) {
                  tmp <- data.frame(table(Maindata$Truetimepattern))
                  colnames(tmp) <- c("Pattern","Number")
                  tmp <- tmp[tmp$Pattern!="",]
                  tmp      
            }            
      })
      
      output$Truetimepatternsummarysave <- downloadHandler(
            filename = function() { "Pattern Summary.csv" },
            content = function(file) {                  
                  tmp <- data.frame(table(Maindata$Truetimepattern))
                  colnames(tmp) <- c("Pattern","Number")
                  tmp <- tmp[tmp$Pattern!="",]                  
                  write.csv(tmp,file=file,quote=F,row.names=F)                  
            }
      )
      
      output$TruetimepatternsummaryGeneselectui <- renderUI({
            if (!is.null(Maindata$Truetimepattern))
                  tagList(
                        checkboxInput("TruetimepatternsummaryGeneselectlistall","List Genes for all Pattern",value=T),
                        conditionalPanel(condition = "input.TruetimepatternsummaryGeneselectlistall==0",
                                         selectInput("TruetimepatternsummaryGeneselect","Select Pattern",sort(unique(Maindata$Truetimepattern)),multiple = T)
                        )
                  )            
      })
      
      output$TruetimepatternsummaryGene <- renderDataTable({
            if (!is.null(input$TruetimepatternsummaryGeneselectlistall)) {
                  if (input$TruetimepatternsummaryGeneselectlistall) {
                        tmp <- Maindata$Truetimepattern
                  } else {
                        tmp <- Maindata$Truetimepattern[Maindata$Truetimepattern %in% input$TruetimepatternsummaryGeneselect]      
                  }            
                  data.frame(Gene=names(tmp),Pattern=tmp)
            }
      })      
      
      output$TruetimepatternsummaryGenesave <- downloadHandler(
            filename = function() { "Gene List.csv" },
            content = function(file) {
                  if (input$TruetimepatternsummaryGeneselectlistall) {
                        tmp <- Maindata$Truetimepattern
                  } else {
                        tmp <- Maindata$Truetimepattern[Maindata$Truetimepattern %in% input$TruetimepatternsummaryGeneselect]      
                  }            
                  tmp <- data.frame(Gene=names(tmp),Pattern=tmp)                  
                  write.csv(tmp,file=file,quote=F,row.names=F)                  
            }
      )
      
      output$Truetimeselectgeneui <- renderUI(            
            selectizeInput("Truetimeselectgene","Select genes of interest",Maindata$Truetimenoconstantname,multiple = TRUE)
      )      
      
      output$Truetimevisualize <- renderPlot({     
            if (!is.null(Maindata$Truetimedata) && !is.null(input$Truetimeselectgene)) {
                  tmptime <- Maindata$Truetimedata[,2]
                  names(tmptime) <- Maindata$Truetimedata[,1]                        
                  truetime <- factor(tmptime)
                  tmp <- data.frame(time=rep(levels(truetime),length(input$Truetimeselectgene)),expmean=0,Gene=rep(input$Truetimeselectgene,each=length(levels(truetime))))
                  tmp$expmean <- as.vector(sapply(input$Truetimeselectgene, function(gene) {
                        e <- Maindata$Truetimeexpr[gene,Maindata$Truetimedata[,1]]
                        sapply(1:length(levels(truetime)), function(i) {
                              mean(e[truetime==levels(truetime)[i]])
                        })
                  }))           
                  ggplot(data = tmp, aes(x=time, y=expmean, colour=Gene)) +
                        geom_line(aes(group=Gene)) +
                        geom_point(size=4) +             
                        xlab("Experiment Time") +
                        ylab("Expression Values") +
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
            
      })
      
      output$Truetimevisualizesave <- downloadHandler(
            filename = function() { "Pattern.pdf" },
            content = function(file) {
                  if (!is.null(Maindata$Truetimedata) && !is.null(input$Truetimeselectgene)) {
                        pdf(file,width=14,height=7)
                        tmptime <- Maindata$Truetimedata[,2]
                        names(tmptime) <- Maindata$Truetimedata[,1]                        
                        truetime <- factor(tmptime)
                        tmp <- data.frame(time=rep(levels(truetime),length(input$Truetimeselectgene)),expmean=0,Gene=rep(input$Truetimeselectgene,each=length(levels(truetime))))
                        tmp$expmean <- as.vector(sapply(input$Truetimeselectgene, function(gene) {
                              e <- Maindata$Truetimeexpr[gene,Maindata$Truetimedata[,1]]
                              sapply(1:length(levels(truetime)), function(i) {
                                    mean(e[truetime==levels(truetime)[i]])
                              })
                        }))           
                        tmp <- ggplot(data = tmp, aes(x=time, y=expmean, colour=Gene)) +
                              geom_line(aes(group=Gene)) +
                              geom_point(size=4) +             
                              xlab("Experiment Time") +
                              ylab("Expression Values") +
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
                        print(tmp)
                        dev.off()
                  }
            }
      )
      
      output$Truetimevisualizetable <- renderDataTable({
            if (!is.null(input$Truetimeselectgene)) {
                  data.frame(Gene=input$Truetimeselectgene,Pattern=Maindata$Truetimepattern[input$Truetimeselectgene])
            }
      })
      
      output$Truetimevisualizetablesave <- downloadHandler(
            filename = function() { "Gene Pattern.csv" },
            content = function(file) {                  
                  tmp <- data.frame(Gene=input$Truetimeselectgene,Pattern=Maindata$Truetimepattern[input$Truetimeselectgene])
                  write.csv(tmp,file=file,quote=F,row.names=F)                                                      
            }
      )
      
      output$TruetimeGOanalysisselectui <- renderUI({
            if (!is.null(Maindata$Truetimepattern))
                  selectInput("TruetimeGOanalysisselect","Select pattern of interest",sort(unique(Maindata$Truetimepattern)),multiple = T)
      })
      
      observe({
            if (input$TruetimeGOanalysisrunbutton > 0) {
                  isolate({
                        if (!is.null(input$TruetimeGOanalysisselect) && Maindata$genenumwithid > 0) {
                              withProgress(message = 'Performing GO analysis...', {
                                    if (input$Readexprchoosespecies=="Human") {                        
                                          mapdb <- "org.Hs.eg.db"
                                    } else {                        
                                          mapdb <- "org.Mm.eg.db"
                                    }                  
                                    allgene <- row.names(Maindata$rawexpr)
                                    if (input$Readexprchooseanno == "ensembl") {
                                          allgene <- sapply(allgene,function(i) strsplit(i,"\\.")[[1]][1])
                                    }                   
                                    inputgene <- allgene[Maindata$Truetimepattern %in% input$TruetimeGOanalysisselect]                              
                                    geneList <- factor(as.integer(allgene %in% inputgene))
                                    names(geneList) <- allgene
                                    Maindata$TruetimeGOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = mapdb, ID = input$Readexprchooseanno)
                                    Maindata$TruetimeGOresultFisher <- runTest(Maindata$TruetimeGOdata, algorithm = "classic", statistic = "fisher")
                              })
                        }                        
                  })
            }            
      })
      
      output$TruetimeGOanalysisresult <- renderDataTable({
            if (!is.null(Maindata$TruetimeGOresultFisher)) {
                  tmp <- GenTable(Maindata$TruetimeGOdata, classicFisher = Maindata$TruetimeGOresultFisher, topNodes = as.numeric(input$TruetimeGOanalysistermnum),orderBy="classicFisher")                  
                  tmp
            }
      })
      
      output$TruetimeGOanalysisresultsave <- downloadHandler(
            filename = function() { "GO Analysis.csv" },
            content = function(file) {
                  if (!is.null(Maindata$TruetimeGOresultFisher)) {
                        tmp <- GenTable(Maindata$TruetimeGOdata, classicFisher = Maindata$TruetimeGOresultFisher, topNodes = as.numeric(input$TruetimeGOanalysistermnum),orderBy="classicFisher")                  
                        tmp$Term <- gsub(","," ",tmp$Term)
                        write.csv(tmp,file=file,quote=F,row.names=F)    
                  }
            }
      )
      
      ### Pseudotime ###
      
      observe({
            if (input$Pseudotimereadin > 0)
                  isolate({
                        FileHandle <- input$PseudotimeFile
                        if (!is.null(FileHandle)) {
                              tmp <- read.table(FileHandle$datapath,header=T,sep=input$Pseudotimesep,as.is=T,blank.lines.skip=TRUE)                              
                              pseudotime <- tmp[,2]
                              names(pseudotime) <- tmp[,1]                                                                        
                              commoncell <- intersect(tmp[,1],colnames(Maindata$expr))
                              Maindata$Pseudotimedata <- pseudotime[commoncell]
                              Maindata$Pseudotimeexpr <- Maindata$expr[,commoncell]
                        }
                  })
      })
      
      output$Pseudotimeshowtime <- renderDataTable({
            data.frame(Cell=names(Maindata$Pseudotimedata),Time=Maindata$Pseudotimedata)
      })
      
      observe({            
            if (!is.null(Maindata$Pseudotimeexpr) && !is.null(Maindata$Pseudotimedata) && input$Pseudotimerunanalysis > 0)
                  isolate({                                                 
                        pseudotime <- Maindata$Pseudotimedata
                        expr <- Maindata$Pseudotimeexpr
                        gap <- as.numeric(input$Pseudotimegapvalue)
                        
                        fastRSS <- function(X, y) {
                              XtX <- crossprod(X)
                              ch <- chol(XtX)
                              V <- chol2inv(ch)
                              sum(drop(X %*% V %*% crossprod(X, y)-y)^2)
                        }                        
                        x <- pseudotime                  
                        breakpoint <- rowMeans(cbind(pseudotime[-1],pseudotime[-length(pseudotime)]))                                          
                        breakpoint <- breakpoint[gap:(length(breakpoint)-gap+1)]                        
                        spmat <- t(sapply(breakpoint,function(i) {
                              ifelse(x > i,x - i,0)
                        })) 
                        # genes that are constantly zero                        
                        expr <- expr[,names(pseudotime)]
                        zerogene <- row.names(expr)[rowSums(expr) == 0]
                        fitexpr <- expr <- expr[rowSums(expr) > 0,]                        
                        # use davies test to identify whether there is a transition point     
                        withProgress(message = 'Determining existence of transition points...', {                                  
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
                        })
                        daviespval[is.na(daviespval)] <- 1
                        daviespval <- p.adjust(daviespval,method="fdr")      
                        notransgene <- names(daviespval)[daviespval >= 0.05]
                        transgene <- names(daviespval)[daviespval < 0.05]
                        withProgress(message = 'Predicting starting positions...', {     
                              transpos <- sapply(transgene,function(gene) {
                                    y <- expr[gene,]
                                    RSS1 <- sapply(1:length(breakpoint),function(i) {
                                          fastRSS(cbind(1,x,spmat[i,]),y)
                                    })            
                                    which.min(RSS1)
                              }) 
                        })
                        pattern <- data.frame(pattern=rep("constant",length(transgene)),transpoint=rep(0,length(transgene)),LCI=rep(0,length(transgene)),UCI=rep(0,length(transgene)),stringsAsFactors = F)      
                        row.names(pattern) <- transgene
                        withProgress(message = 'Calculating transition points...', {     
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
                                    slopecol <- rep("constant",2)
                                    for (j in 1:2) {
                                          if (slope(o.seg1)$x[j,4] * slope(o.seg1)$x[j,5] > 0) {
                                                if (slope(o.seg1)$x[j,4] > 0) {
                                                      slopecol[j] <- "up"
                                                } else {
                                                      slopecol[j] <- "down"
                                                }
                                          }
                                    }
                                    if ((input$Pseudotimesimplify && slopecol[1] == slopecol[2]) || (input$Pseudotimeignoreconst && (slopecol[1]=="constant" || slopecol[2]=="constant"))) {
                                          notransgene <- c(notransgene,transgene[i])
                                    } else {                                    
                                          pattern[i,1] <- paste0(slopecol,collapse = "_")
                                          pattern[i,2] <- confint(o.seg1)$x[1]
                                          pattern[i,3] <- confint(o.seg1)$x[2]
                                          pattern[i,4] <- confint(o.seg1)$x[3]      
                                          fitexpr[transgene[i],] <- fitted(o.seg1)
                                    }         
                                    }
                                    incProgress(1/length(transgene), detail = paste("Calculating transition point for gene", i))                                    
                              }
                              
                        })
                        pattern <- pattern[pattern[,1]!="constant",]                        
                        pattern <- pattern[order(pattern$pattern,pattern$transpoint),]
                        tmp <- pattern[,1]
                        pattern <- pattern[,-1]
                        pattern <- split(pattern,tmp)                                                                  
                        notranspval <- sapply(notransgene,function(i) {
                              set.seed(1234)
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
                        if (!input$Pseudotimeignoreconst)
                              pattern$constant <- sort(abs(constantgene))
                        Maindata$Pseudotimefitexpr <- fitexpr
                        Maindata$Pseudotimepattern <- pattern   
                        tmp <- NULL
                        for (name in names(Maindata$Pseudotimepattern)) {
                              if (typeof(Maindata$Pseudotimepattern[[name]]) == "list") {
                                    tmppattern <- rep(name,length(row.names(Maindata$Pseudotimepattern[[name]])))
                                    names(tmppattern) <- row.names(Maindata$Pseudotimepattern[[name]])
                                    tmp <- c(tmp,tmppattern)
                              } else {
                                    tmppattern <- rep(name,length(names(Maindata$Pseudotimepattern[[name]])))
                                    names(tmppattern) <- names(Maindata$Pseudotimepattern[[name]])
                                    tmp <- c(tmp,tmppattern)
                              }
                        }
                        Maindata$Pseudotimepatternsimple <- tmp
                  })
      })
      
      output$Pseudotimepatternsummary <- renderDataTable({
            if (!is.null(Maindata$Pseudotimepatternsimple)) {
                  tmp <- data.frame(table(Maindata$Pseudotimepatternsimple))
                  colnames(tmp) <- c("Pattern","Number")                  
                  tmp      
            }            
      })
      
      output$Pseudotimepatternsummarysave <- downloadHandler(
            filename = function() { "Pattern Summary.csv" },
            content = function(file) {                  
                  if (!is.null(Maindata$Pseudotimepatternsimple)) {
                        tmp <- data.frame(table(Maindata$Pseudotimepatternsimple))
                        colnames(tmp) <- c("Pattern","Number")                  
                        write.csv(tmp,file=file,quote=F,row.names=F)                  
                  }                     
            }
      )
      
      output$PseudotimepatternsummaryGeneselectui <- renderUI(
            if (!is.null(Maindata$Pseudotimepattern))
                  tagList(
                        checkboxInput("PseudotimepatternsummaryGeneselectlistall","List Genes for all Pattern",value=T),
                        conditionalPanel(condition = "input.PseudotimepatternsummaryGeneselectlistall==0",
                                         selectInput("PseudotimepatternsummaryGeneselect","Select Pattern",names(Maindata$Pseudotimepattern),multiple = T)
                        )
                  )   
      )
      
      output$PseudotimepatternsummaryGene <- renderDataTable({
            if (!is.null(input$PseudotimepatternsummaryGeneselectlistall)) {
                  if (input$PseudotimepatternsummaryGeneselectlistall) {
                        tmp <- Maindata$Pseudotimepatternsimple
                        data.frame(Gene=names(tmp),Pattern=tmp) 
                  } else if (sum(!grepl("_",input$PseudotimepatternsummaryGeneselect)) > 0) {
                        tmp <- Maindata$Pseudotimepatternsimple[Maindata$Pseudotimepatternsimple %in% input$PseudotimepatternsummaryGeneselect]
                        data.frame(Gene=names(tmp),Pattern=tmp)      
                  } else {
                        allres <- NULL
                        for (i in input$PseudotimepatternsummaryGeneselect) {
                              tmp <- cbind(row.names(Maindata$Pseudotimepattern[[i]]),i,Maindata$Pseudotimepattern[[i]])
                              colnames(tmp) <- c("Gene","Pattern","Transition","LCI","UCI")
                              allres <- rbind(allres,tmp)                        
                        }
                        allres
                  }            
            }
      })      
      
      output$PseudotimepatternsummaryGenesave <- downloadHandler(
            filename = function() { "Gene List.csv" },
            content = function(file) {
                  if (input$PseudotimepatternsummaryGeneselectlistall) {
                        tmp <- Maindata$Pseudotimepatternsimple
                        tmp <- data.frame(Gene=names(tmp),Pattern=tmp)      
                  } else if (sum(!grepl("_",input$PseudotimepatternsummaryGeneselect)) > 0) {
                        tmp <- Maindata$Pseudotimepatternsimple[Maindata$Pseudotimepatternsimple %in% input$PseudotimepatternsummaryGeneselect]
                        tmp <- data.frame(Gene=names(tmp),Pattern=tmp)      
                  } else {
                        allres <- NULL
                        for (i in input$PseudotimepatternsummaryGeneselect) {
                              tmp <- cbind(row.names(Maindata$Pseudotimepattern[[i]]),i,Maindata$Pseudotimepattern[[i]])
                              colnames(tmp) <- c("Gene","Pattern","Transition","LCI","UCI")
                              allres <- rbind(allres,tmp)                        
                        }
                        tmp <- allres                        
                  }                
                  write.csv(tmp,file=file,quote=F,row.names=F)                  
            }
      )
      
      observe({
            tmp <- Maindata$Pseudotimepatternsimple
            Maindata$Pseudotimenoconstantname <- names(tmp)[grepl("up",tmp) | grepl("down",tmp)]
      })
      
      output$Pseudotimeselectgeneui <- renderUI(            
            selectizeInput("Pseudotimeselectgene","Select genes of interest",Maindata$Pseudotimenoconstantname,multiple = TRUE)
      )
      
      Pseudotimeplotfunc <- function(pattern,expr,fitexpr,pseudotime,gap,gene,showtrue=F) {            
            findgene <- function(gene) {
                  findres <- sapply(names(pattern), function(pat) {
                        if (grepl("_",pat)) {
                              gene %in% row.names(pattern[[pat]])
                        } else {
                              gene %in% names(pattern[[pat]])
                        }
                  })
                  names(which(findres))
            }
            genename <- gene
            x <- sort(pseudotime)
            gene <- sapply(gene, findgene)      
            if (length(gene) == 1) {
                  if (grepl("_",gene)) {
                        transinfo <- pattern[[gene]][names(gene),]
                        transmean <- transinfo[1,1]                  
                  } else {
                        transinfo <- pattern[[gene]][names(gene)]                  
                  }            
                  y <- expr[names(gene),names(x)]              
                  fity <- fitexpr[names(gene),names(x)]
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
                        geneexpr <- fitexpr[names(gene),names(x)]
                  } else {
                        geneexpr <- expr[names(gene),names(x)]
                  }            
                  colnames(geneexpr) <- NULL
                  geneexpr <- geneexpr[rev(row.names(geneexpr)),]                  
                  linesegmat <- NULL
                  for (g in 1:length(gene)) {
                        if (grepl("_",gene[g])) {
                              tmp <- unlist(pattern[[gene[g]]][names(gene[g]),])
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
                              legend.position = "right"                              
                              )
                  
            }      
      }
      
      output$Pseudotimevisualize <- renderPlot({     
            if (!is.null(Maindata$Pseudotimedata) && !is.null(input$Pseudotimeselectgene)) {
                  pseudotime <- Maindata$Pseudotimedata                                    
                  Pseudotimeplotfunc(Maindata$Pseudotimepattern,Maindata$Pseudotimeexpr,Maindata$Pseudotimefitexpr,pseudotime,as.numeric(input$Pseudotimegapvalue),input$Pseudotimeselectgene,showtrue=input$Pseudotimegeneexpressionshowtrue)            
            }
      })
      
      output$Pseudotimevisualizeui <- renderUI({
            plotOutput("Pseudotimevisualize",width=700,height=ifelse(length(input$Pseudotimeselectgene)==1,700,400))
      })
      
      output$Pseudotimevisualizesave <- downloadHandler(
            filename = function() { "Pattern.pdf" },
            content = function(file) {
                  if (!is.null(Maindata$Pseudotimedata) && !is.null(input$Pseudotimeselectgene)) {
                        pdf(file)
                        pseudotime <- Maindata$Pseudotimedata                                                
                        Pseudotimeplotfunc(Maindata$Pseudotimepattern,Maindata$Pseudotimeexpr,Maindata$Pseudotimefitexpr,pseudotime,as.numeric(input$Pseudotimegapvalue),input$Pseudotimeselectgene,showtrue=F)            
                        dev.off()
                  }
            }
      )
      
      output$Pseudotimevisualizetable <- renderDataTable({
            if (!is.null(input$Pseudotimeselectgene)) {
                  data.frame(Gene=input$Pseudotimeselectgene,Pattern=Maindata$Pseudotimepatternsimple[input$Pseudotimeselectgene])
            }
      })
      
      output$Pseudotimevisualizetablesave <- downloadHandler(
            filename = function() { "Gene Pattern.csv" },
            content = function(file) {                  
                  tmp <- data.frame(Gene=input$Pseudotimeselectgene,Pattern=Maindata$Pseudotimepatternsimple[input$Pseudotimeselectgene])
                  write.csv(tmp,file=file,quote=F,row.names=F)                                                      
            }
      )
      
      output$PseudotimeGOanalysisselectui <- renderUI({
            if (!is.null(Maindata$Pseudotimepatternsimple))
                  selectInput("PseudotimeGOanalysisselect","Select pattern of interest",sort(unique(Maindata$Pseudotimepatternsimple)),multiple = T)
      })
      
      observe({
            if (input$PseudotimeGOanalysisrunbutton > 0) {
                  isolate({
                        if (!is.null(input$PseudotimeGOanalysisselect) && Maindata$genenumwithid > 0) {
                              withProgress(message = 'Performing GO analysis...', value = 0, {
                                    if (input$Readexprchoosespecies=="Human") {                        
                                          mapdb <- "org.Hs.eg.db"
                                    } else {                        
                                          mapdb <- "org.Mm.eg.db"
                                    }                  
                                    allgene <- row.names(Maindata$rawexpr)
                                    if (input$Readexprchooseanno == "ensembl") {
                                          allgene <- sapply(allgene,function(i) strsplit(i,"\\.")[[1]][1])
                                    }                   
                                    inputgene <- allgene[Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisselect]                              
                                    geneList <- factor(as.integer(allgene %in% inputgene))
                                    names(geneList) <- allgene
                                    Maindata$PseudotimeGOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = mapdb, ID = input$Readexprchooseanno)
                                    Maindata$PseudotimeGOresultFisher <- runTest(Maindata$PseudotimeGOdata, algorithm = "classic", statistic = "fisher")
                              })
                        }                        
                  })
            }            
      })
      
      output$PseudotimeGOanalysisresult <- renderDataTable({
            if (!is.null(Maindata$PseudotimeGOresultFisher)) {
                  tmp <- GenTable(Maindata$PseudotimeGOdata, classicFisher = Maindata$PseudotimeGOresultFisher, topNodes = as.numeric(input$PseudotimeGOanalysistermnum),orderBy="classicFisher")                  
                  tmp
            }
      })
      
      output$PseudotimeGOanalysisresultsave <- downloadHandler(
            filename = function() { "GO Analysis.csv" },
            content = function(file) {
                  if (!is.null(Maindata$PseudotimeGOresultFisher)) {
                        tmp <- GenTable(Maindata$PseudotimeGOdata, classicFisher = Maindata$PseudotimeGOresultFisher, topNodes = as.numeric(input$PseudotimeGOanalysistermnum),orderBy="classicFisher")                  
                        tmp$Term <- gsub(","," ",tmp$Term)
                        write.csv(tmp,file=file,quote=F,row.names=F)    
                  }
            }
      )
      
      output$PseudotimeGOanalysisTimeselectui <- renderUI({
            if (!is.null(Maindata$Pseudotimepatternsimple)) {
                  tmp <- sort(unique(Maindata$Pseudotimepatternsimple))
                  tmp <- tmp[grep("_",tmp)]
                  if (length(tmp) == 0) {
                        strong("There is no pattern with transition point!")
                  } else {                        
                        selectInput("PseudotimeGOanalysisTimeselect","Select pattern of interest",tmp,multiple = T)                                                
                  }                  
            }                  
      })
      
      output$PseudotimeGOanalysisTimewinsizesliderui <- renderUI({            
            if (!is.null(input$PseudotimeGOanalysisTimeselect)) {                  
                  tagList(
                        helpText("Select window size (how many genes to be included in each GO analsysis). Notice that the results of GO analyses are not robust if too few genes and included."),
                        sliderInput("PseudotimeGOanalysisTimewinsizeslider","Select window size",min=1,max=sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect),step=1,value=ifelse(sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect) > 50,50,ceiling(sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect)/2)))     
                  )                  
            }                  
      })
      
      output$PseudotimeGOanalysisTimestepsliderui <- renderUI({
            if (!is.null(input$PseudotimeGOanalysisTimewinsizeslider)) 
                  sliderInput("PseudotimeGOanalysisTimestepslider","Select moving step",min=1,max=sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect)-1,step=1,value=ifelse(sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect)-1 > 25,25,ceiling((sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect)-1)/2)))
      })      
      
      observe({
            if (input$PseudotimeGOanalysisTimerunbutton > 0) {
                  isolate({
                        if (!is.null(input$PseudotimeGOanalysisTimeselect) && Maindata$genenumwithid > 0) {
                              GOanalysis <- function(inputgene,allgene,identifier="ENSEMBL",species="Human") {                                        
                                    if (species == "Human") {
                                          mapdb <- "org.Hs.eg.db"      
                                    } else if (species == "Mouse") {
                                          mapdb <- "org.Mm.eg.db"      
                                    }                  
                                    geneList <- factor(as.integer(allgene %in% inputgene))
                                    names(geneList) <- allgene
                                    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = mapdb, ID = identifier)
                                    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")                                    
                                    list(resultFisher=resultFisher,GOdata=GOdata)
                              }
                              allres <- NULL
                              for (i in input$PseudotimeGOanalysisTimeselect) {
                                    tmp <- data.frame(row.names(Maindata$Pseudotimepattern[[i]]),i,Maindata$Pseudotimepattern[[i]],stringsAsFactors = F)
                                    colnames(tmp) <- c("Gene","Pattern","Transition","LCI","UCI")
                                    allres <- rbind(allres,tmp)                        
                              }                              
                              genelist <- allres[order(allres$Transition),1]
                              if (input$Readexprchooseanno!="symbol" && input$Readexpraddsymbol)
                                    genelist <- sapply(genelist,function(i) strsplit(i,"_")[[1]][2])
                              allgene <- row.names(Maindata$rawexpr)
                              if (input$Readexprchooseanno == "ensembl") {
                                    allgene <- sapply(allgene,function(i) strsplit(i,"\\.")[[1]][1])
                                    genelist <- sapply(genelist,function(i) strsplit(i,"\\.")[[1]][1])                                    
                              } 
                              withProgress(message = 'Calculation in Progress...', value = 0, {
                                    genenum <- sum(Maindata$Pseudotimepatternsimple %in% input$PseudotimeGOanalysisTimeselect)
                                    windowsize <- as.numeric(input$PseudotimeGOanalysisTimewinsizeslider)
                                    movesize <- as.numeric(input$PseudotimeGOanalysisTimestepslider)
                                    maxwinnum <- ceiling((genenum-windowsize)/movesize) + 1                                    
                                    res <- list()
                                    for (i in 1:maxwinnum) {
                                          if (i==maxwinnum) {
                                                startid <- genenum-windowsize+1
                                                endid <- genenum
                                          } else {
                                                startid <- 1+movesize*(i-1)
                                                endid <- movesize*(i-1)+windowsize
                                          }                                               
                                          res[[paste0(startid,"-",endid)]] <- GOanalysis(genelist[startid:endid],allgene,identifier=input$Readexprchooseanno,species=input$Readexprchoosespecies)      
                                          incProgress(1/maxwinnum, detail = paste("Performing GO analysis for Window", i))   
                                    }
                                    Maindata$PseudotimewindowGOres <- res                                                                                                            
                              })
                        }                        
                  })
            }            
      })
      
      observe({
            res <- list()
            for (name in names(Maindata$PseudotimewindowGOres)) {
                  res[[name]] <- GenTable(Maindata$PseudotimewindowGOres[[name]][["GOdata"]], classicFisher = Maindata$PseudotimewindowGOres[[name]][["resultFisher"]], topNodes = as.numeric(input$PseudotimeGOanalysisTimetermnum),orderBy="classicFisher")                                          
            }
            Maindata$PseudotimewindowGOtable <- res            
      })
      
      output$PseudotimeGOanalysisTimeselectresultwindowui <- renderUI({
            selectInput("PseudotimeGOanalysisTimeselectresultwindow","Select window",names(Maindata$PseudotimewindowGOtable))
      })
      
      output$PseudotimeGOanalysisTimeshowresultwindow <- renderDataTable({
            if (!is.null(input$PseudotimeGOanalysisTimeselectresultwindow))
                  Maindata$PseudotimewindowGOtable[[input$PseudotimeGOanalysisTimeselectresultwindow]]
      })
      
      output$PseudotimeGOanalysisTimeshowresultwindowsave <- downloadHandler(
            filename = function() { "GO Analysis.csv" },
            content = function(file) {
                  if (!is.null(input$PseudotimeGOanalysisTimeselectresultwindow))                        
                        write.csv(Maindata$PseudotimewindowGOtable[[input$PseudotimeGOanalysisTimeselectresultwindow]],file=file,quote=F,row.names=F)    
            }            
      )
      
      output$PseudotimeGOanalysisTimeselectresulttimetrendselectmethodui <- renderUI({
            if (input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=='Specific') {
                  selectInput("PseudotimeGOanalysisTimeselectresulttimetrendspecific","Select GO term",unique(as.vector(sapply(Maindata$PseudotimewindowGOtable,function(i) i[,1]))),multiple = T)
            } else {
                  sliderInput("PseudotimeGOanalysisTimeselectresulttimetrendtop","Select Number of Top GO Terms",1,as.numeric(input$PseudotimeGOanalysisTimetermnum),1,1)
            }
      })
      
      output$PseudotimeGOanalysisTimeshowresulttimetrend <- renderPlot({
            if (!is.null(Maindata$PseudotimewindowGOtable) && ((input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=="Specific" && !is.null(input$PseudotimeGOanalysisTimeselectresulttimetrendspecific)) || (input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=="Top" && !is.null(input$PseudotimeGOanalysisTimeselectresulttimetrendtop)))) {
                  GOres <- Maindata$PseudotimewindowGOtable
                  termnum <- as.numeric(input$PseudotimeGOanalysisTimetermnum)      
                  if (input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=="Specific") {
                        GOterm <- input$PseudotimeGOanalysisTimeselectresulttimetrendspecific
                  } else {
                        GOterm <- unique(as.vector(sapply(GOres,function(i) i[1:as.numeric(input$PseudotimeGOanalysisTimeselectresulttimetrendtop),1])))
                  }
                  GOdes <- NULL
                  for (i in GOres) {
                        GOdes <- rbind(GOdes,as.matrix(i[,1:2]))
                  }
                  GOdes <- unique(GOdes)
                  rankres <- sapply(GOterm,function(term) {
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
                  if (input$PseudotimeGOanalysisTimeselectresulttimetrendshowheatmap) {
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
      })
      
output$PseudotimeGOanalysisTimeshowresultwindowplotsave <- downloadHandler(
      filename = function() { "GO pattern.pdf" },
      content = function(file) {
            if (!is.null(Maindata$PseudotimewindowGOtable) && ((input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=="Specific" && !is.null(input$PseudotimeGOanalysisTimeselectresulttimetrendspecific)) || (input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=="Top" && !is.null(input$PseudotimeGOanalysisTimeselectresulttimetrendtop)))) {
                  pdf(file,width=14,height=7)
                  GOres <- Maindata$PseudotimewindowGOtable
                  termnum <- as.numeric(input$PseudotimeGOanalysisTimetermnum)      
                  if (input$PseudotimeGOanalysisTimeselectresulttimetrendselectmethod=="Specific") {
                        GOterm <- input$PseudotimeGOanalysisTimeselectresulttimetrendspecific
                  } else {
                        GOterm <- unique(as.vector(sapply(GOres,function(i) i[1:as.numeric(input$PseudotimeGOanalysisTimeselectresulttimetrendtop),1])))
                  }
                  GOdes <- NULL
                  for (i in GOres) {
                        GOdes <- rbind(GOdes,as.matrix(i[,1:2]))
                  }
                  GOdes <- unique(GOdes)
                  rankres <- sapply(GOterm,function(term) {
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
                  if (input$PseudotimeGOanalysisTimeselectresulttimetrendshowheatmap) {
                        p <- ggplot(data=rankres, aes(x=Var1, y=GOTerm)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient2(low = "blue",high = "red",mid="white")
                  } else {
                        p <- ggplot(data = rankres, aes(x=Var1, y=value, colour=GOTerm)) +
                              geom_line(aes(group=GOTerm)) + tmpyset +
                              geom_point(size=4)
                  }
                  p <- p + xlab("Interval") +
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
                  print(p)
            }
                  dev.off()
            }     
)      


      ### Comparison ###
      
      observe({
            if (is.null(Maindata$Comparisonrawdata))
                  Maindata$Comparisonrawdata <- list()
      })
      
      output$Comparisonlistnameui <- renderUI({
            if (!is.null(Maindata$Comparisonrawdata))
                  textInput("Comparisonlistname","List Name",paste0("List",length(Maindata$Comparisonrawdata) + 1))
      })
      
      observe({
            if (input$Comparisonreadin > 0)
                  isolate({
                        FileHandle <- input$ComparisonFile
                        if (!is.null(FileHandle)) {
                              Maindata$Comparisonrawdata[[input$Comparisonlistname]] <- read.csv(FileHandle$datapath,header=T,as.is=T,blank.lines.skip=TRUE)
                        }
                  })
      })
      
      observe({
            if (!input$Comparisonsimplify && !input$Comparisonignoreconst) {
                  Maindata$Comparisondata <- Maindata$Comparisonrawdata
            } else {
                  Maindata$Comparisondata <- lapply(Maindata$Comparisonrawdata, function(i) {
                        i[,2] <- sapply(i[,2], function(pattern) {
                              pattern <- strsplit(pattern,"_")[[1]]
                              if (input$Comparisonignoreconst) 
                                    pattern <- pattern[pattern!="constant"]
                              if (input$Comparisonsimplify)
                                    pattern <- rle(pattern)$values    
                              if (length(pattern) == 0) {
                                    "NA"
                              } else {
                                    paste(pattern,collapse = "_")      
                              }                              
                        })  
                        i
                  })          
            }
            
      })
      
      output$Comparisonchoosefileui <- renderUI({
            selectInput("Comparisonchoosefile","Choose file to display",names(Maindata$Comparisonrawdata),names(Maindata$Comparisonrawdata)[length(Maindata$Comparisonrawdata)])
      })
      
      output$Comparisonshowfile <- renderDataTable({
            if (!is.null(input$Comparisonchoosefile))
                  Maindata$Comparisondata[[input$Comparisonchoosefile]]
      })
      
      observe({
            if (length(Maindata$Comparisondata) > 1) {
                  tmp <- Reduce(function(x, y) merge(x, y, by="Gene"), Maindata$Comparisondata)
                  colnames(tmp) <- c("Gene",names(Maindata$Comparisondata))
                  Maindata$Comparisonsumdata <- tmp      
            }
      })
      
      output$Comparisonshowcomparisonresults <- renderDataTable({
            if (!is.null(Maindata$Comparisonsumdata)) {
                  if (input$Comparisononlyshowdifferentgene=="diffpattern") {
                        Maindata$Comparisonsumdata[apply(Maindata$Comparisonsumdata[,-1],1,function(x) length(unique(x))) > 1,]
                  } else if (input$Comparisononlyshowdifferentgene=="samepattern") {
                        Maindata$Comparisonsumdata[apply(Maindata$Comparisonsumdata[,-1],1,function(x) length(unique(x))) == 1,]
                  } else {
                        Maindata$Comparisonsumdata
                  }      
            }                              
      })
      
      
      
      
      
})


