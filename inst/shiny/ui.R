######################################################
##                       SEPA                       ##
##             Interactive User Interface           ##
##                     UI File                      ##
##           Author:Zhicheng Ji, Hongkai Ji         ##
##       Maintainer:Zhicheng Ji (zji4@jhu.edu)      ##
######################################################

shinyUI(      
      pageWithSidebar(
            
            headerPanel('SEPA: Single-Cell Gene Expression Pattern Analysis'),
            
            sidebarPanel(
                  
                  wellPanel(
                        helpText(a("Youtube short video demo",href="",target="_blank")),
                        radioButtons("MainMenu","Main Menu",
                                     list("Reading in Gene Expression Profile"="Readexp",
                                          "Analysis for True Experimental Time"="Truetime",
                                          "Analysis for Pseudo Temporal Cell Ordering"="Pseudotime",
                                          "Pattern Comparison"="Comparison",
                                          "About"="About")
                        )
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Readexp'",
                                   radioButtons("Readexpselectstep","Select Step",list("Step 1: Input Gene Expression Data"="S1","Step 2: Choose species and identifier"="S2","Step 3: Take logarithm (optional)"="S3","Step 4: Filter genes and cells (optional)"="S4")),
                                   wellPanel(
                                         conditionalPanel(condition="input.Readexpselectstep=='S1'",
                                                          h4("Step 1: Input Gene Expression Data (See Instructions on the right!)"),
                                                          fileInput('ReadexpFile', 'Choose File'),
                                                          p(actionButton("Readexpreadin","Read in")),
                                                          radioButtons('Readexpsep', 'Separator',
                                                                       c('Tab'='\t',
                                                                         'Space'=' ',
                                                                         'Comma(csv)'=',',
                                                                         'Semicolon'=';'
                                                                       ),
                                                                       '\t')
                                         ),
                                         conditionalPanel(condition="input.Readexpselectstep=='S2'",
                                                          h4("Step 2: Choose species and identifier"),
                                                          strong("Important: Make sure the species and gene identifiers are correctly specified. Otherwise GO enrichement analysis cannot be performed."),
                                                          radioButtons("Readexprchoosespecies","Choose Species",c("Human","Mouse")),
                                                          selectizeInput("Readexprchooseanno","Choose Gene Identifier", c("entrez", "ensembl", "symbol", "unigene")),
                                                          uiOutput("Readexpraddsymbolui")                                                          
                                         ),
                                         conditionalPanel(condition="input.Readexpselectstep=='S3'",
                                                          h4("Step 3: Take logarithm (optional)"),
                                                          checkboxInput("Readexplogtf","Take log of current data",value = T),
                                                          conditionalPanel(condition="input.Readexplogtf==1",                                         
                                                                           radioButtons("Readexplogbase","Choose log base",choices = c("2","10","e")),
                                                                           textInput("Readexplogpseudocount","Enter pseudo count added when taking log",value = 1)
                                                          )           
                                         ),
                                         conditionalPanel(condition="input.Readexpselectstep=='S4'",
                                                          h4("Step 4: Filter genes and cells (optional)"),
                                                          h4("Filter Cells"),
                                                          uiOutput("Readexpcellui"),
                                                          h4("Filter Genes"),
                                                          h5("Genes with constant zero expressions across all samples are automatically removed"),
                                                          h5("1. Keep genes with expressions larger than X in more than Y percent of all cells"),
                                                          textInput("Readexpgeneexpcutoff","X: Expression cutoff",0),
                                                          textInput("Readexpgenepercentcutoff","Y: Percent cutoff",0),
                                                          h5("2. Keep genes with coefficient of variation (CV, sd/mean) larger than the cutoff"),
                                                          textInput("Readexpgenecvcutoff","CV cutoff",0)                                                                                             
                                         )
                                   )
                  ),
                  conditionalPanel(condition="input.MainMenu=='Truetime'",
                                   radioButtons("Truetimeselectstep","Select Step",list("Step 1: Input Time Information"="S1","Step 2: Identify Pattern"="S2","Step 3: Visualize Pattern (optional)"="S3","Step 4: GO Analysis (optional)"="S4")),
                                   wellPanel(
                                         conditionalPanel(condition="input.Truetimeselectstep=='S1'",
                                                          h4("Step 1: Input True Experiment Time Information (See Instructions on the right!)"),
                                                          fileInput('TruetimeFile', 'Choose File'),
                                                          p(actionButton("Truetimereadin","Read in")),
                                                          radioButtons('Truetimesep', 'Separator',
                                                                       c('Tab'='\t',
                                                                         'Space'=' ',
                                                                         'Comma(csv)'=',',
                                                                         'Semicolon'=';'
                                                                       ),
                                                                       '\t')
                                         ),
                                         conditionalPanel(condition="input.Truetimeselectstep=='S2'",
                                                          h4("Step 2: Identify Pattern"),
                                                          helpText("The calculation could take a long time depending on the size of the dataset."),
                                                          textInput("Truetimepvalcutoff","Choose cutoff for adjusted p-values of t-tests",0.05),                                                          
                                                          actionButton("Truetimerunanalysis","Run Analysis")                                                          
                                         ),
                                         conditionalPanel(condition="input.Truetimeselectstep=='S3'",
                                                          h4("Step 3: Visualize Gene Expression Pattern (optional)"),
                                                          uiOutput("Truetimeselectgeneui")
                                         ),
                                         conditionalPanel(condition="input.Truetimeselectstep=='S4'",
                                                          h4("Step 4: GO Enrichment Analysis (optional)"),
                                                          helpText("GO analysis could take a long time for a large gene set."),
                                                          uiOutput("TruetimeGOanalysisselectui"),
                                                          actionButton("TruetimeGOanalysisrunbutton","Run GO analysis"),
                                                          sliderInput("TruetimeGOanalysistermnum","Select number of terms to display",1,200,20,1)
                                         )
                                         
                                   ),                                   
                                   wellPanel(
                                         h4("Other Options"),
                                         checkboxInput("Truetimesimplify","Simplify expression patterns. (up_up_constant will be simplified to up_constant)",value=TRUE),
                                         checkboxInput("Truetimeignoreconst","Ignore constant patterns. (up_up_constant will be converted to up_up)",value=FALSE)
                                   )
                  ),
                  
                  
                  conditionalPanel(condition="input.MainMenu=='Pseudotime'",
                                   radioButtons("Pseudotimeselectstep","Select Step",list("Step 1: Input Time Information"="S1","Step 2: Identify Pattern"="S2","Step 3: Visualize Pattern (optional)"="S3","Step 4: GO Analysis (optional)"="S4")),
                                   wellPanel(
                                         conditionalPanel(condition="input.Pseudotimeselectstep=='S1'",
                                                          h4("Step 1: Input True Experiment Time Information (See Instructions on the right!)"),
                                                          fileInput('PseudotimeFile', 'Choose File'),
                                                          p(actionButton("Pseudotimereadin","Read in")),
                                                          radioButtons('Pseudotimesep', 'Separator',
                                                                       c('Tab'='\t',
                                                                         'Space'=' ',
                                                                         'Comma(csv)'=',',
                                                                         'Semicolon'=';'
                                                                       ),
                                                                       '\t')
                                         ),
                                         conditionalPanel(condition="input.Pseudotimeselectstep=='S2'",
                                                          h4("Step 2: Identify Pattern"),
                                                          helpText("The calculation could take a long time depending on the size of the dataset."),
                                                          textInput("Pseudotimepvalcutoff","Choose cutoff for adjusted p-values of t-tests",0.05),  
                                                          textInput("Pseudotimegapvalue","Choose number of gap cells",10),
                                                          checkboxInput("Pseudotimesimplify","Simplify expression patterns. (up_up_constant will be simplified to up_constant)",value=TRUE),
                                                          checkboxInput("Pseudotimeignoreconst","Ignore constant patterns. (up_up_constant will be converted to up_up)",value=FALSE),
                                                          actionButton("Pseudotimerunanalysis","Run Analysis")                                                          
                                         ),
                                         conditionalPanel(condition="input.Pseudotimeselectstep=='S3'",
                                                          h4("Step 3: Visualize Gene Expression Pattern (optional)"),
                                                          uiOutput("Pseudotimeselectgeneui"),
                                                          checkboxInput("Pseudotimegeneexpressionshowtrue","Show true expression",value=F)
                                         ),
                                         conditionalPanel(condition="input.Pseudotimeselectstep=='S4'",
                                                          h4("Step 4: GO Enrichment Analysis (optional)"),
                                                          radioButtons("PseudotimeGOanalysischoosetype","Choose Analysis Type",list("Gene Expression Pattern"="Pattern","Transition Point Time Window"="Time")),
                                                          conditionalPanel(condition = "input.PseudotimeGOanalysischoosetype=='Pattern'",                                                                           
                                                                           helpText("Perform GO analysis on all genes for certain patterns. It could take a long time for a large gene set."),
                                                                           uiOutput("PseudotimeGOanalysisselectui"),
                                                                           actionButton("PseudotimeGOanalysisrunbutton","Run GO analysis"),
                                                                           sliderInput("PseudotimeGOanalysistermnum","Select number of terms to display",1,200,20,1)                                                                          
                                                          ),
                                                          conditionalPanel(condition = "input.PseudotimeGOanalysischoosetype=='Time'",                                                                           
                                                                           helpText("This analysis type is specifically designed for genes with transition points. Genes are ordered according to the transition points and GO analyses are performed repeatedly on genes within certain quantiles. It could take a long time for a large gene set and small window size."),
                                                                           uiOutput("PseudotimeGOanalysisTimeselectui"),
                                                                           uiOutput("PseudotimeGOanalysisTimewinsizesliderui"),
                                                                           uiOutput("PseudotimeGOanalysisTimestepsliderui"),
                                                                           actionButton("PseudotimeGOanalysisTimerunbutton","Run GO analysis"),
                                                                           sliderInput("PseudotimeGOanalysisTimetermnum","Select number of terms to display",1,200,20,1)                 
                                                                           
                                                          )                                                          
                                         )
                                         
                                   )
                  ),
                  conditionalPanel(condition="input.MainMenu=='Comparison'",                                
                                                          h4("Read in pattern files generated by SEPA. The files should be in csv format."),
                                                          fileInput('ComparisonFile', 'Choose File'),
                                                          uiOutput("Comparisonlistnameui"),
                                                          p(actionButton("Comparisonreadin","Read in")),                                                          
                                                          checkboxInput("Comparisonsimplify","Simplify expression patterns for all pattern files. (up_up_constant will be simplified to up_constant)",value=TRUE),
                                                          checkboxInput("Comparisonignoreconst","Ignore constant patterns for all pattern files. (up_up_constant will be converted to up_up)",value=FALSE)
                  )                                    
            ),
            mainPanel(
                  
                  uiOutput("showbusybar"),
                  
                  conditionalPanel(condition="input.MainMenu=='Readexp'",
                                   conditionalPanel("input.Readexpselectstep=='S1'",
                                                    checkboxInput("Readexpshowinstructiontf","Show instructions",value=T)
                                   ),
                                   conditionalPanel(condition="input.Readexpshowinstructiontf==1 && input.Readexpselectstep=='S1'",
                                                    wellPanel(
                                                          h5("Instructions:"),
                                                          p("Single cell data should be prepared in a matrix-like data format. Each row corresponds to a gene/feature and each column corresponds to a single cell."),
                                                          p("Notice that each row should have the same number of entries, especially for the header (first row, see example below)"),
                                                          p("Remove all the single or double quote in the expression file"),
                                                          p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to read in different file formats."),
                                                          h5("A typical example of tab-delimited file format:"),
                                                          tags$pre(p("Gene\tT0_A01\tT0_B01\tT1_G02\tT1_G01\nSOX2\t0.455\t0.543\t0.000\t2.188\nPAT1\t0.231\t2.792\t1.222\t0.000"))                                                          
                                                    )
                                   ),
                                   uiOutput("Readexpshowsummaryui")
                  ),
                  conditionalPanel(condition="input.MainMenu=='Truetime'",
                                   conditionalPanel(condition="input.Truetimeselectstep=='S1'",
                                                    checkboxInput("Truetimeshowinstructiontf","Show instructions",value=T),
                                                    conditionalPanel(condition="input.Truetimeshowinstructiontf==1",
                                                                     wellPanel(
                                                                           h5("Instructions:"),
                                                                           p("True experiment time data should be prepared in a matrix-like data format."),
                                                                           p("Each row corresponds to a cell. First column is the cell name. Second column is the true experiment time"),
                                                                           p("Please make sure the cell names here are exactly the same as the cells names in the expression file."),
                                                                           p("Remove all the single or double quote in the file"),
                                                                           p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to read in different file formats."),
                                                                           h5("A typical example of tab-delimited file format:"),                                                          
                                                                           tags$pre(p("Cell\tTime\nT0_A01\tT0\nT0_B01\tT0\nT1_G01\tT1\nT1_G02\tT1"))                                                                   
                                                                     )
                                                    ),
                                                    dataTableOutput("Truetimeshowtime")
                                   ),
                                   conditionalPanel(condition="input.Truetimeselectstep=='S2'",
                                                    tabsetPanel(
                                                          tabPanel("List All Genes",wellPanel(downloadButton("TruetimepatternsummaryGenesave","Save Gene List")),uiOutput("TruetimepatternsummaryGeneselectui"),dataTableOutput("TruetimepatternsummaryGene")),
                                                          tabPanel("Summary Table",wellPanel(downloadButton("Truetimepatternsummarysave","Save Summary Table")),dataTableOutput("Truetimepatternsummary"))                                                          
                                                    )
                                   ),
                                   conditionalPanel(condition="input.Truetimeselectstep=='S3'",
                                                    wellPanel(
                                                          downloadButton("Truetimevisualizesave","Save Gene Expression Plot"),
                                                          downloadButton("Truetimevisualizetablesave","Save Gene Pattern Table")
                                                    ),
                                                    plotOutput("Truetimevisualize"),
                                                    dataTableOutput("Truetimevisualizetable")
                                   ),
                                   conditionalPanel(condition="input.Truetimeselectstep=='S4'",
                                                    wellPanel(downloadButton("TruetimeGOanalysisresultsave","Save GO Analysis Results")),
                                                    dataTableOutput("TruetimeGOanalysisresult")
                                   )                                                                                                         
                  ),
                  
                  conditionalPanel(condition="input.MainMenu=='Pseudotime'",
                                   conditionalPanel(condition="input.Pseudotimeselectstep=='S1'",
                                                    checkboxInput("Pseudotimeshowinstructiontf","Show instructions",value=T),
                                                    conditionalPanel(condition="input.Pseudotimeshowinstructiontf==1",
                                                                     wellPanel(
                                                                           h5("Instructions:"),
                                                                           p("Pseudo-time ordering data should be prepared in a matrix-like data format."),
                                                                           p("Each row corresponds to a cell. First column is the cell name. Second column is the pseudo-time."),
                                                                           p("Please make sure the cell names here are exactly the same as the cells names in the expression file."),
                                                                           p("Remove all the single or double quote in the file"),
                                                                           p("Please make sure the data is correctly read in before any further analysis is conducted. Adjust the options on the left to read in different file formats."),
                                                                           h5("A typical example of tab-delimited file format:"),                                                          
                                                                           tags$pre(p("Cell\tTime\nT0_A01\t0\nT0_B01\t23.4\nT1_G01\t68.9\nT1_G02\t100"))                                                                   
                                                                     )
                                                    ),
                                                    dataTableOutput("Pseudotimeshowtime")
                                   ),
                                   conditionalPanel(condition="input.Pseudotimeselectstep=='S2'",
                                                    tabsetPanel(                                                          
                                                          tabPanel("List All Genes",wellPanel(downloadButton("PseudotimepatternsummaryGenesave","Save Gene List")),uiOutput("PseudotimepatternsummaryGeneselectui"),dataTableOutput("PseudotimepatternsummaryGene")),
                                                          tabPanel("Summary Table",wellPanel(downloadButton("Pseudotimepatternsummarysave","Save Summary Table")),dataTableOutput("Pseudotimepatternsummary"))
                                                    )
                                   ),
                                   conditionalPanel(condition="input.Pseudotimeselectstep=='S3'",
                                                    wellPanel(
                                                          downloadButton("Pseudotimevisualizesave","Save Gene Expression Plot"),
                                                          downloadButton("Pseudotimevisualizetablesave","Save Gene Pattern Table")
                                                    ),
                                                    uiOutput("Pseudotimevisualizeui"),
                                                    dataTableOutput("Pseudotimevisualizetable")
                                   ),
                                   conditionalPanel(condition="input.Pseudotimeselectstep=='S4'",
                                                    conditionalPanel(condition = "input.PseudotimeGOanalysischoosetype=='Pattern'",                                                                           
                                                            wellPanel(downloadButton("PseudotimeGOanalysisresultsave","Save GO Analysis Results")),
                                                          dataTableOutput("PseudotimeGOanalysisresult")),
                                                    conditionalPanel(condition = "input.PseudotimeGOanalysischoosetype=='Time'", 
                                                                     tabsetPanel(
                                                                           tabPanel("GO Terms Table",
                                                                                    wellPanel(downloadButton("PseudotimeGOanalysisTimeshowresultwindowsave","Save GO Analysis Results")),                      
                                                                                    uiOutput("PseudotimeGOanalysisTimeselectresultwindowui"),
                                                                                    dataTableOutput("PseudotimeGOanalysisTimeshowresultwindow")
                                                                                    ),
                                                                           
                                                                           tabPanel("Single GO terms",
                                                                                    wellPanel(downloadButton("PseudotimeGOanalysisTimeshowresultwindowplotsave","Save Plot")),
                                                                                    radioButtons("PseudotimeGOanalysisTimeselectresulttimetrendselectmethod","",list("Specific GO terms"="Specific","Top GO terms"="Top")),
                                                                                    checkboxInput("PseudotimeGOanalysisTimeselectresulttimetrendshowheatmap","Show results with heatmap",value=F),
                                                                                    uiOutput("PseudotimeGOanalysisTimeselectresulttimetrendselectmethodui"),
                                                                                    plotOutput("PseudotimeGOanalysisTimeshowresulttimetrend")
                                                                                    )
                                                                           )
                                                    )
                                   )              
                  ),
                  conditionalPanel(condition="input.MainMenu=='Comparison'",
                                   tabsetPanel(
                                         tabPanel("Pattern Files",
                                                  uiOutput("Comparisonchoosefileui"),
                                                  dataTableOutput("Comparisonshowfile")
                                                  ),                                          
                                         tabPanel("Comparsion Results",
                                                  radioButtons("Comparisononlyshowdifferentgene","",list("Only show genes with different patterns"="diffpattern","Only show genes with same patterns"="samepattern","Show all genes"="all")),
                                                  dataTableOutput("Comparisonshowcomparisonresults")
                                                  )                                         
                                         )
                  ),
                  conditionalPanel(condition="input.MainMenu=='About'",
                                   p('SEPA: Single-Cell Gene Expression Pattern Analysis'),
                                   p('Current Version: 0.99.0'),
                                   p('Release Date: 2015-04-17'),
                                   p('Author: Zhicheng Ji,Hongkai Ji'),
                                   p('Maintainer: Zhicheng Ji <zji4@jhu.edu>'),
                                   p(a("Visit my homepage",href="http://www.biostat.jhsph.edu/~zji4/",target="_blank")),
                                   p(a("Visit web page of our lab",href="http://www.biostat.jhsph.edu/~hji/",target="_blank"))                                   
                  )
                  
                  
            )
            
      ))