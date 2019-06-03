#
#
####### minfiShinyl App #######
##
##UI.R code


# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


#Loading libraries
library(shiny)
library(minfi)
library(minfiData)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggfortify)
library(Hmisc)
library(ComplexHeatmap)
library(circlize) 
library(DT)
library(dplyr)
library("IlluminaHumanMethylationEPICmanifest")
library(pvclust)



# Define UI for application 
shinyUI(pageWithSidebar(

  
  # Sidebar with a slider input for number of bins 
  headerPanel("minfiShinyl"),
  sidebarPanel(
    fileInput('targetsFile', 'Choose your targets file',
              accept=c('csv', 'comma-separated-values', '.csv')),
    p(strong("Input"), ": targets file", br(), 
      "refers to a CSV file with 12 datacamps:"),
    em('Sample_Name,  Sample_Well, Sample_Plate, Sample_Group, 
       Pool_ID, person_id, the covariables, status, Array, Slide, 
       Basename, appPath'),
    tags$hr(),
    p("Take a look to", a(href="targets_example.html", target="_blank",
                          "targets file"), "input example"),
    checkboxInput('header', 'Header', TRUE),
    radioButtons('sep', 'Separator', 
                 c(Comma=",",
                   Semicolon=";",
                   Tab="\t"),
                 ","
                 ),
    radioButtons('quote', 'Quote',
                 c(None='',
                   'Double Quote' = '"',
                   'Single Quote' = "'"),
                 ""
    )
    ),
  
  
  # Show the outputs in the main panel
  mainPanel(
    em("minfiShinyl: a Bioinformatics and Biostatistics Master Project, June 2019"),
    br(),
    br(),
    img(src="minfiShinyl (3).png", height=120, width=90),
    img(src="UOC.jpg", height= 70, width=200),
    img(src="RStudio.jpg", height=40, width=150),
    br(),
    h3("minfiShinyl", align="center"),
    h4("The interactive methyation analysis", em("Shiny App"), 
       "interface for", em("HumanMethylationEPIC arrays"), 
       align="center"),
    code("Click to see the", 
         a("minfi", href= "https://bioconductor.org/packages/release/bioc/html/minfi.html/",
           target="_blank"), "R Bioconductor package"),
    br(),
    code("Visit the", a("Shiny", href="https://shiny.rstudio.com/", 
                        target="_blank"), "RStudio product"),
    br(),
    br(),
    br(),
    br(),
    tabsetPanel(type="tabs",
                tabPanel("Targets Table", dataTableOutput('contents')),
                tabPanel("QC Samples Plot", plotOutput('QCplot'), 
                         downloadButton('downloadPlot', 'Download Plot'),
                         textOutput('warningBadSample')),
                tabPanel("Beta values Density Plot", plotOutput('DensPlot'),
                         downloadButton('downloadDensPlot', 
                                        'Download Density Plot')),
                tabPanel("PCA plot Distribution", plotOutput('PCAPlot'),
                         downloadButton('downloadPCAPlot',
                                        'Download PCA distribution Plot')),
                tabPanel("Heatmap", plotOutput('Heatmap'),
                         downloadButton('downloadHeatmap',
                                        'Download Heatmap')),
                tabPanel("Cluster", plotOutput('cluster'),
                         downloadButton('downloadCluster',
                                        'Download cluster')),
                tabPanel("Beta values", tableOutput('betas'),
                         downloadButton('downloadBetas',
                         'Download beta values')),
                tabPanel("M values", tableOutput('M'),
                         downloadButton('downloadM',
                         'Download M values')),
                tabPanel("DMR Table", tableOutput('DMRtab'),
                         downloadButton('downloadDMRtable',
                         'Download DMRs table')),
                tabPanel("DMR Coef",tableOutput('DMRcoef'),
                         downloadButton('downloadDMRcoef',
                          'Download DMRs coef')),
                tabPanel("dmp", tableOutput('dmp'),
                         downloadButton('downloaddmp',
                         'Download dmps')),
                tabPanel("Annotation", tableOutput('annotation'),
                         downloadButton('downloadAnnotation',
                                        'Download annotation data'))
    ),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br()
    
  )
  ))


