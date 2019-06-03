#
#
####### minfiShinyl App #######
##
##Server code

#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#



# Loading libraries
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


# The Server Shiny Function #

shinyServer(function(input, output) {

    
  ##########Data Loading Functions 
  
  #Reading CSV file
  targetsFunc <- function(input){
    inFile <- input$targetsFile
    
    targets <- read.csv(inFile$datapath, header=input$header, sep=input$sep,
                        quote = input$quote, stringsAsFactors = F
    )
    return(targets)

  }
  
  #RGset loading
  RGsetFunc <- function(input){
    
    targets <- targetsFunc(input)
    
    RGset <- read.metharray.exp(targets=targets, force=T)
    
    return(RGset)
  }
  
  #MSet loading
  MSetFunc <- function(input){
    
    targets <- targetsFunc(input)
    
    RGset <- read.metharray.exp(targets=targets, force=T)
    
    MSet <- preprocessRaw(RGset)
    
    return(MSet)
  }
  
  #pd loading
  pdFunc <- function(input){
    
    targets <- targetsFunc(input)
    
    RGset <- read.metharray.exp(targets=targets, force=T)
    
    pd = pData(RGset)
    return(pd)
  }
  
 
  #MSet Normalizing
  MSet.QNorm.filt.Func <- function(input){
    
    targets <- targetsFunc(input)
    RGset <- read.metharray.exp(targets=targets, force=T)
    
    pd <- pData(RGset)
    
    detP <-detectionP(RGset)
    
    keep <- colMeans(detP)<0.05
    
    RGset <- RGset[,keep]
    
    targets <- targets[keep,]
    
    MSet <- preprocessRaw(RGset)
    
    MSet.QNorm <- preprocessQuantile(MSet)
    
    #Filtering
    
    detP <- detP[match(featureNames(MSet.QNorm), rownames(detP)),]
    
    #removing probes failed in one or more samples
    keep <- rowSums(detP<0.01) == ncol(MSet.QNorm)
    
    MSet.QNorm.filt <- MSet.QNorm[keep,]
    
    #removing probes on the sex chromosome
    annotEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    
    keep <- !(featureNames(
      MSet.QNorm.filt) %in% annotEPIC$Name[annotEPIC$chr %in% c("chrX", "chrY")])
    
    MSet.QNorm.filt <- MSet.QNorm.filt[keep,]
    
    #removing probes with SNPs at CpG site:
    
    MSet.QNorm.filt <- dropLociWithSnps(MSet.QNorm.filt)
    
    #exclude cross-reactive probes and polymorphic probes-> El ususario debe tener este documento! o debe estar subido a algún sitio para que el pipeline lo lea automáticamente desde dónde sea
    reactprobesDir <- "C:/Users/alumne/Desktop/TFM/data/shinyApp/"
    reactive.probes <- read.table(file=file.path(reactprobesDir, "reactive_probes.txt"))
    
    keep <- !(featureNames(MSet.QNorm.filt) %in% reactive.probes$V1)
    
    MSet.QNorm.filt <- MSet.QNorm.filt[keep,]
    
    return(MSet.QNorm.filt)
  }
  
  
  #Beta Vals Loading
  getBetaValsFunc <- reactive({
    pd <- pdFunc(input)
    MSet.QNorm.filt <- MSet.QNorm.filt.Func(input)
    betaVals <- getBeta(MSet.QNorm.filt)
    colnames(betaVals) <- pd$Sample_Name
    data.frame(head(betaVals))
  })
    
  #M values output  
  getMValsFunc <- reactive({
    pd <- pdFunc(input)
    MSet.QNorm.filt <- MSet.QNorm.filt.Func(input)
    MVals <- getM(MSet.QNorm.filt)
    colnames(MVals) <- pd$Sample_Name
    data.frame(head(MVals))
  })
  
  
  # DMR table output
  DMRFuncTab <- reactive ({
    targets <- input$targetsFile
    pd <- pdFunc(input)
    RGset <- RGsetFunc(input)
    
    detP <-detectionP(RGset)
    keep <- colMeans(detP)<0.05
    RGset <- RGset[,keep]
    
    GRset.funnorm <- preprocessFunnorm(RGset)
    
    designMatrix <- model.matrix(~pd$status)
    dmrs <- bumphunter(GRset.funnorm, design = designMatrix,
                       coef=1,
                       cutoff=0.5,
                       type="Beta")
    head(dmrs$table)
  })
  
  
  #DMR coef output
  DMRFuncCoef <- reactive ({
    targets <- input$targetsFile
    pd <- pdFunc(input)
    RGset <- RGsetFunc(input)
    
    detP <-detectionP(RGset)
    keep <- colMeans(detP)<0.05
    RGset <- RGset[,keep]
    
    GRset.funnorm <- preprocessFunnorm(RGset)
    
    designMatrix <- model.matrix(~pd$status)
    dmrs <- bumphunter(GRset.funnorm, design = designMatrix,
                       coef=1,
                       cutoff=0.5,
                       type="Beta")
    head(dmrs$coef)
  })
  
  
  #dmp output  
  dmpFunc <- reactive({
    MSet.QNorm.filt <- MSet.QNorm.filt.Func(input)
    betaVals <- getBeta(MSet.QNorm.filt)
    
    pd = pdFunc(input)
    
    dmp = dmpFinder(betaVals, pheno = pd$status,
                    type = "categorical")
    head(dmp)
  })
  
  
  #annotation function
  annotFunc <- reactive({
    betaVals <- getBetaValsFunc(input)
    
    annotEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    
    annotEPIC.df <- as.data.frame(annotEPIC)
    
    betaVals.annot <- merge(
      betaVals,
      annotEPIC.df,
      by.x = "X",
      by.y = "Name"
    )
    genes.list.annot <- betaVals.annot$UCSC_RefGene_Name
    
    head(genes.list.annot)
  })
  
  
  
  ######Outputs
  
  #(1) Targets table output
  output$contents <- renderDataTable({
    
    inFile <- input$targetsFile
    
    if(is.null(inFile)){
      return(NULL)
    }
    
    read.csv(inFile$datapath, header=input$header, sep=input$sep,
             quote = input$quote, stringsAsFactors = F
    )
  })
  
  
  # 2) QC plot
  output$QCplot <- renderPlot({

   RGset = RGsetFunc(input)

   MSet = MSetFunc(input)

   pd = pdFunc(input)

   qc <- getQC(MSet)

   plotQC(qc)

   #warning signal for "bad" samples
   #modifying plotQC() function from github::hansenlab/minfi
   warningBadSample <- function(qc, badSampleCutoff = 10.5) {
     meds <- (qc$mMed + qc$uMed)/2
     whichBad <- which(meds < badSampleCutoff)
     if (length(whichBad) > 0) {
       print(paste("Warning:",
                   sprintf("Warning: There are %s samples with bad quality to be considered during the following analyses",
                   length(whichBad)))
       )
     }

   }
   output$warningBadSample <- renderText({
     warningBadSample(qc)
   })

   output$downloadPlot <- downloadHandler(
     filename = function() {
       paste("Download QC plot", "png", sep=".")
     },
     content = function(file) {
       png(file)
       print(plotQC(qc))
       dev.off()
     }
   )
  })

    
  #3) Density Plot
  output$DensPlot <- renderPlot({

    MSet = MSetFunc(input)

    pd = pdFunc(input)

    densityPlot(MSet, sampGroups = pd$Sample_Group)

    output$downloadDensPlot <- downloadHandler(
     filename = function() {
      paste("Download Density plot", "png", sep=".")
      },
     content = function(file) {
       png(file)
       print(densityPlot(MSet, sampGroups = pd$Sample_Group,
             main="Beta Density Plot"))
       dev.off()
     },
     contentType='image/png'
    )
  })

  
  #4) PCA Plot
  output$PCAPlot <- renderPlot({

    MSet.QNorm.filt = MSet.QNorm.filt.Func(input)

    pcs <- prcomp(t(getM(MSet.QNorm.filt)))

    library(ggfortify)
    print(autoplot(pcs, label=TRUE))

    output$downloadPCAPlot <- downloadHandler(
      filename = function() {
        paste("Download PCA plot", "png", sep=".")
      },
      content = function(file) {
        png(file)
        print(autoplot(pcs, label=T))
        dev.off()
        # library(ggplot2)
        # device = function(..., width, height) grDevices::png(..., width=width, height = height)
        # ggsave(file, densPlot, device=device)
      },
      contentType='image/png'
    )
  })

  
  
  #5) Heatmap
  output$Heatmap <- renderPlot ({

    betaVals <- getBetaValsFunc(input)
    print(
    Heatmap(t(betaVals[1:100, 1:6]),
            name = "CpGs beta values",
            column_title_gp = gpar(fontsize=9),
            col=colorRamp2(c(0, 0.5, 1),c("blue", "white", "red")),
            show_column_names = F,
            cluster_rows = T,
            cluster_columns = T,
            column_dend_side = "top",
            #column_names_side = "bottom",
            column_names_gp = gpar(fontsize = 4),
            clustering_distance_rows = "pearson",
            clustering_distance_columns = "euclidean",
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2",
            #top_annotation = ha,
            column_title = "Heatmap of CpGs beta values"
    )
    )
    output$downloadHeatmap <- downloadHandler(
      filename = function() {
        paste("Download heatmap", "png", sep=".")
      },
      content = function(file) {
        png(file)
        print(
          heatmap
        )
        dev.off()
      },
      contentType='image/png'
    )
  })
    
  
  #6) Cluster plot
  output$Cluster <- renderPlot({
    betaVals <- getBetaValsFunc(input)

    fit <- pvclust(betaVals[1:500, 1:6], method.hclust="ward",
                   method.dist="euclidean"
           )
    plot(fit)

    output$downloadCluster <- downloadHandler(
          filename = function() {
            paste("Download hierarchical Cluster", "png", sep=".")
          },
          content = function(file) {
            png(file)
            print(
              plot(fit)
            )
            dev.off()
          },
          contentType='image/png'
        )
  })



  #7) Beta Values Table
    output$betas <- renderTable({

      getBetaValsFunc()

      output$downloadBetas <- downloadHandler(
        filename = function() {
          paste("betas", ".csv", sep="")
        },

        content = function(file) {
          write.csv(betas, file)
        }
      )

    })

  
  #8) M Values Table
  output$M <- renderTable({
  
   getMValsFunc()
  
   output$downloadM <- downloadHandler(
     filename = function() {
       paste("Mvals", ".csv", sep="")
     },
  
     content = function(file) {
       write.csv(Ms, file)
     }
   )
  })
  
  
  #9) DMR table output
  output$DMRtab <- renderTable({

    DMRFuncTab()

    output$downloadDMRtable <- downloadHandler(
      filename = function() {
        paste("DMRs_table", ".csv", sep="")
      },

      content = function(file) {
        write.csv(dmrs$table, file)
      }
    )
    })

  
  #10) DMR coef output  
  output$DMRcoef <- renderTable({

    DMRFuncTab()

    output$downloadDMRcoef <- downloadHandler(
      filename = function() {
        paste("DMRs_coef", ".csv", sep="")
      },
      content = function(file) {
        write.csv(dmrs$coef, file)
      }
    )
  })
    
  
  #11) dmp table output
  output$dmp <- renderTable({
    dmpFunc()

    output$downloaddmp <- downloadHandler(
         filename = function() {
           paste("dmp", ".csv", sep="")
         },

         content = function(file) {
           write.csv(dmp, file)
         }
       )
  })


  #12) Annotation table output
  output$annotation <- renderTable({

    annotFunc()

    output$downloadAnnotation <- downloadHandler(
             filename = function() {
               paste("genes_annotations", ".csv", sep="")
             },

             content = function(file) {
               write.csv(genes.list.annot, file)
             }
    )
  })

})
  

