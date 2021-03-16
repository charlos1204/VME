# Load packages ----
library(shiny)
library(shinydashboard)
library(ccrepe) # ccrepe is from bioconductor
library(visNetwork)
library(vegan)
library(ape)
library(cluster)
library(fpc)
library(GUniFrac)
library(geiger)
library(ggplot2)
library(dichromat)
library(reshape2)
library(ggiraph)
library(cowplot)
library(DT)
library(shinycssloaders)
library(stringr)

library(MASS)

cat("To exit press ctrl + C\n")
cat("Copy and paste the link into a browser:\n")

# Source helpers ----
#source("/home/carlos/Documents/pulications/drafts/2020/Community_Explorer/SuppMat/Docker/app/processingfunctions.R")
source("/root/app/processingfunctions.R") # change to the correct path
tooltip_css <- "background-color:gray;color:white;font-style:italic;padding:10px;border-radius:5px;"

# User interface ----

header <- dashboardHeader(title = "Community Explorer")
sidebar <- dashboardSidebar(
  tags$head(tags$style(HTML(".sidebar { height: 90vh; overflow-y: auto; }"))),
  # sidebarSearchForm(textId = "searchText", buttonId = "searchTextButton", label = "Search..."),
  sidebarMenu(
    menuItem("Example Data", tabName = "example_data", icon = icon("vials"),
             menuItem(radioButtons("microbiomedatasets", "",
              choices = c(Endesfelder="Endesfelder", Wegner="Wegner", Bazanella="Bazanella"), selected = c("None selected" = "")))),
    menuItem("Import Data", tabName = "import_data", icon = icon("arrow-circle-up"),
      menuItem(fileInput("otutable", "OTU table", multiple = FALSE,
        accept = c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))),
      menuItem(fileInput("metadatainfo", "Metada file", multiple = FALSE,
        accept = c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))),
      menuItem(fileInput("treefile", "Load tree file", multiple = FALSE,
         accept = c(".tree", ".nwk", ".txt"))),
      menuItem(fileInput("loadcomm","Community File",multiple = FALSE,
         accept = c("text/tsv", "text/tab-separated-values,text/plain", ".tsv")))
    ),

    menuItem("Export Data", tabName = "export_data", icon = icon("arrow-circle-down"),
       menuItem(downloadButton(outputId="exportadjmtx", label="Export Adjacency Matrix")),
       menuItem(downloadButton(outputId="exportdata", label="Export Data"))
    ),
    
    menuItem("Taxonomy Level", tabName = "tax_level",  startExpanded = TRUE,
       menuItem(radioButtons("level", "",choices = c(Kingdom = "Kingdom", Phylum = "Phylum",
           Class = "Class", Order = "Order", Family = "Family", Genus = "Genus",
           Species = "Species", OTU="OTU"), selected = "Genus"))
        ),
    menuItem("Settings", tabName = "setts", icon = icon("cogs"),
         menuItem(textInput("errthrshld", label = "ccrepe error treshold", value = "1e-4")),
         menuItem(textInput("pvalsccrepe", label = "ccrepe pvalue filtering", value = "0.03")),
         menuItem(textInput("scoreccrepe", label = "ccrepe score filtering", value = "0.2"))
    )
  ),
  selectizeInput('topbacteria', 'Top Bacteria', choices = c(6),
    options = list(create = TRUE)
  ),
  h5(strong("Metadata")),
  uiOutput("metadataItems"),
  hr(),
  
  actionButton("processall", "run"),
  actionButton("processmetadata", "Filter metadata")
)

body <- dashboardBody(
  fluidRow(
    #column(6, girafeOutput("pcoaplot_stackedbarplot")),
    column(6, girafeOutput("pcoaplot",height = "400px"), girafeOutput("stackedbarplot", height = "400px")),
    column(6, visNetworkOutput("network", width = "100%", height = "300px"), selectizeInput('taxonomiccomposition', 'Taxonomic composition', #
           choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), selected = "Order"),
           girafeOutput("microbiomeplot",height = "400px"))
    )
)
ui <- dashboardPage(header,sidebar,body)

# Server logic ----

server <- function(input, output, session) {
  #-----------------------------Reactive values------------------
  example_data <- reactiveValues(exmicro=NULL)
  data_pts <- reactiveValues(datapts=NULL, ldngs=NULL, top_bac=NULL, thetop=6)
  alldata <- reactiveValues(total.ab.matrix=NULL, relabumtx=NULL, otuDfRaw=NULL, taxonomy=NULL, otulist.mtx=NULL, 
                            tracking=NULL, metadata=NULL, theTree=NULL, pam_commdw=NULL, num_clustdw = NULL,
                             dw=NULL, pie_mtx=NULL, haschecks=FALSE, hassliders=FALSE, clstr_clr=NULL)
  copydata <- reactiveValues(relabumtx.copy=NULL, total.ab.matrix.copy=NULL)
  nds <- reactiveValues(nodes = NULL)
  communities <- reactiveValues(comm=NULL, hascommunities=FALSE)

  #---------------------------------------ObserveEvents---------------------------------------

  # the main observeEvent 
  observeEvent(input$processall, {
    example_data$exmicro <- "No example"
    all <- get.alldata(input$otutable$datapath, input$level, input$metadatainfo$datapath, input$treefile$datapath, example_data$exmicro)
    alldata$total.ab.matrix <- all@total.ab.matrix
    alldata$relabumtx <- all@relabumtx
    alldata$otuDfRaw <- all@otuDfRaw
    alldata$taxonomy <- all@taxonomy
    alldata$otulist.mtx <- all@otulist.mtx
    alldata$tracking <- all@tracking
    alldata$metadata <- all@metadata
    alldata$theTree <- all@theTree
    alldata$pam_commdw <- all@pam_commdw
    alldata$num_clustdw <- all@num_clustdw
    alldata$dw <- all@dw
    alldata$pie_mtx <- all@pie_mtx

    # Estimation of the network graph
    alldata$adjmtx <- get.adjMtx(alldata$relabumtx, as.numeric(input$errthrshld), as.numeric(input$pvalsccrepe), as.numeric(input$scoreccrepe), example_data$exmicro)
    nds$nodes <- data.frame(id = 1:ncol(alldata$adjmtx), label = colnames(alldata$adjmtx))

    # PCoA and stacked bar plot
    dw_pcoa <- pcoa(all@dw)
    tmp_pcoa <- get.pcoa(alldata$pam_commdw, alldata$num_clustdw, alldata$tracking, example_data$exmicro, dw_pcoa, alldata$relabumtx, dir.axis2 = -1)
    data_pts$datapts <- as.data.frame(tmp_pcoa[1])
    data_pts$ldngs <-tmp_pcoa[2]
    # to be implemented
    # xy <- metaMDS(all@dw,k=2,trymax=100)
    # data_pts$datapts[,'Axis.1'] <- xy$points[,1]
    # data_pts$datapts[,'Axis.2'] <- xy$points[,2]

    if(example_data$exmicro == "Bazanella"){
      metalbl <- unique(alldata$metadata[, 'Group'])
      for(lbl in metalbl){
        subjects <- names(alldata$metadata[which(alldata$metadata == lbl),])
        data_pts$datapts[subjects,"cluster"] <- lbl
      }

      data_pts$top_bac <- get.stackedbarplot.metadata(alldata$relabumtx, data_pts$datapts, data_pts$thetop)
      alldata$clstr_clr <- ggplotColours(length(metalbl))
    } else {
      data_pts$top_bac <- get.stackedbarplot(alldata$relabumtx, alldata$pam_commdw, alldata$num_clustdw, data_pts$thetop)
      alldata$clstr_clr <- ggplotColours(alldata$num_clustdw)
    }

    # copy objects
    copydata$relabumtx.copy <- makeCopyMtx(alldata$relabumtx)
    copydata$total.ab.matrix.copy <- makeCopyMtx(alldata$total.ab.matrix)
  })

  # ObserveEvent for filter by metadata
  observeEvent(input$processmetadata, {
    idx <- NULL
    #checks
    if(!is.null(input$checkmetadataItems)){
      idx <- get.filter.check.idx.metadata(input$checkmetadataItems, alldata$metadata)  
    }

    #sliders
    if (alldata$hassliders) {
      if(!is.null(idx)){
        filter.metadata <- alldata$metadata[idx,]
        sldrs <- elementsitems()$"sliders"
        idx <- get.filter.slider.metadata(sldrs, filter.metadata, input, alldata$tracking)
      } else {
        sldrs <- elementsitems()$"sliders"
        idx <- get.filter.slider.metadata(sldrs, alldata$metadata, input, alldata$tracking)
      }
    }

    if(length(idx) > nrow(copydata$total.ab.matrix.copy)){
      idx <- intersect(idx, rownames(copydata$total.ab.matrix.copy))
    }

    filter.total.ab.matrix <- copydata$total.ab.matrix.copy[intersect(idx,rownames(copydata$total.ab.matrix.copy)), , drop=FALSE]
    if(input$level != "OTU"){
      alldata$relabumtx <- get.relabumtx.metadata(filter.total.ab.matrix, copydata$total.ab.matrix.copy)
    }else{
      if(example_data$exmicro == "Wegner"){
        alldata$relabumtx <- wisconsin(sqrt(filter.total.ab.matrix))
      } else {
        alldata$relabumtx <- get.relabumtx.metadata.OTU(filter.total.ab.matrix)
      }
    }

    if(alldata$theTree != "no tree"){
      if(ncol(alldata$relabumtx) != nrow(alldata$otulist.mtx)){
        otulist.mtx <- update.otulist(alldata$otulist.mtx, colnames(alldata$relabumtx))
      }else{
        otulist.mtx <- alldata$otulist.mtx
      }
      the_tree <- read.tree(file=alldata$theTree)
      mtx_raw <- alldata$relabumtx
      colnames(mtx_raw) <- rownames(otulist.mtx)
      # samples as columns
      mtx_raw <- t(mtx_raw)
      mtx_the_tree <- NULL
      tryCatch({
        toremove <- name.check(the_tree, mtx_raw)$tree_not_data
        mtx_the_tree <- drop.tip(the_tree, toremove)
      }, error=function(e){})
      
      if(is.null(mtx_the_tree)){
        mtx_the_tree <- the_tree
      }

      if(example_data$exmicro == "Endesfelder"){
        # return samples as rows
        mtx_raw <- t(mtx_raw)
        mtx_unifracDW <- unifrac2(mtx_raw, mtx_the_tree)
        mtx_unifracDW[is.na(mtx_unifracDW)] = 0
        alldata$dw <- mtx_unifracDW
      } else {
        # return samples as rows
        mtx_raw <- t(mtx_raw)
        mtx_unifracDW <- GUniFrac(mtx_raw, mtx_the_tree, alpha=c(0.5,1))$unifracs
        mtx_unifracDW[is.na(mtx_unifracDW)] = 0
        alldata$dw <- mtx_unifracDW[,,"d_1"]
      }
    } else {
      mtx_raw <- alldata$relabumtx
      alldata$dw <- vegdist(mtx_raw, method = "horn")
      alldata$dw[is.na(alldata$dw)] <- 0
      alldata$dw <- as.matrix(alldata$dw)
    }

    # calculating the number of cluster and clustering with pam
    if(example_data$exmicro == "Wegner"){
      alldata$num_clustdw <-pamk(alldata$dw, krange = 1:6, criterion = "ch")$nc
      alldata$pam_commdw <- pam(alldata$dw, k = alldata$num_clustdw , metric = "euclidean")$cluster
    } else {
      alldata$num_clustdw <-pamk(alldata$dw, krange = 1:6, criterion = "ch")$nc
      alldata$pam_commdw <- pam(alldata$dw, k = alldata$num_clustdw)$cluster
    }

    # PCoA and stacked bar plot
    dw_pcoa <- pcoa(alldata$dw)
    tmp_pcoa <- get.pcoa(alldata$pam_commdw, alldata$num_clustdw, alldata$tracking, example_data$exmicro, dw_pcoa, alldata$relabumtx, dir.axis2 = -1)
    data_pts$datapts <- as.data.frame(tmp_pcoa[1])
    data_pts$ldngs <-tmp_pcoa[2]
    #to be implemented
    # xy <- metaMDS(alldata$dw,k=2,trymax=100)
    # data_pts$datapts[,'Axis.1'] <- xy$points[,1]
    # data_pts$datapts[,'Axis.2'] <- xy$points[,2]

    data_pts$top_bac <- get.stackedbarplot(alldata$relabumtx, alldata$pam_commdw, alldata$num_clustdw, data_pts$thetop)
    alldata$clstr_clr <- ggplotColours(alldata$num_clustdw)
    
    if(!communities$hascommunities){
      alldata$adjmtx <- get.adjMtx(alldata$relabumtx, as.numeric(input$errthrshld), as.numeric(input$pvalsccrepe), as.numeric(input$scoreccrepe), example_data$exmicro)
      
      if(sum(rowSums(alldata$adjmtx)) > 2){
        nds$nodes <- data.frame(id = 1:ncol(alldata$adjmtx), label = colnames(alldata$adjmtx))
        edges <- as.data.frame(makeEdges(alldata$adjmtx))
        names(edges)[1] <- "from"
        names(edges)[2] <- "to"
        visNetworkProxy("network") %>%
          visUpdateNodes(nds$nodes, updateOptions = FALSE) %>%
          visUpdateEdges(edges)
      }
    }
  
    bacteria_filtered <- as.matrix(colSums(alldata$relabumtx))
    bacteria_filtered <- bacteria_filtered[bacteria_filtered[,1] > 0.01,,drop=FALSE]
    alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, rownames(bacteria_filtered), input$level, input$taxonomiccomposition, "Taxonomic composition")
  })

  ##########################
  #     Demo datasets      #
  ##########################
  observeEvent(input$microbiomedatasets, {
    example_data$exmicro <- "No example"
    if(input$microbiomedatasets == "Endesfelder"){
      example_data$exmicro <- "Endesfelder"
      exdatapth <- "/root/datasets/Endesfelder/otu_data.tsv"
      exmetadatapth <- "/root/datasets/Endesfelder/metadata.txt"
      extreepth <- "/root/datasets/Endesfelder/97_otus.tree"
      #exdatapth <- "datasets/Endesfelder/otu_data.tsv"
      #exmetadatapth <- "datasets/Endesfelder/metadata.txt"
      #extreepth <- "datasets/Endesfelder/97_otus.tree"
      all <- get.alldata(exdatapth, input$level, exmetadatapth, extreepth, example_data$exmicro)
    } else  if(input$microbiomedatasets == "Wegner"){
      example_data$exmicro <- "Wegner"
      exdatapth <- "/root/datasets/Wegner/data_set.tsv"
      exmetadatapth <- "/root/datasets/Wegner/metadata.txt"
      #exdatapth <- "datasets/Wegner/data_set.tsv"
      #exmetadatapth <- "datasets/Wegner/metadata.txt"
      extreepth <- input$treefile$datapath
      all <- get.alldata(exdatapth, input$level, exmetadatapth, extreepth, example_data$exmicro)
    } else  if(input$microbiomedatasets == "Bazanella"){
      example_data$exmicro <- "Bazanella"
      exdatapth <- "/root/datasets/Bazanella/bazanella_mth12.tsv"
      exmetadatapth <- "/root/datasets/Bazanella/Metafile_Month12.txt"
      extreepth <- "/root/datasets/Bazanella/Rooted_tree_Month12.nwk"
      #exdatapth <- "datasets/Bazanella/bazanella_mth12.tsv"
      #exmetadatapth <- "datasets/Bazanella/Metafile_Month12.txt"
      #extreepth <- "datasets/Bazanella/Rooted_tree_Month12.nwk"
      all <- get.alldata(exdatapth, input$level, exmetadatapth, extreepth, example_data$exmicro)
    }
    alldata$total.ab.matrix <- all@total.ab.matrix
    alldata$relabumtx <- all@relabumtx
    alldata$otuDfRaw <- all@otuDfRaw
    alldata$taxonomy <- all@taxonomy
    alldata$otulist.mtx <- all@otulist.mtx
    alldata$tracking <- all@tracking
    alldata$metadata <- all@metadata
    alldata$theTree <- all@theTree
    alldata$pam_commdw <- all@pam_commdw
    alldata$num_clustdw <- all@num_clustdw
    alldata$dw <- all@dw
    alldata$pie_mtx <- all@pie_mtx
    
    # Estimation of the network graph
    alldata$adjmtx <- get.adjMtx(alldata$relabumtx, as.numeric(input$errthrshld), as.numeric(input$pvalsccrepe), as.numeric(input$scoreccrepe), example_data$exmicro)
    nds$nodes <- data.frame(id = 1:ncol(alldata$adjmtx), label = colnames(alldata$adjmtx))
    
    # PCoA and stacked bar plot
    dw_pcoa <- pcoa(all@dw)
    tmp_pcoa <- get.pcoa(alldata$pam_commdw, alldata$num_clustdw, alldata$tracking, example_data$exmicro, dw_pcoa, alldata$relabumtx, dir.axis2 = -1)
    data_pts$datapts <- as.data.frame(tmp_pcoa[1])
    data_pts$ldngs <-tmp_pcoa[2]
    # to be implemented
    # xy <- metaMDS(all@dw,k=2,trymax=100)
    # data_pts$datapts[,'Axis.1'] <- xy$points[,1]
    # data_pts$datapts[,'Axis.2'] <- xy$points[,2]
    
    if(example_data$exmicro == "Bazanella"){
      metalbl <- unique(alldata$metadata[, 'Group'])
      for(lbl in metalbl){
        subjects <- names(alldata$metadata[which(alldata$metadata == lbl),])
        data_pts$datapts[subjects,"cluster"] <- lbl
      }
      
      data_pts$top_bac <- get.stackedbarplot.metadata(alldata$relabumtx, data_pts$datapts, data_pts$thetop)
      alldata$clstr_clr <- ggplotColours(length(metalbl))
    } else {
      data_pts$top_bac <- get.stackedbarplot(alldata$relabumtx, alldata$pam_commdw, alldata$num_clustdw, data_pts$thetop)
      alldata$clstr_clr <- ggplotColours(alldata$num_clustdw)
    }
    
    # copy objects
    copydata$relabumtx.copy <- makeCopyMtx(alldata$relabumtx)
    copydata$total.ab.matrix.copy <- makeCopyMtx(alldata$total.ab.matrix)
  })
  
  #obsereEvent for network community
  observeEvent(input$loadcomm, {
    communities$comm <- load.communities(input$loadcomm$datapath, alldata$adjmtx)
    communities$hascommunities <- TRUE
    comm_color <- communities$comm[,"color"]
    gpo <- communities$comm[,"gpo"]
    nds$nodes["Community"] <- gpo
    nds$nodes["color"] <- comm_color
    visNetworkProxy("network") %>%
      visUpdateNodes(nds$nodes, updateOptions = FALSE)
  
    tmp_comm <- sort(unique(gpo))
    mtx_pie <- c()
    for(cm in tmp_comm){
      comm_idx <- which(communities$comm[, "gpo"] == cm)
      filter_com <- rownames(communities$comm)[comm_idx]
      tx <- get.pie.matrix(alldata$taxonomy, filter_com, input$level, input$taxonomiccomposition, cm)
      mtx_pie <- rbind(mtx_pie, tx)
    }

    alldata$pie_mtx <- mtx_pie
  })

  # ObserveEvent of metadata checks to set sliders
  observeEvent(input$checkmetadataItems, {
    mtd <- alldata$metadata
    if(length(input$checkmetadataItems) > 1) {
      v <- rownames(which(mtd==input$checkmetadataItems[1], arr.ind = TRUE))
      for (el in input$checkmetadataItems) {
        v <- intersect(rownames(which(mtd == el, arr.ind = TRUE)), v)
      }

      if (length(v) == 0) {
        for (el in input$checkmetadataItems) {
          v <- c(rownames(which(mtd == el, arr.ind = TRUE)), v)
        }

        #sliders update
        if (alldata$hassliders) {
          sldrs <- elementsitems()$"sliders"
          nms <- names(sldrs)
          mtd2 <- mtd[v,]
          for (nn in nms) {
            nmid <- paste("slider", nn, sep = "_")
            vmax <- max(as.numeric(mtd2[!is.na(mtd2[, nn]), nn]))
            vmin <- min(as.numeric(mtd2[!is.na(mtd2[, nn]), nn]))
            updateSliderInput(session = session, nmid, max = vmax, min = vmin, value = c(vmin, vmax))
          }
        }
      } else {
        #sliders update
        if(alldata$hassliders) {
          sldrs <- elementsitems()$"sliders"
          nms <- names(sldrs)
          mtd2 <- mtd[v,]
          for(nn in nms){
            nmid <- paste("slider", nn, sep="_")
            vmax <- max(as.numeric(mtd2[!is.na(mtd2[, nn]), nn]))
            vmin <- min(as.numeric(mtd2[!is.na(mtd2[, nn]), nn]))
            updateSliderInput(session = session, nmid, max=vmax, min = vmin, value = c(vmin, vmax))
          }
        }
      }
    } else {
      #sliders update
      if(alldata$hassliders){
        sldrs <- elementsitems()$"sliders"
        nms <- names(sldrs)
        idxs <- rownames(which(mtd==input$checkmetadataItems, arr.ind = TRUE))
        mtd2 <- mtd[idxs,]
        for(nn in nms){
          nmid <- paste("slider", nn, sep="_")
          vmax <- max(as.numeric(mtd2[!is.na(mtd2[, nn]), nn]))
          vmin <- min(as.numeric(mtd2[!is.na(mtd2[, nn]), nn]))
          updateSliderInput(session = session, nmid, max=vmax, min = vmin, value = c(vmin, vmax))
        }
      }
    }
  })

  # ObserveEvent network community
  observeEvent(input$network_selectedBy, {
    if (input$network_selectedBy != "") {
      the_tree <- read.tree(file = alldata$theTree) #
      comm_idx <- which(communities$comm[, "gpo"] == input$network_selectedBy)
      filter_com <- rownames(communities$comm)[comm_idx]
      filter_relabumtx <- alldata$relabumtx[, filter_com, drop=FALSE]
      
      cls <- colnames(filter_relabumtx)

      if(ncol(filter_relabumtx) != nrow(alldata$otulist.mtx)){
        otulist.mtx <- update.otulist(alldata$otulist.mtx, colnames(filter_relabumtx))
      }

      colnames(filter_relabumtx) <- rownames(otulist.mtx)
      # samples as columns
      filter_relabumtx <- t(filter_relabumtx)
      mtx_the_tree <- NULL
      tryCatch({
        toremove <- name.check(the_tree, filter_relabumtx)$tree_not_data
        mtx_the_tree <- drop.tip(the_tree, toremove)
      }, error=function(e){})

      if(is.null(mtx_the_tree)){
        mtx_the_tree <- the_tree
      }

      # return samples as rows
      dw <- NULL
      if(example_data$exmicro == "Endesfelder"){
        filter_relabumtx <- t(filter_relabumtx)
        mtx_unifracDW <- unifrac2(filter_relabumtx, mtx_the_tree)
        mtx_unifracDW[is.na(mtx_unifracDW)] = 0
        dw <- mtx_unifracDW
      } else {
        filter_relabumtx <- t(filter_relabumtx)
        tryCatch({
          mtx_unifracDW <- suppressWarnings(GUniFrac(filter_relabumtx, mtx_the_tree, alpha=c(0, 0.5,1))$unifracs)
          mtx_unifracDW[is.na(mtx_unifracDW)] = 0
          dw <- mtx_unifracDW[,,"d_1"]
        }, error=function(e){})
      }

      # calculating the number of cluster and clustering with pam
      if(example_data$exmicro == "Wegner"){
        alldata$num_clustdw <- pamk(dw, krange = 1:6, criterion = "ch")$nc
        alldata$pam_commdw <- pam(dw, k = alldata$num_clustdw, metric = "euclidean")$cluster
      } else {
        alldata$num_clustdw <- pamk(dw, krange = 1:6, criterion = "ch")$nc
        alldata$pam_commdw <- pam(dw, k = alldata$num_clustdw)$cluster
      }

      # PCoA and stacked bar plot
      colnames(filter_relabumtx) <- cls
      dw_pcoa <- pcoa(dw)
      tmp_pcoa <- get.pcoa(alldata$pam_commdw, alldata$num_clustdw, alldata$tracking, example_data$exmicro, dw_pcoa, filter_relabumtx, dir.axis2 = -1) #-1
      data_pts$datapts <- as.data.frame(tmp_pcoa[[1]])
      data_pts$ldngs <-tmp_pcoa[[2]]

      data_pts$top_bac <- get.stackedbarplot(filter_relabumtx, alldata$pam_commdw, alldata$num_clustdw, data_pts$thetop)
      alldata$clstr_clr <- ggplotColours(alldata$num_clustdw)

      alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, colnames(filter_relabumtx), input$level, input$taxonomiccomposition, input$network_selectedBy)
    }
  })

  # ObserveEvent pcoa to network and Stackedbar plot
  observeEvent(selected_arrow(), {
    if(startsWith(selected_arrow(),"bacteria_") & !startsWith(selected_arrow(),"track_")) {
      bac = gsub("bacteria_","", selected_arrow())
      session$sendCustomMessage(type = 'stackedbarplot_set', message = bac)
      idx = which(colnames(alldata$adjmtx) == bac)
      visNetworkProxy("network") %>%
        visSelectNodes(id=idx) %>%
        visFocus(id = idx, scale = 4)
    } else if(selected_arrow() == "unselected_arrow" & !startsWith(selected_arrow(),"track_")) {
      if(!is.null(alldata$adjmtx)){
        session$sendCustomMessage(type = 'stackedbarplot_set', message = character(0))
      }
    }
  })

  # ObserveEvent Stackedbar to network and pcoa plot
  observeEvent(selected_bar(), {
    if(startsWith(selected_bar(),"bacteria_")) {
      bac = gsub("bacteria_","",selected_bar())
      level_incomposition <- get.taxa.in.taxonomy(alldata$taxonomy, bac, input$level, input$taxonomiccomposition)
      session$sendCustomMessage(type = 'pcoaplot_set', message = bac)
      session$sendCustomMessage(type = 'microbiomeplot_set', message = level_incomposition)
      idx = which(colnames(alldata$adjmtx) == bac)
      visNetworkProxy("network") %>%
        visSelectNodes(id=idx) %>%
        visFocus(id = idx, scale = 4)
    } else if(selected_bar() == "unselected_bar" & !startsWith(selected_bar(),"track_")) {
      if(!is.null(alldata$adjmtx)){
        session$sendCustomMessage(type = 'pcoaplot_set', message = character(0))
        session$sendCustomMessage(type = 'microbiomeplot_set', message = character(0))
      }
    }
  })

  # ObserveEvent net to pcoa and Stackedbar plot
  observeEvent(input$current_node_id, {
    bacteria <- colnames(alldata$adjmtx)
    if(length(input$current_node_id$nodes) > 0){
      bac <- bacteria[input$current_node_id$nodes[[1]]]
      session$sendCustomMessage(type = 'pcoaplot_set', message = bac)
      session$sendCustomMessage(type = 'stackedbarplot_set', message = bac)
      level_incomposition <- get.taxa.in.taxonomy(alldata$taxonomy, bac, input$level, input$taxonomiccomposition)
      session$sendCustomMessage(type = 'pcoaplot_set', message = bac)
    } else {
      session$sendCustomMessage(type = 'pcoaplot_set', message = character(0))
      session$sendCustomMessage(type = 'stackedbarplot_set', message = character(0))
      session$sendCustomMessage(type = 'microbiomeplot_set', message = character(0))
    }
  })

  # ObserveEvent for tracking, this clicks on elements for tracking
  observeEvent(selected_track(), {
    if(startsWith(selected_track(),"track_")) {
      no_selected <- rownames(data_pts$datapts[-which(data_pts$datapts$track_ID == selected_track()), , drop=FALSE])
      data_pts$datapts["alpha_trans"]  <- rep(1.0, nrow(data_pts$datapts))
      data_pts$datapts[no_selected, "alpha_trans"] <- 0.1
      #test for taxonomy composition
      subjects <- rownames(data_pts$datapts[which(data_pts$datapts$track_ID == selected_track()), , drop=FALSE])
      submtx <- alldata$relabumtx[subjects, , drop=FALSE]
      submtx <- submtx[, colSums(submtx) > 0, drop=FALSE]
      alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, colnames(submtx), input$level, input$taxonomiccomposition, selected_track())
    }else if(selected_track() == "unselected_tracking"){
      if(!is.null(data_pts$datapts)){
        data_pts$datapts["alpha_trans"]  <- rep(1.0, nrow(data_pts$datapts))
        bacteria_filtered <- as.matrix(colSums(alldata$relabumtx))
        bacteria_filtered <- bacteria_filtered[bacteria_filtered[,1] > 0.01,,drop=FALSE]
        alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, rownames(bacteria_filtered), input$level, input$taxonomiccomposition, "Taxonomic composition")
      }
    }
  })

  observeEvent(selected_cluster(),{
    if(!is.null(alldata$pie_mtx)) {
      if(selected_cluster() == "unselected_cluster") {
        alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, colnames(alldata$relabumtx), input$level, input$taxonomiccomposition, "Taxonomic composition")
      } else {
          subjects <- rownames(data_pts$datapts[which(data_pts$datapts[,"cluster"] == selected_cluster()), , drop=FALSE])
          submtx <- alldata$relabumtx[subjects, , drop=FALSE]
          submtx <- submtx[, colSums(submtx) > 0.01, drop=FALSE]
          alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, colnames(submtx), input$level, input$taxonomiccomposition, selected_cluster())
      }
    }
  })

  #observeEvent for the top bacteria
  observeEvent(input$topbacteria, {
    if(!is.null(data_pts$datapts)){
      if(input$topbacteria!=""){
        data_pts$thetop <- as.integer(input$topbacteria)
        data_pts$top_bac <- get.stackedbarplot(alldata$relabumtx, alldata$pam_commdw, alldata$num_clustdw, data_pts$thetop)
      }
    }
  })

  # observeEvent for the Taxonomic composition
  observeEvent(input$taxonomiccomposition, {
    if(!is.null(alldata$pie_mtx)) {
      if(communities$hascommunities){
        if (input$network_selectedBy != "") {
          # selected community
          comm_idx <- which(communities$comm[, "gpo"] == input$network_selectedBy)
          filter_com <- rownames(communities$comm)[comm_idx]
          alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, filter_com, input$level, input$taxonomiccomposition, input$network_selectedBy)
        } else if (input$network_selectedBy == "") {
          # all communities
          tmp_comm <- sort(unique(communities$comm[,"gpo"]))
          mtx_pie <- c()
          for(cm in tmp_comm){
            comm_idx <- which(communities$comm[, "gpo"] == cm)
            filter_com <- rownames(communities$comm)[comm_idx]
            tx <- get.pie.matrix(alldata$taxonomy, filter_com, input$level, input$taxonomiccomposition, cm)
            mtx_pie <- rbind(mtx_pie, tx)
          }
          alldata$pie_mtx <- mtx_pie
        }
      } else {
        # cluster is selected
        if(selected_cluster() == "unselected_cluster"){
            # on starting
            bacteria_filtered <- as.matrix(colSums(alldata$relabumtx))
            bacteria_filtered <- bacteria_filtered[bacteria_filtered[,1] > 0.01,,drop=FALSE]
            alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, rownames(bacteria_filtered), input$level, input$taxonomiccomposition, "Taxonomic composition")
        }else{
            subjects <- rownames(data_pts$datapts[which(data_pts$datapts[,"cluster"] == selected_cluster()), , drop=FALSE])
            submtx <- alldata$relabumtx[subjects, , drop=FALSE]
            submtx <- submtx[, colSums(submtx) > 0.01, drop=FALSE]
            alldata$pie_mtx <- get.pie.matrix(alldata$taxonomy, colnames(submtx), input$level, input$taxonomiccomposition, selected_cluster())
        }
      }
    }
  })

  observeEvent(input$level, {
    if(get.taxonomy.position(input$level) < get.taxonomy.position("Order")){
      updateSelectInput(session, "taxonomiccomposition", selected = input$level)
    }
  })
  
  #---for unselect the example data
  observeEvent(input$otutable, {
      updateRadioButtons(session, "microbiomedatasets", 
        choices = c(Endesfelder="Endesfelder",Wegner="Wegner",Bazanella="Bazanella"), 
        selected = c("None selected" = ""))
  })
  
  observeEvent(input$metadatainfo, {
    updateRadioButtons(session, "microbiomedatasets", 
                       choices = c(Endesfelder="Endesfelder",Wegner="Wegner",Bazanella="Bazanella"), 
                       selected = c("None selected" = ""))
  })
  
  observeEvent(input$treefile, {
    updateRadioButtons(session, "microbiomedatasets", 
                       choices = c(Endesfelder="Endesfelder",Wegner="Wegner",Bazanella="Bazanella"), 
                       selected = c("None selected" = ""))
  })
  
  # # observeEvent for the text search to be define
  # observeEvent(input$searchTextButton ,{
  #   print(input$searchText) 
  # })
  
  #---------------------------------------Reactives---------------------------------------
  # selection of lines
  selected_arrow <- reactive({
    if(length(input$pcoaplot_selected) > 0){
      if(!startsWith(input$pcoaplot_selected, "track_")){
        return(paste("bacteria", input$pcoaplot_selected, sep="_"))
      }
    } else {
      return("unselected_arrow")
    }
  })
  
  selected_bar <- reactive({
    if(length(input$stackedbarplot_selected) > 0){
      return(paste("bacteria",input$stackedbarplot_selected, sep="_"))
    }else{
      return("unselected_bar")
    }
  })
  
  # selection for tracking
  selected_track <- reactive({
    if(length(input$pcoaplot_selected) > 0){
      return(paste(input$pcoaplot_selected))
    } else {
      return("unselected_tracking")
    }
  })
  
  # metadata items
  elementsitems <- reactive({
    if(!is.null(alldata$metadata)){
      mtd <- alldata$metadata
      get.metadata.items(mtd)
    }
  })

  # selection of cluster
  selected_cluster <- reactive({
    if(length(input$pcoaplot_key_selected) > 0){
      return(input$pcoaplot_key_selected)
    } else {
      return("unselected_cluster")
    }
  })
  #---------------------------------------Renders---------------------------------------
  #metadata items
  output$metadataItems <- renderUI({

    chcs_chk <- elementsitems()$"checkbox"
    sldrs <- elementsitems()$"sliders"

    if(length(chcs_chk) > 0){
      alldata$haschecks = TRUE
    }

    if(length(sldrs) > 0){
      alldata$hassliders=TRUE
    }

    elements <- list()
    if(alldata$haschecks){
      elements[[1]] <- checkboxGroupInput("checkmetadataItems", "", choices = chcs_chk)
    }

    if(alldata$hassliders){
      nms <- names(sldrs)
      cnt <- 2
      for(nn in nms){
        tmp <- sldrs[[nn]]
        sld_max <- tmp[1]
        sld_min <- tmp[2]
        nmid <- paste("slider", nn, sep="_")
        elements[[cnt]] <- sliderInput(nmid, nn, min=sld_min, max=sld_max, value=c(sld_min,sld_max))
        cnt <- cnt + 1
      }
    }

    tagList(elements)
  })
  
  # network interaction plot
  output$network <- renderVisNetwork({
    if(!is.null(alldata$adjmtx)){
      edges <- as.data.frame(makeEdges(alldata$adjmtx))
      names(edges)[1] <- "from"
      names(edges)[2] <- "to"

      visNetwork(nds$nodes, edges, height = "300px") %>% #
       visOptions(highlightNearest = TRUE) %>%
       visOptions(selectedBy = list(variable = "Community")) %>%
       visEvents(click = "function(nodes){
          Shiny.onInputChange('current_node_id', nodes);;}") %>%
       visPhysics(stabilization = FALSE)%>%
       visNodes(size = 10) %>%
       visEdges(smooth = FALSE)
    }
  })

  # pcoa plots
  output$pcoaplot<- renderGirafe({
    if(!is.null(data_pts$datapts)){
      datapts <- as.data.frame(data_pts$datapts)
      data_pts$ldngs <- as.data.frame(data_pts$ldngs)

      vals <- setNames(alldata$clstr_clr, unique(datapts[,'cluster']))
      legendId <- setNames(names(vals), unique(datapts[,'cluster']))
      legendTooltip <- setNames(names(vals), unique(datapts[,'cluster']))
      pltpcoa <- ggplot(data = datapts, aes(x=Axis.1,y=Axis.2)) +
        isolate(geom_point_interactive(data = datapts, aes(x=Axis.1,y=Axis.2, tooltip = track_label, 
        data_id= track_ID, colour = datapts[,'cluster']), alpha=datapts$alpha_trans, size = 3)) +
        labs(color = "Cluster") + xlab("PC1") + ylab("PC2") +
        scale_color_manual_interactive(
          data_id = legendId,
          tooltip = legendTooltip,
          values=vals) +
        isolate(geom_segment_interactive(data= data_pts$ldngs, aes(x = 0,y = 0, xend = data_pts$ldngs[, 'Axis.1'],
        yend = data_pts$ldngs[, 'Axis.2'], tooltip = data_pts$ldngs[, 'variable'], data_id=data_pts$ldngs[, 'variable']),
        arrow = arrow(length = unit(1 / 2, "picas")), show.legend = F, size=1.3)) +
        theme_minimal(base_size = 20)

      girafe(ggobj = pltpcoa,
             width_svg = 12, height_svg = 6.8,
             options = list(opts_toolbar(saveaspng = FALSE),
                            opts_tooltip(css = tooltip_css),
                            opts_hover(css = "fill:orange;stroke:orange;cursor:pointer"),
                            opts_selection(type = "single", css = "fill:red;stroke:red;cursor:pointer"),
                            opts_selection_key(css = "stroke:black;r:5pt;"),
                            opts_hover_key(css = "stroke:black;r:5pt;cursor:pointer;")
             ),
      )
    }
  })
  
  # bar plots
  output$stackedbarplot <- renderGirafe({
    if(!is.null(data_pts$datapts)){
      top6.colors <- colorRampPalette(c("darkblue",  "lightgray"))
      clr <- top6.colors(data_pts$thetop)
      pltstkbar <- ggplot(data=data_pts$top_bac, aes(x=cluster, y=value, fill=bacteria)) + 
        geom_bar(stat="identity", position = "stack") +
        scale_fill_manual(values=clr) + labs(fill = "Bacteria", y="Relative Abundance") +
        geom_col_interactive(aes(tooltip = bacteria, data_id = bacteria))  + 
        theme_minimal(base_size = 20)

      girafe(ggobj = pltstkbar,
             width_svg = 12, height_svg = 6.8,
             options = list(opts_toolbar(saveaspng = FALSE),
                            opts_tooltip(css = tooltip_css),
                            opts_hover(css = "fill:orange;stroke:orange;cursor:pointer"),
                            opts_selection(type = "single", css = "fill:red;stroke:red;cursor:pointer")
             ),
      )

    }
  })
  
  # pie char of bacteria
  output$microbiomeplot <- renderGirafe({
    if(!is.null(alldata$pie_mtx)) {
      coms <- sort(unique(alldata$pie_mtx[,"community"]))
      if(length(coms) > 1) {
        microbpie <- list()
        for(com in coms){
          microb <- as.data.frame(alldata$pie_mtx[which(alldata$pie_mtx[,"community"]==com), ,drop=FALSE])
          plt <- ggplot(data=microb, aes(x="", y=value, fill=variable)) +
            geom_bar(stat="identity", position = "stack") +
            geom_col_interactive(aes(tooltip = variable, data_id = variable)) +
            coord_polar("y", start = 0) + labs(fill = "") +  xlab("") + ylab("") +
            ggtitle(com) +
            theme(axis.text.x=element_blank(), legend.position="none")
          microbpie[[com]] <- plt
        }
        girafe(ggobj=plot_grid(plotlist=microbpie, ncol = 2), 
               width_svg = 12, height_svg = 6.8,
               options = list(opts_toolbar(saveaspng = FALSE),
                              opts_tooltip(css = tooltip_css),
                              opts_hover(css = "fill:red;stroke:red;cursor:pointer")
               ))
      } else {
        microb <- as.data.frame(alldata$pie_mtx)
        com <- unique(alldata$pie_mtx[,"community"])
        microbpie <- ggplot(data=microb, aes(x="", y=value, fill=variable)) +
          geom_bar(stat="identity", position = "stack") +
          geom_col_interactive(aes(tooltip = variable, data_id = variable)) +
          coord_polar("y", start = 0) + labs(fill = "") +  xlab("") + ylab("") +
          ggtitle(com) +
          theme(axis.text.x=element_blank(), legend.position="none")
        
        girafe(ggobj=microbpie, 
               width_svg = 12, height_svg = 6.8,
               options = list(opts_toolbar(saveaspng = FALSE),
                              opts_tooltip(css = tooltip_css),
                              opts_hover(css = "fill:red;stroke:red;cursor:pointer")
               ))
      }
    }
  })

  #------------------Export Section-------------------
    
  output$exportadjmtx <- downloadHandler(
      filename = function() {
       paste0("Adjacency_matrix", ".tsv",sep="")
      },
      content = function(file) {
        write.table(file=file, alldata$adjmtx, row.names=FALSE, sep = "\t", quote=FALSE)
      }
  )
    
  output$exportdata <- downloadHandler(
      filename = function() {
        paste0("data_table", ".tsv",sep="")
      },
      content = function(file) {
        write.table(file=file, alldata$tat_mtx, row.names=FALSE, sep = "\t", quote=FALSE)
      }
  )
        
}


# Run the app ----
shinyApp(ui, server)

