#                               antes tat  
setClass("tabledata", slots = c(total.ab.matrix="matrix", relabumtx="matrix", otuDfRaw="matrix", taxonomy="matrix",
                                otulist.mtx="matrix", tracking="matrix", metadata="matrix", theTree="character",
                                pam_commdw="numeric", num_clustdw="numeric", dw="matrix", pie_mtx="matrix"))
community_colors <- c("#000000","#FF6105", "#0082C8", "#FFE105", "#3CB4AF", "#F58230", "#911EB4", "#46F0F0", "#F032E6", "#D2F53C", "#E6194B", "#008080", "#E6BEFF", "#AA6E28", "#FFFAC8", "#800000", "#AAFFC3", "#808000", "#FFD7B4")
#                        -1         1          2         3           4         5          6          7

get.alldata <- function(otutable, level, metadata_table, theTree, example_data){
  # ---Processing of the OTU table and taxonomy
  otuDfRaw <- as.matrix(read.table(otutable, sep="\t", header = TRUE, check.names=FALSE))
  rownames(otuDfRaw) <- trim(otuDfRaw[,1])
  otuDfRaw <- otuDfRaw[,-1]
  
  taxonomy <- as.matrix(otuDfRaw[,ncol(otuDfRaw)])
  otuDfRaw <- otuDfRaw[,-ncol(otuDfRaw)]
  
  class(otuDfRaw) <- "numeric"
  
  taxonomy <- gsub("k__","",taxonomy)    # 1 kg
  taxonomy <- gsub("; p__",";",taxonomy) # 2 ph
  taxonomy <- gsub(";p__",";",taxonomy)
  taxonomy <- gsub("; c__",";",taxonomy) # 3 cls
  taxonomy <- gsub(";c__",";",taxonomy)
  taxonomy <- gsub("; o__",";",taxonomy) # 4 or
  taxonomy <- gsub(";o__",";",taxonomy)
  taxonomy <- gsub("; f__",";",taxonomy) # 5 fam
  taxonomy <- gsub(";f__",";",taxonomy)
  taxonomy <- gsub("; g__",";",taxonomy) # 6 gen
  taxonomy <- gsub(";g__",";",taxonomy)
  taxonomy <- gsub("; s__",";",taxonomy) # 7 sp
  taxonomy <- gsub(";s__",";",taxonomy)
  taxonomy <- gsub("\\[","",taxonomy)
  taxonomy <- gsub("\\]","",taxonomy)
  
  taxonomy <- get.tax.to.levels(taxonomy)
  otulist.mtx <- get.taxonomy.list(taxonomy, level)
  
  total.ab.matrix <- NULL
  relabumtx <- NULL
  if(level == "OTU"){
    #total.ab.matrix <- t(otuDfRaw)
    relabumtx <- get.relabumtx.OTU(t(otuDfRaw))
  } else {
    total.ab.matrix <- get.total.abundance.mtx(otuDfRaw, otulist.mtx, taxonomy)
    if(example_data == "Wegner"){
      relabumtx <- wisconsin(sqrt(total.ab.matrix))
    } else {
      relabumtx <- get.relabumtx(otuDfRaw, otulist.mtx, taxonomy)
    }
  }

  if(ncol(relabumtx) != nrow(otulist.mtx)){
    otulist.mtx <- update.otulist(otulist.mtx, colnames(relabumtx))
  }
  
  # ---End of processing of the OTU table and taxonomy
  
  # ---Calculating the pie matrix
  level2 <- "Order"
  if(get.taxonomy.position(level) < get.taxonomy.position("Order")){
    level2 <- level
  }
  
  pie_mtx <- get.pie.matrix(taxonomy, colnames(relabumtx), level, level2, "Taxonomic composition")
  # ---End calculating the pie matrix
  
  # ---processing of the metadata and tracking
  metadata <- as.matrix(read.table(metadata_table, sep="\t", header = TRUE, check.names=FALSE))
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]

  metadata <- metadata[intersect(rownames(relabumtx), rownames(metadata)),,drop=FALSE]

  tracking <- NULL
  if("track_ID" %in% colnames(metadata)){
    tracking <- metadata[,1:2]
    metadata <- metadata[,-c(1,2),drop=FALSE]
  }
  
  # ---End processing of the metadata and tracking

  # initialization of the pcoa and stacked bar plot
  dw <- NULL
  if(length(theTree) != 0){
    the_tree <- read.tree(file=theTree)
    mtx_raw <- relabumtx
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

    if(example_data == "Endesfelder"){
      # return samples as rows
      mtx_raw <- t(mtx_raw)
      mtx_unifracDW <- unifrac2(mtx_raw, mtx_the_tree)
      mtx_unifracDW[is.na(mtx_unifracDW)] = 0
      dw <- mtx_unifracDW
    } else {
      # return samples as rows
      mtx_raw <- t(mtx_raw)
      mtx_unifracDW <- GUniFrac(mtx_raw, mtx_the_tree, alpha=c(0.5,1))$unifracs
      mtx_unifracDW[is.na(mtx_unifracDW)] = 0
      dw <- mtx_unifracDW[,,"d_1"]
    }
  } else {
    theTree <- "no tree"
    set.seed(2)
    mtx_raw <- relabumtx
    dw <- vegdist(mtx_raw, method = "horn")
    dw[is.na(dw)] <- 0
    dw <- as.matrix(dw)
  }

  # calculating the number of cluster and clustering with pam
  if(example_data == "Wegner"){
    num_clustdw <- pamk(dw, krange = 1:6, criterion = "ch")$nc
    pam_commdw <- pam(dw, k = num_clustdw, metric = "euclidean")$cluster
  }else{
    num_clustdw <- pamk(dw, krange = 1:6, criterion = "ch")$nc
    pam_commdw <- pam(dw, k = num_clustdw)$cluster 
  }

  alldata <- new("tabledata", total.ab.matrix=total.ab.matrix, relabumtx=relabumtx, 
                 otuDfRaw=otuDfRaw, taxonomy=taxonomy, otulist.mtx=otulist.mtx, 
                 tracking=tracking, metadata=metadata, theTree=theTree,
                 pam_commdw=pam_commdw, num_clustdw=num_clustdw, dw=dw,
                 pie_mtx=pie_mtx)
  return(alldata)
}

update.otulist <- function(otulist.mtx, relabumtx){
  tmp <- otulist.mtx[which(otulist.mtx[,1] %in% relabumtx), , drop=FALSE]
  otulist.mtx <- tmp
  return(otulist.mtx)
}

get.tax.to.levels <- function(taxonomy){
  tmp <- matrix(NA, nrow = nrow(taxonomy), ncol = 8)
  rownames(tmp) <- rownames(taxonomy)
  colnames(tmp) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")
  
  for(i in 1:nrow(taxonomy)){
    cols <- unlist(strsplit(taxonomy[i,],';'))
    for(j in 1:length(cols)){
      if(cols[j] != ""){
        tmp[i,j] <- cols[j]
      } 
    }
  }
  
  tmp[,8] <- rownames(tmp)

  return(tmp)
}

get.taxonomy.list <- function(taxonomy, level){
  taxlist <- as.matrix(unique(taxonomy[!is.na(taxonomy[, level]), level, drop=FALSE]))
  colnames(taxlist) <- level

  return(taxlist)
}

get.relabumtx <- function(the_raw, taxlist, taxonomy){
  relabumtx <- matrix(0L,nrow = nrow(taxlist), ncol = ncol(the_raw))
  rownames(relabumtx) <- taxlist[,1]
  colnames(relabumtx) <- colnames(the_raw)
  
  for(bac in taxlist[,1]){
    b <- rownames(taxonomy[which(taxonomy[, colnames(taxlist)] %in% bac), ,drop=FALSE])
    if(length(b) >= 1){
      relabumtx[bac, ] <- colSums(the_raw[b,,drop=FALSE])
    }
  }
  
  # transpose the matrix, row-names are the samples, column-names are the bacteria
  relabumtx <- t(relabumtx)

  # make the matrix relative
  mtx2 <- relabumtx[, colSums(relabumtx, na.rm = TRUE)/sum(colSums(relabumtx, na.rm = TRUE))*100 > 0.01,drop=FALSE]
  mtx2 <- mtx2[, colSums(mtx2) != 0, drop=FALSE]
  mtx2 <- sweep(mtx2, MARGIN = 1, FUN = "/", rowSums(mtx2, na.rm = TRUE))
  relabumtx <- mtx2

  return(relabumtx)
}

get.relabumtx.OTU <- function(the_raw){
  relabumtx <- sweep(the_raw, MARGIN = 1, FUN = "/", rowSums(the_raw, na.rm = TRUE))

  return(relabumtx)
}

get.total.abundance.mtx <- function(the_raw, taxlist, taxonomy){
  total.ab.matrix <- matrix(0L,nrow = nrow(taxlist), ncol = ncol(the_raw))
  rownames(total.ab.matrix) <- taxlist[,1]
  colnames(total.ab.matrix) <- colnames(the_raw)
  
  for(bac in taxlist[,1]){
    b <- rownames(taxonomy[which(taxonomy[, colnames(taxlist)] %in% bac), ,drop=FALSE])
    if(length(b) >= 1){
      total.ab.matrix[bac, ] <- colSums(the_raw[b,,drop=FALSE])
    }
  }
  
  # transpose the matrix, row-names are the samples, column-names are the bacteria
  total.ab.matrix <- t(total.ab.matrix)
  mtx2 <- total.ab.matrix[, colSums(total.ab.matrix, na.rm = TRUE)/sum(colSums(total.ab.matrix, na.rm = TRUE))*100 > 0.01,drop=FALSE]
  mtx2 <- mtx2[, colSums(mtx2) != 0, drop=FALSE]
  total.ab.matrix <- mtx2

  return(total.ab.matrix)
}

get.pie.matrix <- function(taxonomy, bacteria, first_level, second_level, label){
  level = switch(first_level, "Kingdom" = c("Kingdom",1), "Phylum" = c("Phylum",2), "Class" = c("Class",3), 
                 "Order" = c("Order",4), "Family" = c("Family",5), "Genus" = c("Genus",6), "Species" = c("Species",7), "OTU" = c("OTU", 8))

  subtax <- taxonomy[which(taxonomy[,first_level] %in% bacteria),,drop=FALSE]
  subtax <- unique(subtax[, 1:level[2]])
  subtax <- subtax[,second_level]
  subtax <- as.matrix(table(subtax))
  total <- colSums(subtax)# tal vez remover el is.na  
  #   cols = value, community, variable, tooltip
  subtax <- cbind(subtax, rep(label, nrow(subtax)), rownames(subtax))#, paste(rownames(subtax), "\n", sprintf((subtax/total)*100, fmt = '%#.3f'),"%"))
  colnames(subtax) <- c("value", "community", "variable")#, "tool_tip")

  return(subtax)
}

get.taxa.in.taxonomy <- function(taxonomy, bacteria, first_level, second_level){
  level = switch(first_level, "Kingdom" = c("Kingdom",1), "Phylum" = c("Phylum",2), "Class" = c("Class",3), 
                 "Order" = c("Order",4), "Family" = c("Family",5), "Genus" = c("Genus",6), "Species" = c("Species",7), "OTU" = c("OTU", 8))
  subtax <- taxonomy[which(taxonomy[,first_level] %in% bacteria),,drop=FALSE]
  
  subtax <- unique(subtax[,second_level])
  return(subtax)
}

get.adjMtx <- function(relabumtx, errthreshold, pvalsccrepe, scoreccrepe, example_data){
  if(example_data == "Wegner"){
    subjs <- 2
  } else {
    subjs <- 10
  }

  ccrepe_out <- NULL
  tryCatch({
    ccrepe_out <- ccrepe(relabumtx, iterations= 1000, min.subj= subjs, sim.score.args= list(method= "spearman", use= "complete.obs"), errthresh= errthreshold)
  }, error=function(e){})
  
  if(is.null(ccrepe_out)){
    subjs <- 5
    ccrepe_out <- ccrepe(relabumtx, iterations= 1000, min.subj= subjs, sim.score.args= list(method= "spearman", use= "complete.obs"), errthresh= errthreshold)
  }
  
  pvals_ccrepe_out <- ccrepe_out$p.values
  pvals_ccrepe_out[is.na(pvals_ccrepe_out)] <- 1
  scores_ccrepe_out <- ccrepe_out$sim.score
  scores_ccrepe_out[is.na(scores_ccrepe_out)] <- 0
  adj_mtx <- matrix(0, nr = nrow(pvals_ccrepe_out), nc = ncol(pvals_ccrepe_out))
  rownames(adj_mtx) <- colnames(adj_mtx) <- rownames(pvals_ccrepe_out)
  adj_mtx[pvals_ccrepe_out < pvalsccrepe & scores_ccrepe_out > scoreccrepe] <- 1
  
  return(adj_mtx)
}

makeEdges <- function(adjmtx){
  edg <- c()
  for(i in 1:ncol(adjmtx)){
    for (j in 1:ncol(adjmtx)){
      if(i!=j){
        opt <- adjmtx[i,j]
        if(opt == 1){
          pair <- c(i, j)
          pair <- sort(pair)
          edg <- rbind(pair, edg)
        }
      }
    }
  }
  
  edg <- as.matrix(unique(edg))
  edg <- edg[order(edg[,1]),]
  rownames(edg) <- seq(1,nrow(edg))
  
  return(edg)
}

get.metadata.items <- function(info) {
  cols <- list()
  for(i in 1:ncol(info)){
    if(suppressWarnings(!is.na(as.numeric(info[1,i])))){
      #is numeric a slider
      cols[i] <-"slider"
    } else {
      #is text a check box
      cols[i] <- "check_box"
    }
  }
  
  chcs_chk <- list()
  sldrs <- list()
  clnnames <- colnames(info)
  
  for(i in 1:length(cols)){
    if(cols[i]=="check_box"){
      chkbx <- unique(info[,i])
      for(chk in chkbx){
        if(!is.na(chk)){
          chcs_chk[chk] <- chk
        }
      }
    }else{
      nm <- clnnames[i]
      values <- as.numeric(info[!is.na(info[,i]),i])
      sldrs[[nm]] <- c(max(values), min(values))
    }
  }
  
  all_elements <- list()
  all_elements[["checkbox"]] <- chcs_chk
  all_elements[["sliders"]] <- sldrs
  
  return(all_elements)
}

get.pcoa <- function (pam_commdw, num_clustdw, tracking, example_data, x, Y = NULL, plot.axes = c(1, 2), dir.axis1 = 1, 
                      dir.axis2 = 1, rn = NULL, main = NULL,...) {
    if (!inherits(x, "pcoa")) 
        stop("Object of class 'pcoa' expected")
    pr.coo <- x$vectors
    if (x$correction[2] > 1) 
        pr.coo <- x$vectors.cor
    k <- ncol(pr.coo)
    if (k < 2) 
        stop("There is a single eigenvalue. No plot can be produced.")
    if (k < plot.axes[1]) 
        stop("Axis", plot.axes[1], "does not exist.")
    if (k < plot.axes[2]) 
        stop("Axis", plot.axes[2], "does not exist.")
    if (!is.null(rn)) 
        rownames(pr.coo) <- rn
    labels = colnames(pr.coo[, plot.axes])
    diag.dir <- diag(c(dir.axis1, dir.axis2))
    pr.coo[, plot.axes] <- pr.coo[, plot.axes] %*% diag.dir
    if (is.null(Y)) {
        print("sin el tatmtx")
        limits <- apply(pr.coo[, plot.axes], 2, range)
        ran.x <- limits[2, 1] - limits[1, 1]
        ran.y <- limits[2, 2] - limits[1, 2]
        xlim <- c((limits[1, 1] - ran.x/10), (limits[2, 1] + 
            ran.x/5))
        ylim <- c((limits[1, 2] - ran.y/10), (limits[2, 2] + 
            ran.y/10))
        par(mai = c(1, 1, 1, 0.5))
        plot(pr.coo[, plot.axes], xlab = labels[1], ylab = labels[2], 
            xlim = xlim, ylim = ylim, asp = 1)
        text(pr.coo[, plot.axes], labels = rownames(pr.coo), 
            pos = 4, cex = 1, offset = 0.5)
        if (is.null(main)) {
            title(main = "PCoA ordination", line = 2)
        }
        else {
            title(main = main, family = "serif", line = 2)
        }
    } else {
        n <- nrow(Y)
        points.stand <- scale(pr.coo[, plot.axes])
        S <- cov(Y, points.stand)
        U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n - 1))^(-0.5))
        colnames(U) <- colnames(pr.coo[, plot.axes])

        # my code
        data_pts <- as.data.frame(pr.coo)
        U <- inrange(pr.coo[,1:2], U)
        ldngs <- as.data.frame(U)
        ldngs["variable"] <- rownames(ldngs)

        data_pts <- data_pts[,1:2]
        data_pts["cluster"] <- rep(0, nrow(data_pts))
        data_pts["alpha_trans"]  <- rep(1.0, nrow(data_pts))

        for(row in 1:num_clustdw){
          data_pts[which(pam_commdw == row), "cluster"] <- paste("Cluster",row)
        }
        data_pts[,c("track_ID", "track_label")] <- NA

        for(row in rownames(data_pts)){
          if(row %in% rownames(tracking)){
            r <- tracking[row,]
            data_pts[row,"track_ID"] <- paste0("track_",gsub(" ","",r[1]))
            if(str_detect(r[2],";")){
              new_label <- str_replace_all(r[2],";","\n")
              data_pts[row,"track_label"] <- new_label
            }else{
              data_pts[row,"track_label"] <- r[2]  
            }
            if(example_data == "Wegner"){
              data_pts[row,"beds"] <- r[1]
            }
          }else{
            data_pts[row,"track_label"] <- paste0(row,"\n","Data not available")
          }

        }
        
        all_DFs <- list(data_pts, ldngs)
        return(all_DFs)
    }
    invisible()
}

get.stackedbarplot <- function(relabumtx, pam_commdw, num_clustdw, top=6) {
  if (length(colnames(relabumtx)) < top) {
    top <- length(colnames(relabumtx))
  }

  top6 <- names(sort(colSums(relabumtx), decreasing = TRUE))[1:top]
  top6.mtx <- relabumtx[, top6]

  all_clrs <- c()
  for(cl in 1:num_clustdw){
    subjects <- names(pam_commdw[which(pam_commdw == cl)])
    submtx <- top6.mtx[subjects, , drop=FALSE]
    submtx <- t(as.matrix(colMeans(submtx)))
    rownames(submtx) <- paste("cluster", cl)
    all_clrs <- rbind(all_clrs, submtx)
  }

  all_clrs <- melt(all_clrs, id.var = "cluster")
  colnames(all_clrs) <- c("cluster","bacteria", "value")
  
  return(all_clrs)
}

load.communities <- function(fncomm, adjmtx){
  comms <- as.matrix(read.table(fncomm, sep="\t", header = FALSE, check.names=FALSE))
  rownames(comms) <- comms[,1]
  comms <- comms[,-1]
  comms <- as.matrix(comms)
  class(comms) <- "numeric"
  colnames(comms) <- c("community")

  idx <- sort(unique(comms[,1]))
  comms <- as.data.frame(comms)
  comms["color"] <- rep("A", nrow(comms))
  comms["gpo"] <- rep("B", nrow(comms))
  
  for(i in idx){
    if(i >= 0){
      comms[which(comms == i), "color"] <- community_colors[i+2]
      comms[which(comms == i), "gpo"] <- paste("community", i+1)
    }else{
      comms[which(comms == i), "color"] <- community_colors[1]
      comms[which(comms == i), "gpo"] <- paste("No community")
    }
  }
  return(comms)
}

get.filter.check.idx.metadata <- function(selected_checks, metadata) {
  if (length(selected_checks) > 1) {
    v <- rownames(which(metadata == selected_checks[1], arr.ind = TRUE))
    for (el in selected_checks) {
      v <- intersect(rownames(which(metadata == el, arr.ind = TRUE)), v)
    }
    
    if (length(v) == 0) {
      for (el in selected_checks) {
        v <- c(rownames(which(metadata == el, arr.ind = TRUE)), v)
      }
    }
    
    return(v)
  } else {
    v <- rownames(which(metadata == selected_checks, arr.ind = TRUE))
    return(v)
  }
}

get.filter.slider.metadata <- function(sliders_choice, metadata, input, tracking){
  nms <- names(sliders_choice)
  nmid2 <- paste("slider", nms[1], sep="_")
  min_max_vls <- input[[nmid2]]
  cl <- metadata[,nms[1],drop=FALSE]
  
  v <- NULL
  subjects <- unique(tracking[rownames(metadata),"track_ID"])
  
  for(subject in subjects){
    sub_idx <- rownames(tracking[tracking[,"track_ID"] %in% subject, ,drop=FALSE])
    #sub <- cl[sub_idx, ,drop=FALSE]
    sub <- cl[which(rownames(cl) %in% sub_idx),,drop=FALSE]
    idx <- sub[which(as.numeric(sub) >= min_max_vls[1] & as.numeric(sub) <= min_max_vls[2]), ,drop=FALSE]
    if(nrow(idx)==0){
      idx <- rownames(sub[which.min(as.numeric(sub)), ,drop=FALSE])
    }else{
      idx <- rownames(idx[which.min(as.numeric(idx)), ,drop=FALSE])
    }
    
    v <- c(idx,v)
  }

  return(v)
}

get.relabumtx.metadata <- function(the.total.filtered, the.total.ab.mtx){
  meta.relabumtx <-the.total.filtered[, colSums(the.total.ab.mtx, na.rm = TRUE)/sum(colSums(the.total.ab.mtx, na.rm = TRUE))*100 > 0.01, drop=FALSE]
  meta.relabumtx <- meta.relabumtx[, colSums(meta.relabumtx) != 0]
  meta.relabumtx <- sweep(meta.relabumtx, MARGIN = 1, FUN = "/", rowSums(meta.relabumtx, na.rm = TRUE))

  return(meta.relabumtx)
}

get.relabumtx.metadata.OTU <- function(the.total.filtered){
  meta.relabumtx <- sweep(the.total.filtered, MARGIN = 1, FUN = "/", rowSums(the.total.filtered, na.rm = TRUE))
  return(meta.relabumtx)
}

get.taxonomy.position <- function(tax){
  level = switch(tax, "Kingdom" = c("Kingdom",1), "Phylum" = c("Phylum",2), "Class" = c("Class",3), 
                 "Order" = c("Order",4), "Family" = c("Family",5), "Genus" = c("Genus",6), "Species" = c("Species",7), "OTU" = c("OTU", 8))
  return(level[2])
}

#---------------paper code---------------

get.stackedbarplot.metadata <- function(relabumtx, datapts, top=6) {
  if (length(colnames(relabumtx)) < top) {
    top <- length(colnames(relabumtx))
  }

  toremove <- "Bifidobacteriaceae"

  top6 <- NULL
  if(toremove %in% colnames(relabumtx)){
    idx.toremove <- which(colnames(relabumtx) == toremove)
    relabumtx_tmp <- relabumtx[, -idx.toremove]
    top6 <- names(sort(colSums(relabumtx_tmp), decreasing = TRUE))[1:top]
  }else{
    top6 <- names(sort(colSums(relabumtx), decreasing = TRUE))[1:top]
  }
  top6.mtx <- relabumtx[, top6]

  all_clrs <- c()

  metalbl <- c("Placebo Formula", "Intervention Formula","Breast Fed")
  for(lbl in metalbl){
    subjects <- rownames(datapts[which(datapts[,'cluster'] == lbl), ])
    submtx <- top6.mtx[subjects, , drop=FALSE]
    submtx <- t(as.matrix(colMeans(submtx)))
    rownames(submtx) <- lbl
    all_clrs <- rbind(all_clrs, submtx)
  }
  
  all_clrs <- melt(all_clrs, id.var = "cluster")
  colnames(all_clrs) <- c("cluster","bacteria", "value")
  
  return(all_clrs)
}

unifrac2 <- function(otu.tab, tree){
  if (!is.rooted(tree))
    stop("Rooted phylogenetic tree required!")
  otu.tab <- as.matrix(otu.tab)
  n <- nrow(otu.tab)
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste("comm", 1:n, sep = "_")
  }
  unifracs <- array(NA, c(n, n), dimnames = list(rownames(otu.tab),
                                                 rownames(otu.tab)))
  for (j in 1:n) {
    unifracs[j, j] <- 0
  }
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU
names\n\t\t\t\t\tin the OTU table
and the tree should match!")
  }
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)
  edge <- tree$edge
  edge2 <- edge[, 2]
  br.len <- tree$edge.length
  cum <- matrix(0, nbr, n)
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]
    node <- edge[tip.loc, 1]
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      cum1 <- cum[, i]
      cum2 <- cum[, j]
      ind <- (cum1 + cum2) != 0
      cum1 <- cum1[ind]
      cum2 <- cum2[ind]
      br.len2 <- br.len[ind]
      diff <- abs(cum1 - cum2)/(cum1 + cum2)
      w <- br.len2 * (cum1 + cum2)
      unifracs[i, j] <- unifracs[j, i] <- sum(diff * w)/sum(w)
    }
  }
  return(unifracs)
}

#-------Tools-------------------
makeCopyMtx <- function(mtx1){
  mtx1 <- as.matrix(mtx1)
  
  mtx2 <- matrix(nrow= nrow(mtx1), ncol= ncol(mtx1), byrow= TRUE)
  for(i in 1:nrow(mtx1)){
    for(j in 1:ncol(mtx1)){
      mtx2[i,j] <- mtx1[i,j]
    }
  }
  
  rownames(mtx2) <- rownames(mtx1)
  colnames(mtx2) <- colnames(mtx1)
  return(mtx2)
}

trim <- function(x) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

ggplotColours <- function(n=6, h=c(0, 360) + 15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

inrange <- function(dots, arrows){
  max_val <- max(abs(dots[,1]))
  min_val <- min(abs(dots[,1]))
  max_data <- max(abs(arrows[,1]))
  min_data <- min(abs(arrows[,1]))
  
  arrows[,1] <- min_val + (max_val - min_val) * ((arrows[,1] - min_data)/(max_data - min_data))
  
  max_val <- max(abs(dots[,2]))
  min_val <- min(abs(dots[,2]))
  max_data <- max(abs(arrows[,2]))
  min_data <- min(abs(arrows[,2]))
  arrows[,2] <- min_val + (max_val - min_val) * ((arrows[,2] - min_data)/(max_data - min_data))
  return(arrows)
}