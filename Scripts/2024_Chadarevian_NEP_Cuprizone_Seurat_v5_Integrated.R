# Detaching all packages
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T,unload = T,force = T))

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(sctransform)
  library(data.table)
  library(ggpubr)
  library(glmGamPoi)
  library(ggrepel)
  library(presto)
  
})
# Setting a random seed to get consistent results across analysis days
set.seed(030989)
# Setting the color palettes for the output files. Tcols currently has colors for up to 10 clusters.
Tcols <- c("#39b600", "#d89000", "#f8766d", "#00b0f6", "#a19f08", "orchid3", "#9590ff", "firebrick1",
           "gray40", 'purple', "darkorange4", "lightskyblue1", "midnightblue")
colfunc <- colorRampPalette(c("midnightblue", "lightskyblue1", "orange", "red"))

# =============================================================================================================
# ====================== Redefined Feature Plot Function ======================================================
# =============================================================================================================
NewFeaturePlot <- function (object, features, dims = c(1, 2), cells = NULL, pt.size = NULL,
                            alpha = 1, order = FALSE, min.cutoff = NA, max.cutoff = NA,
                            reduction = NULL, split.by = NULL, keep.scale = "feature",
                            shape.by = NULL, layer = "data", blend = FALSE, blend.threshold = 0.5,
                            label = FALSE, label.size = 4, label.color = "black", repel = FALSE,
                            ncol = NULL, coord.fixed = FALSE, by.col = TRUE,
                            sort.cell = lifecycle::deprecated(), interactive = FALSE, combine = TRUE,
                            raster = NULL, raster.dpi = c(512, 512),
                            cols = if (blend) {c("lightgrey", "#ff0000", "#00ff00")}
                            else {c("lightgrey", "blue")}) 
{
  if (lifecycle::is_present(arg = sort.cell)) {
    deprecate_stop(when = "4.9.0", what = "FeaturePlot(sort.cell = )",
                   with = "FeaturePlot(order = )")
  }
  if (isTRUE(x = interactive)) {
    return(IFeaturePlot(object = object, feature = features[1],
                        dims = dims, reduction = reduction, layer = layer))
  }
  if (!is.null(x = keep.scale)) {
    keep.scale <- rlang::arg_match0(arg = keep.scale, values = c("feature",
                                                                 "all"))
  }
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
                                                                                           size = 14, margin = margin(r = 7)))
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (!rlang::is_integerish(x = dims, n = 2L, finite = TRUE) && !all(dims > 0L)) {
    rlang::abort(message = "'dims' must be a two-length integer vector")
  }
  if (isTRUE(x = blend) && length(x = features) != 2) {
    rlang::abort(message = "Blending feature plots only works with two features")
  }
  if (isTRUE(x = blend)) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = {
                     warn(message = "No colors provided, using default colors")
                     default.colors
                   }, `1` = {
                     warn(message = paste("Only one color provided, assuming",
                                          sQuote(x = cols), "is double-negative and augmenting with default colors"))
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warn(message = paste("Only two colors provided, assuming specified are for features and agumenting with",
                                          sQuote(default.colors[1]), "for double-negatives",
                     ))
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warn(message = "More than three colors provided, using only first three")
                     cols[1:3]
                   })
  }
  if (isTRUE(x = blend) && length(x = cols) != 3) {
    rlang::abort("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% Cells(x = object[[reduction]])
  data <- FetchData(object = object, vars = c(dims, "ident",
                                              features), cells = cells, layer = layer)
  if (ncol(x = data) < 4) {
    rlang::abort(message = paste("None of the requested features were found:",
                                 paste(features, collapse = ", "), "in layer ", layer))
  }
  else if (!all(dims %in% colnames(x = data))) {
    rlang::abort(message = "The dimensions requested were not found")
  }
  features <- setdiff(x = names(x = data), y = c(dims, "ident"))
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[,
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[,
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    rlang::abort(message = "There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  names(x = min.cutoff) <- names(x = max.cutoff) <- features
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
  ]$maxcolors, no = length(x = cols))
  for (i in seq_along(along.with = features)) {
    f <- features[i]
    data.feature <- data[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[f], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[f], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      }
      else {
        as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                                         breaks = brewer.gran)))
      }
    }
    data[[f]] <- data.feature
  }
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  }
  else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells,
                                                            drop = TRUE], object[[split.by, drop = TRUE]][cells,
                                                                                                          drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[,
                                                                   dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[,
                                                               dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3],
                                col.threshold = blend.threshold, negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ],
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident,
                      , drop = FALSE]
    if (isTRUE(x = blend)) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[,
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        rlang::abort(message = paste("The following features have no value:",
                                     paste(no.expression, collapse = ", ")))
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")],
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (isTRUE(x = blend)) {
        cols.use <- as.numeric(x = as.character(x = data.plot[,
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      }
      else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature,
                                   shape.by)]
      plot <- SingleDimPlot(data = data.single, dims = dims,
                            col.by = feature, order = order, pt.size = pt.size,
                            alpha = alpha, cols = cols.use, shape.by = shape.by,
                            label = FALSE, raster = raster, raster.dpi = raster.dpi) +
        scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
        theme_cowplot() + CenterTitle()
      if (isTRUE(x = label)) {
        plot <- LabelClusters(plot = plot, id = "ident",
                              repel = repel, size = label.size, color = label.color)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA,
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident),
                                                                    limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(),
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                               axis.title.x = element_blank())
        }
      }
      else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        }
        else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warn(message = paste0("All cells have the same value (",
                                  unique.feature.exp, ") of ", dQuote(x = feature)))
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad,
                                                                       guide = "colorbar"))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale ==
          "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot & scale_color_gradientn(colors = cols,
                                                              limits = c(min.feature.value, max.feature.value)))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (isTRUE(x = blend)) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots,
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")),
                                                                                              expand = c(0, 0)) + labs(x = features[1],
                                                                                                                       y = features[2], title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || isTRUE(x = blend),
                 yes = ncol, no = length(x = features))
  legend <- if (isTRUE(x = blend)) {
    "none"
  }
  else {
    split.by %iff% "none"
  }
  if (isTRUE(x = combine)) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""),
                                                                   limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]),
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) ==
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      }
      else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots),
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- patchwork::wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      plots <- patchwork::wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                                       length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" &&
        !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots & scale_color_gradientn(colors = cols,
                                                              limits = c(min.feature.value, max.feature.value)))
    }
  }
  return(plots)
}

scVolc <- function (obj = NULL, metaName = NULL, ID1 = NULL, ID2 = NULL,
                    assay = "RNA", LFC = 0.3219, FDR = 0.05, numLabs = NULL, labList = NULL)
{
  if (is.null(obj)) {
    stop("You must supply a Seurat object.")
  }
  if (is.null(metaName)) {
    stop("You must supply a metadata slot.")
  }
  if (is.null(ID1) || is.null(ID2)) {
    stop("You must supply both identities that you want to compare.")
  }
  if (is.null(obj@meta.data[[metaName]])) {
    stop(paste0(metaName, " is not present in the supplied object."))
  }
  if (!(dir.exists(file.path(outDir, "Volcano")))) {
    dir.create(file.path(outDir, "Volcano"))
    dir.create(file.path(outDir, "Volcano", "Volc_Genes"))
  }
  Idents(obj) <- metaName
  DefaultAssay(obj) <- assay
  markers <- FindMarkers(object = obj,
                         ident.1 = ID1,
                         ident.2 = ID2,
                         logfc.threshold = -Inf,
                         min.pct = 0.1,
                         only.pos = F)
  markers$Gene <- rownames(markers)
  markers <- markers[, c(6, 1:5)]
  markers <- markers[order(markers$avg_log2FC, decreasing = T),
  ]
  ID1 <- I(unlist(lapply(list(ID1), paste, collapse = "")))
  ID2 <- I(unlist(lapply(list(ID2), paste, collapse = "")))
  
  write.table(markers, file = file.path(outDir, "Volcano", "Volc_Genes",
                                        paste0("VolcGenes_", Exp, "_", ID, "_", ID1,
                                               "vs", ID2, "_Assay_", assay, ".tsv")), quote = F,
              sep = "\t", row.names = F)
  
  
  # Replacing any padj=0 with a value one order of magnitude lower than the lowest padj
  lowP <- min(subset(markers$p_val_adj, -log10(markers$p_val_adj) < Inf))
  markers$p_val_adj[markers$p_val_adj==0] <- as.numeric(lowP*0.1)
  markers$PercDiff <- abs(markers$pct.1 - markers$pct.2)
  # Apply an FDR cutoff to the results
  pcutoff_res <- markers[(markers$p_val_adj<FDR),]
  # Apply a log-fold change cutoff to the results
  folded_res <- pcutoff_res[(pcutoff_res$avg_log2FC >= LFC | pcutoff_res$avg_log2FC <= -LFC),]
  
  # Calculating the ideal limit values for the x and y axes
  logMin <- min(markers$avg_log2FC)
  logMax <- max(markers$avg_log2FC)
  pMax <- max(-log10(markers$p_val_adj))
  
  print(paste0("logMin: ", logMin))
  print(paste0("logMax: ", logMax))
  print(paste0("pMax: ", pMax))
  
  logRound <- 1
  pRound <- 10
  
  xMin <- (logRound*round((logMin-(logRound/2))/logRound))
  xMax <- (logRound*round((logMax+(logRound/2))/logRound))
  xLim <- max(abs(xMin), xMax)
  xDivide <- 1
  xExpand <- xLim * 0.05
  yMax <- (pRound*round((pMax+(pRound/2))/pRound))
  #yDivide <- round(yMax/10)
  yDivide <- 20
  yExpand <- yMax * 0.015
  
  if (TRUE) {
    # Generating and saving the volcano plot
    mutateddf <- mutate(as.data.frame(markers),
                        Cutoff=rownames(markers) %in% rownames(folded_res),
                        Color=ifelse(markers$avg_log2FC>0, '#AD0505', '#0C0C80'),
                        Combo=markers$avg_log2FC*-log10(markers$p_val_adj)) #Will have different colors depending on significance
    mutateddf <- within(mutateddf, Color[Cutoff == "FALSE"] <- 'gray80')
    rownames(mutateddf) <- rownames(markers)
    
    input <- cbind(gene=rownames(mutateddf), mutateddf) #convert the rownames to a column
    input <- input[order(input$Combo, decreasing = T),]
    
    geneLabels <- subset(input, input$Cutoff == "TRUE")
    geneUp <- subset(geneLabels, geneLabels$avg_log2FC > 0) # Keeping only positive LFC genes
    geneDown <- subset(geneLabels, geneLabels$avg_log2FC < 0) # Keeping only negative LFC genes
    
    # Setting the color scale based on which subsets of genes are present
    if (nrow(geneUp)== 0 && nrow(geneDown) == 0){
      cols <- c("gray80")
    } else if (nrow(geneUp) == 0 && nrow(geneDown) > 0){
      cols <- c("#0C0C80", "gray80")
    } else if (nrow(geneUp) >0 && nrow(geneDown) == 0) {
      cols <- c("#AD0505", "gray80")
    } else {
      cols <- c("#0C0C80", "#AD0505", "gray80")
    }
    
    if (is.null(labList)){
      # Checking if there are enough positive LFC genes to label the number set
      # by numLabs. If not, using all genes in the list
      if (nrow(geneUp) >= numLabs) {
        geneUp <- head(geneUp, numLabs)
      }
      
      # Checking if there are enough negative LFC genes to label the number set
      # by numLabs. If not, using all genes in the list
      if (nrow(geneDown) >= numLabs){
        geneDown <- tail(geneDown, numLabs)
      }
      
      finLabs <- rbind(geneUp, geneDown)
      
    } else {
      
      finLabs <- subset(geneLabels, geneLabels$gene %in% labList)
      
    }
    
    volc <- ggplot(input, aes(avg_log2FC, -log10(p_val_adj))) + #volcano plot with log2Foldchange versus pvalue
      theme_classic() +
      theme(axis.line = element_line(linewidth=1),
            axis.text = element_text(color = "black", size = 18, face = "bold"),
            axis.ticks = element_line(color = "black", linewidth = 1),
            axis.title = element_text(face = 'bold', size = 18),
            plot.title = element_text(face = 'bold'),
            legend.position = "NA") +
      geom_point(data=input, aes(col=Color, size=PercDiff)) + #add points colored by significance
      scale_color_manual(values = cols) +
      xlim(-xLim, xLim) +
      ylim(0, yMax) +
      scale_x_continuous(limits = c(-xLim, xLim), breaks = seq(-xLim, xLim, by=xDivide), expand = c(0, xExpand)) +
      scale_y_continuous(limits = c(0, yMax), breaks = seq(0, yMax, by=yDivide), expand = c(0, yExpand)) +
      geom_hline(yintercept = -log10(FDR), linetype = "dashed", linewidth = 1, color = "gray30") +
      geom_vline(xintercept = c(LFC, -LFC), linetype = "dashed", linewidth = 1, color = "gray30") +
      #ggtitle(names(comparisons)[k]) + #e.g. 'Volcanoplot DESeq2'
      geom_label_repel(data=finLabs, color="white", size=8, segment.color = 'black',
                       force = 200, aes(label=gene, fill = alpha(c(Color), 0.7))) +
      #geom_label_repel(data=geneDown, color="white", size=8, segment.color = 'black',
      #                aes(label=gene, fill = alpha(c(Color), 0.7))) +
      scale_fill_identity()
    
    # Saving the volcano plot
    png(filename = file.path(outDir, "Volcano",
                              paste0("Volcano_", Exp, "_", ID, "_", ID1, "vs", ID2,"_FDR", FDR, "_LFC", LFC, "_Assay_", assay, ".png")),
         res=300, units='in', width=6, height=6)
    
    plot(volc)
    
    dev.off()
    
  }
}

ClusterDistributions <- function (obj = NULL, metaName = NULL, assay = NULL, cluster.names = NULL) 
{
  if (is.null(levels(obj@meta.data[[metaName]]))) {
    stop(paste0("The data in the ", metaName, " slot must stored as a factor."))
  }
  else if (length(levels(obj@meta.data[[metaName]])) >= 1) {
    maxClust <- (max(as.integer(obj[[paste0(cluster.names)]][, 1])) - 1)
    row_base <- c("Total_Cells", "Percent_of_Cluster", "Percent_of_Total")
    rows <- as.list(sapply(1:length(levels(obj@meta.data[[metaName]])), 
                           function(x) {
                             paste(levels(obj@meta.data[[metaName]])[x], 
                                   row_base, sep = "_")
                           }))
    cellDist <- data.frame(matrix(ncol = 1, nrow = ((3 * 
                                                      length(levels(obj@meta.data[[metaName]])))+1), dimnames = list(c(rows, "FoldChange"),
                                                                                                                  c("Rows"))))
    cellDist[, 1] <- rownames(cellDist)
    for (i in levels(Idents(obj))) {
      clusterIdents <- subset(obj@meta.data[[metaName]], 
                              rownames(obj@meta.data) %in% WhichCells(obj, 
                                                                      ident = i))
      clustSize <- length(clusterIdents)
      print("-----------------------------------------------")
      print(paste0(i))
      for (j in 1:length(levels(obj@meta.data[[metaName]]))) {
        metaIdent <- levels(obj@meta.data[[metaName]])[[j]]
        totIdent <- length(grep(paste0("^", metaIdent, 
                                       "$"), obj@meta.data[[metaName]]))
        numClust <- length(grep(paste0("^", metaIdent, 
                                       "$"), clusterIdents))
        percClust <- format(round((numClust/clustSize) * 
                                    100, 1), nsmall = 1)
        percTot <- format(round(numClust/(totIdent) * 
                                  100, 1), nsmall = 1)
        cellDist[((j * 3) - 2), (i)] <- numClust
        cellDist[((j * 3) - 1), (i)] <- percClust
        cellDist[((j * 3)), (i)] <- percTot
        print(paste("There are", numClust, metaIdent, 
                    "cells in the", i, "cluster", sep = " "))
        print(paste(percClust, "percent of the", i, 
                    "cluster is", metaIdent, "cells", sep = " "))
        print(paste("The", i, "cluster contains", percTot, 
                    "percent of the total", metaIdent, "Cells", 
                    sep = " "))
      }
    }
    
    for (ii in 1:(ncol(cellDist)-1)){
      cellDist[nrow(cellDist),ii+1] <- (as.numeric(cellDist[6,ii+1])-as.numeric(cellDist[3,ii+1]))/as.numeric(cellDist[3,ii+1])
    }
    
    write.table(cellDist, file.path(outDir, paste0(Exp, 
                                                   "_", ID, "_", metaName, "_ClusterStats.tsv")), sep = "\t", 
                quote = F, row.names = F)
    return(cellDist)
  }
  else {
    stop(paste0("The metadata slot ", metaName, " does not exists in the object."))
    return(NULL)
  }
}

HorizontalBarPlot <- function(dat=NULL, cols=NULL, obj=NULL, metaName=NULL)
{
  # Calculates the required axis height to the nearest 5 percent
  # x is the melted data
  maxAxis <- function(x) {
    highVal <- floor(max(x$Percent)/0.05)
    maxWhole <- (highVal * 5) + 5
    maxDec <- maxWhole / 100
    return(maxDec)
  }
  # Keep only the rows and columns containing the cluster percentages
  percs <- dat[seq(3, nrow(dat), by=3), 2:ncol(dat)]
  # Convert the values in percs to numerics and then percentage decimals
  percs <- mutate_all(percs, function(x) as.numeric(x)/100)
  # Set the identities that will classify each sample
  identities <- levels(obj@meta.data[[metaName]])
  # Add the identities to the percs file
  percs$ID <- identities
  # Melt the data into long format
  percs2 <- data.table::melt(data = as.data.table(x = percs), 
                             id.vars = "ID",
                             variable.name = "Cluster", 
                             value.name = c("Percent"))
  
  # Matches the factor order to the order of percs$Genotype input above
  # This will arrange the bars in the same order as percs$Genotype
  percs2$ID <- factor(x = percs2$ID, levels = unique(identities))
  
  # Calculating saved image dimensions to be proportionate to axis size
  # Width = 0.269in per 5 percent
  # Height = 0.4in per cluster
  imageWidth <- 3
  imageHeight <- ((ncol(percs)-1) * 0.35)
  
  ggplot(percs2, aes(x=Cluster, y=Percent, fill = ID)) + 
    theme_linedraw() +
    geom_bar(data=percs2,
             stat = 'identity',
             width = 0.9,
             position = position_dodge2(width = 0.9, padding = 0, reverse = T),
             group= 'ID') +
    scale_fill_manual(values = cols) +
    theme(legend.position = c(0.8, 0.22), 
          legend.title = element_text(face = 'bold'),
          legend.text = element_text(face = 'bold'),
          legend.key.size = unit(1, "line")) + 
    theme(axis.line = element_line(size=0.5, colour = "black"),
          axis.ticks = element_line(size=0.5, colour = "black"),
          plot.margin = margin(r=0.1, b=0.1, unit = 'in'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          panel.background = element_blank()) + 
    xlab(NULL) + 
    ylab("Percent (%)") +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.x=element_text(colour = "black", size = 8, face = "bold", angle = 45, hjust = 0.05),
          axis.text.y=element_text(colour="black", size = 10, face = 'bold'),
          axis.title.x = element_text(face = 'bold', size = 10)) + 
    theme(axis.ticks.y = element_blank()) +
    scale_y_continuous(limits=c(0, maxAxis(percs2)), 
                       breaks=seq(0.00, maxAxis(percs2), by=0.05), 
                       labels = scales::percent_format(accuracy = 1L, suffix = ""), 
                       expand = c(0,0), 
                       position = "right") +
    scale_x_discrete(expand = c(0,0), limits = rev(levels(percs2$Cluster))) +
    coord_flip()
  
  ggsave(filename=file.path(outDir, paste0("BarPlot_", Exp, "_", ID, "_", metaName, "_ClustersStats.png")),
         width=(imageWidth+1.25), height=(imageHeight+0.5), units="in", dpi=300)
  
}

UpdateClusterNames <- function(obj=NULL, idents=NULL)
{
  # Only run if the correct number of identities are supplied
  if (length(idents) == length(levels(Idents(obj)))){
    # Add the names of the current identities to the idents object
    names(idents) <- levels(obj)
    # Use the idents list to rename the object IDs
    obj <- RenameIdents(obj, idents)
    # Print the updated ID levels for the obj
    levels(obj)
    # Save the cluster IDs to a new metadata slot
    obj$Clusters <- Idents(obj)
    
  } else {
    # If the list of identities is not the same length as the current object identities
    print(paste0("Wrong number of identities supplied. Current object has ", length(levels(Idents(obj)))), " identities.")  
    
  }
  
  return(obj)  
}

CreateIntegratedDirectories <- function(outdir=NULL)
{
  # Create the Seurat output directory if it doesn't exists
  dir.create(file.path(outDir), showWarnings = F, recursive = T)
  # Creating the necessary Contamination and Preliminary analysis output directories
  # Warnings are silenced to quiet the process in case directories already exist
  for(dir in c("Barcodes", "Gene_Plots", "Gene_Expression", "Signature_Plots")){
    # Creating the above directories in the base.path location
    dir.create(file.path(outDir, dir), recursive = T, showWarnings = F)
  }
  for(dir in c("UMAP", "Violin")){
    # Creating the above directories in the Gene_Plots location
    dir.create(file.path(outDir, "Gene_Plots", dir), showWarnings = F)
  }
  
}

CreateHeatmaps <- function(obj, geneList, assay = 'RNA', metaName=NULL)
{
  # Set the default assay to the specified assay
  DefaultAssay(obj) <- assay
  # Set the identities to the specified metadata slot
  Idents(obj) <- metaName
  # Make sure the identities provided match the identities on the gene list
  if (!(identical(levels(Idents(obj)), levels(geneList[['cluster']])))){
    stop(paste0("Identities in ", metaName, " do not match identities in the gene list."))
    }
  # Colors to use for the cluster annotation in DoHeatmap
  Tcols <- c("#39b600", "#d89000", "#f8766d", "#00b0f6", "#a19f08", "orchid3", "#9590ff", "firebrick1", "gray40", 'purple')
  # Create Seurat heatmap using the DoHeatmap function
  # If the object contains more than 30,000 cells, it will be downsampled to 20,000 cells to allow the 
  # heatmap to print properly
  if (ncol(obj) >= 30000){
    DoHeatmap(object = subset(obj, downsample = 20000), features = geneList$gene, group.colors = Tcols) + NoLegend()
  } else {
    DoHeatmap(object = obj, features = geneList$gene, group.colors = Tcols) + NoLegend()
  }
  # Save the heatmap
  ggsave(filename=file.path(outDir, paste0("Heatmap_", Exp, "_", ID, ".png")), 
         width=30, height=12, dpi=300, units='in')
  # Create a dot plot of the supplied genes
  DotPlot(obj, 
          features = unique(geneList$gene),
          idents = levels(Idents(obj)),
          cols = c("gray100", "darkgreen")) +
    theme_classic() +
    theme(axis.line = element_line(size=1), 
          axis.text = element_text(color = "black", size = 11, face = "bold"),
          axis.ticks = element_line(color = "black", size = 0),
          axis.title = element_text(face = 'bold', size=11),
          plot.title = element_text(face = 'bold'),
          legend.position = "right",
          legend.text=element_text(face = 'bold', size = 11)) +
    guides(x = guide_axis(angle = 45)) +
    ylab(NULL) + xlab(NULL)
  # Save the dot plot
  ggsave(filename=file.path(outDir, paste0("DotPlot_", Exp, "_", ID, ".png")), 
         width=(3.2+(0.185*length(unique(geneList$gene)))), height=5, dpi=300, units='in')
}

CompareIdentMarkers <- function(obj=NULL, metaName=NULL, ID1=NULL, ID2=NULL, assay='RNA', log.thresh=0.1, FDR=0.05)
{
  # Check the all necessary variables exist
  if (is.null(obj)){
    stop("You must supply a Seurat object.")
  }
  if (is.null(metaName)){
    stop("You must supply a metadata slot.")
  }  
  if (is.null(ID1) || is.null(ID2)){
    stop("You must supply both identities that you want to compare.")
  }
  if (is.null(obj@meta.data[[metaName]])){
    stop(paste0(metaName, " is not present in the supplied object."))
  }
  # If an output directory for the metadata comparisons doesn't exist, create it
  if (!(dir.exists(file.path(outDir, "Gene_Expression", metaName)))){
    dir.create(file.path(outDir, "Gene_Expression", metaName))
  }
  # Set the object identities to the specified metadata slot
  Idents(obj) <- metaName
  # Set the default assay to the input assay
  DefaultAssay(obj) <- assay
  # Find the markers between the supplied IDs
  markers <- FindMarkers(object = obj, ident.1 = ID1, ident.2 = ID2,
                         logfc.threshold = log.thresh, min.pct = 0.1, only.pos = F)
  # Subset to only keep the genes that meet the FDR cutoff
  markers <- subset(markers, markers$p_val_adj <= FDR)
  # Create a column for the gene names
  markers$Gene <- rownames(markers)
  # Reorder the columns so the gene names come first
  markers <- markers[,c(6,1:5)]
  # Order the log2fc column in decreasing order
  markers <- markers[order(markers$avg_log2FC, decreasing=T),]
  # Collapse the ID1 variable in case it is longer than one item
  ID1 <- I(unlist(lapply(list(ID1),paste,collapse="")))
  # Collapse the ID2 variable in case it is longer than one item
  ID2 <- I(unlist(lapply(list(ID2),paste,collapse="")))
  # Save the gene expression table
  write.table(markers, 
              file= file.path(outDir, "Gene_Expression", metaName,
                              paste0("SigGenes_", Exp, "_", ID, "_", ID1, "vs", ID2, "_Assay_", assay, ".tsv")),
              quote=F, sep='\t', row.names = F)
  
  return(markers)
}

# ====================================== Read in previous RDS object ===========================================
# User inputs the basic  the location and sample ID for the file being processed
storageDir <- "/path/to/main/storage/directory"
study <- "2024_Chadarevian_NEP"
Exp <- "Cuprizone"
Analysis <- "Preliminary"
samps <- c("Sublibrary1", "Sublibrary2")
ID <- "Sublibrary_Combined_CPZ"

# Generate a pointer to the input directory of the analysis
inDir <- file.path(storageDir, study, Exp, "Seurat", Analysis)

# Generate the pointer to the output directory for the analysis
outDir <- file.path(storageDir, study, Exp, "Seurat", paste0(day, "_Analysis"), ID)

# Create the necessary output directories
CreateIntegratedDirectories(outDir)

# Import and merge sublibraries
objList <- list()
for (i in 1:length(samps)){
  # Read in the object
  objList[[i]] <- readRDS(file.path(inDir, samps[i],
                                    paste0(Exp, "_", samps[i], "_Seurat_Object.rds")))

  # Set the default assay to Spatial
  DefaultAssay(objList[[i]]) <- "RNA"
  # Remove the data and scale.data slots
  objList[[i]]@assays$RNA$data <- NULL
  objList[[i]]@assays$RNA$scale.data <- NULL
}

# Merge the objects
scRNA <- merge(objList[[1]], objList[[2]])

# Create a new metadata column of merged values
scRNA$Treatment <- factor(scRNA$Treatment, levels = c("Normal_Chow", "Cuprizone"))

# Saving the Seurat object
saveRDS(scRNA, file = file.path(outDir, paste0(Exp, "_", ID, "_Seurat_Object.rds")))

# Remove unneeded variables and release memory
rm(objList)
gc()

# =====================================================================================================================
# ============================================ Normalization and Scaling ==============================================
# =====================================================================================================================
# Normalize the data
scRNA <- NormalizeData(scRNA)

# Calculate the percentage of mitochondrial reads
scRNA <- PercentageFeatureSet(scRNA, pattern = '^MT-', col.name = "percent.mito")

# Find variable features
scRNA <- FindVariableFeatures(scRNA, assay = "RNA", nfeatures = 3000)

# Update variable features to not include mitochondrial genes or ribosomal genes for PCA
PCAFeats <- data.frame(Gene = VariableFeatures(scRNA))
PCAFeats <- subset(PCAFeats, !(PCAFeats$Gene %in% grep("^MT-", PCAFeats$Gene, value=T)))
PCAFeats <- subset(PCAFeats, !(PCAFeats$Gene %in% grep("^RP[SL]", PCAFeats$Gene, value=T)))

# Set the variables that will be regressed out
regVars <- c("percent.mito", "percent.ribo", "nCount_RNA")

scRNA <- ScaleData(scRNA, 
                   vars.to.regress = regVars)
scRNA <- RunPCA(scRNA,
                assay = "RNA",
                features = PCAFeats$Gene, npcs = 100)

# Printing an elbow plot for PC selection
ElbowPlot(scRNA, reduction = "pca", ndims = 50)

# Save the elbow plot
ggsave(filename= file.path(outDir, paste0("ElbowPlot_", Exp, "_", ID, ".png")),
       width=8, height=5, units="in", dpi=300, bg="white")

# Integrate and join layers
scRNA <- IntegrateLayers(scRNA, method = RPCAIntegration, orig.reduction = "pca", 
                         new.reduction = "integrated.rpca", verbose = F)
scRNA[["RNA"]] <- JoinLayers(scRNA[["RNA"]])

# Saving the Seurat object
saveRDS(scRNA, file = file.path(outDir, paste0(Exp, "_", ID, "_Seurat_Object_PostIntegration.rds")))

# ===========================================================================================================
# ======================================== Clustering the cells =============================================
# ===========================================================================================================

# Set the PCs and resolution for the clustering
PCs <- c(1:14)
Res <- 0.25 

# Nearest neighbor analysis, clustering, and dimensional reduction
scRNA <- FindNeighbors(scRNA, dims = PCs, reduction = "integrated.rpca")
scRNA <- FindClusters(scRNA, resolution = Res, cluster.name = "rpca_clusters")
scRNA <- RunUMAP(scRNA, dims = PCs, reduction = "integrated.rpca", reduction.name = "umap.rpca")

# Generating a plot
DimPlot(scRNA,
        reduction = "umap.rpca",
        label = T,
        pt.size = 1,
        group.by = c("Treatment", "rpca_clusters"))

# =======================================================================================================
# ===================================== Differential Gene Expression ====================================
# =======================================================================================================

# Set the assay to perform analysis
DefaultAssay(scRNA) <- 'RNA'
AnalysisAssay <- 'RNA'

# Set the metadata column for barplot
metaName <- "Treatment"
scRNA$Treatment <- factor(scRNA$Treatment, levels = c("Normal_Chow", "Cuprizone"))

# Create a list of new cluster IDs to be added to the object
new.cluster.ids <- c("Homeostatic", "MHCII", "RhoGTPase", "IFN")

# This function takes a Seurat object and a list of cluster names as input and replaces
# the current identity classes with the names provided by new.cluster.ids. The new names
# will be stored in a metadata column titled 'Clusters' for easy reference.
scRNA <- UpdateClusterNames(obj=scRNA,
                            idents=new.cluster.ids)

# Find markers for every cluster compared to all remaining cells
scRNA.markers <- FindAllMarkers(object = scRNA,
                                assay = AnalysisAssay,
                                slot = "data",
                                logfc.threshold = 0.1,
                                only.pos = F)

# Subset the list to only contain significant gene hits and reorder the columns
scRNA.markers <- subset(scRNA.markers, scRNA.markers$p_val_adj <= 0.01)
scRNA.markers <- scRNA.markers[, c(7,1:6)]

# Sort by decreasing log2FC
scRNA.markers <- scRNA.markers[order(scRNA.markers$avg_log2FC, decreasing=T),]

# Reorder by cluster
scRNA.markers <- scRNA.markers[order(scRNA.markers$cluster),]

# Save the DEG data
write.table(scRNA.markers,
            file = file.path(outDir,
                             "Gene_Expression",
                             paste0("SigGenes_", Exp, "_", ID, "_Assay_", AnalysisAssay,
                                    "_AllGroups.tsv")),
            quote = F,
            sep = "\t",
            row.names = F)



# Collect the top 10 markers (or all markers if less than 10) for each cluster.
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top10$cluster <- factor(top10$cluster, levels = new.cluster.ids)

# This function generates and saves a traditional Seurat heatmap and a dotplot version of the genes supplied
CreateHeatmaps(scRNA,
               top10,
               assay = AnalysisAssay,
               metaName = 'Clusters')

# Saving the Seurat object
saveRDS(scRNA,
        file = file.path(outDir, paste0(Exp, "_", ID, "_Seurat_Object_PostClustering.rds")))

# Measure the cluster cell distributions for the identities in the metadata slot
cellDists <- ClusterDistributions(scRNA,
                                  metaName = metaName,
                                  assay = AnalysisAssay,
                                  cluster.names = "Clusters")

# Set the colors for the barplot
treatCols <- c( '#0C0C80', '#AD0505')

# Create a barplot of the cluster cell distributions
HorizontalBarPlot(dat = cellDists,
                  cols = treatCols,
                  obj = scRNA,
                  metaName = metaName)

# Set the IDs to compare gene expression for
# Results will be ID1 vs ID2
ID1 <- c("Cuprizone")
ID2 <- c("Normal_Chow")

# Set additional variable values for comparison
metaName <- "Treatment"
AnalysisAssay <- "RNA"

# Compare the gene expression in the above clusters
comp.markers <- CompareIdentMarkers(scRNA,
                                    metaName = metaName,
                                    ID1 = ID1,
                                    ID2 = ID2,
                                    log.thresh = 0.01,
                                    FDR = 0.05,
                                    assay = AnalysisAssay)

# Create a volcano plot for the above comparison
geneList <- c("IFI44L", "IFI6", "MX1", "OAS2", "IFIT1", "EYA2", "IGF2BP2")

scVolc(scRNA,
       metaName = "Treatment",
       ID1 = ID1,
       ID2 = ID2,
       LFC = 1,
       assay = AnalysisAssay,
       numLabs = 5,
       labList = NULL)

# ============================================================================================================
# ==================================== Save Publication Quality Plots ========================================
# ============================================================================================================

# Split the plots by a specified metadata slot
metaName <- 'Treatment'

# ========================================== UMAP Clusters ===================================================
# UMAP plot with labels
DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap.rpca') +
  theme_classic() +
  scale_color_manual(values=Tcols) +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks = element_line(color='black', linewidth=3),
        axis.ticks.length = unit(0, 'cm'),
        axis.text = element_blank(),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename=file.path(outDir, paste0("UMAP_", Exp, "_", ID, "_ClustersLabeled.png")),
       width=5, height=5, units="in", dpi=300)

# UMAP plot without labels
DimPlot(object = scRNA, pt.size = 1, label = F, reduction = 'umap.rpca') + theme_classic() +
  scale_color_manual(values=Tcols) +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks = element_line(color='black', linewidth=3),
        axis.ticks.length = unit(0, 'cm'),
        axis.text = element_blank(),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename= file.path(outDir, paste0("UMAP_", Exp, "_", ID, ".png")),
       width=5, height=5, units="in", dpi=300)

if (length(unique(scRNA@meta.data[[metaName]])) > 1){
  
  DimPlot(object = scRNA, split.by=metaName,  pt.size = 1, label = T, reduction = 'umap.rpca') +
    theme_classic() +
    scale_color_manual(values = Tcols) +
    theme(axis.line=element_line(linewidth=1),
          axis.ticks = element_line(color='black', linewidth=3),
          axis.ticks.length = unit(0, 'cm'),
          axis.text = element_blank(),
          text = element_text(face='bold', color ='black', size=24),
          panel.border = element_blank(),
          plot.margin = unit(c(0,1,0,0), "cm"),
          legend.position='none') + ggtitle(NULL) +
    xlab("UMAP 1") + ylab("UMAP 2")
  
  ggsave(filename = file.path(outDir, paste0("UMAP_", Exp, "_", ID, "_SplitBy_", metaName, ".png")),
         width=15, height=5, units="in", dpi=300)
  
  for (i in unique(scRNA@meta.data[[metaName]])) {
    
    temp <- scRNA
    Idents(temp) <- temp@meta.data[[metaName]]
    sub_obj <- subset(temp, idents = i)
    Idents(sub_obj) <- sub_obj@meta.data[["Clusters"]]
    
    DimPlot(object = sub_obj,  pt.size = 1, label = F, reduction = 'umap.rpca') +
      theme_classic() +
      scale_color_manual(values = Tcols) +
      theme(axis.line=element_line(linewidth=1),
            axis.ticks = element_line(color='black', linewidth=3),
            axis.ticks.length = unit(0, 'cm'),
            axis.text = element_blank(),
            text = element_text(face='bold', color ='black', size=24),
            panel.border = element_blank(),
            plot.margin = unit(c(0,1,0,0), "cm"),
            legend.position='none') + ggtitle(NULL) +
      xlab("UMAP 1") + ylab("UMAP 2")
    
    ggsave(filename = file.path(outDir, paste0("UMAP_", Exp, "_", ID, "_Only_", i, ".png")),
           width=5, height=5, units="in", dpi=300)
    
    rm(temp, sub_obj)
    gc()
    
  }
}

# =============================== Individual Genes =============================================
AnalysisAssay <- 'RNA'
DefaultAssay(scRNA) <- AnalysisAssay

genes <- c("P2RY12", "P2RY13", "CX3CR1", "SELPLG", "CST3", "OLFML3", "VSIR", "SPI1", "TMEM119",
           "LYVE1", "C1orf56", "CSF1R", "C1QA", "MRC1", "AIF1", "C3", "MERTK", "CCL3", "SPP1",
           "CD9", "CD83", "CCL4", "CXCL10", "OLR1", "LIPA", "LGALS1", "SGK1", "IL1B", "APOC1",
           "SOCS6", "MAFB", "CXCR4", "CTSD", "HLA-DQB1", "HLA-DRA", "LGALS3", "TREM2", "APOE",
           "LPL", "CLEC7A", "CD74", "ITGAX", "ITGAM", "MS4A7", "MS4A6A", "MSR1", "HLA-DRB1", "HLA-DRB4",
           "FOS", "HEXA", "AXL", "SERINC3", "LGMN", "TMEM173", "ABCG2", "HLA-DRB5", "TYROBP", "NFKB2",
           "ITM2C", "MS4A4A", "BAX", "AEN", "ATM", "FAS", "MAPK8", "MDM2", "TP53", "XPC", "GPNMB",
           "PHLDA3", "GADD45A", "DDB2", "RPS27L", "ISG15", "IRF7", "IRF8", "STAT1", "IFITM3",
           "HLA-DQA1", "HLA-DPB1", "CCL2", "CCL8", "CD163", "ANXA3", "PLAC8", "AGR2", "CD44",
           "CCDC173", "PLA2G7", "ZMAT3", "PHPT1", "IFI44L", "CD40", "CHI3L1", "CIITA", "CD36",
           "SELENOP", "RNASE1", "HLA-DPA1", "CD81", "ACTB", "CYBA", "B2M", "IFI27", "MTRNR2L6",
           "SREBF2", "MTRNR2L10", "MX1", "IL3RA", "IFIT1", "IFIT2", "MX2", "IFI6", "PPARG", "LDHA",
           "ACLY", "RXRA", "HMGCS1", "HGMCR1", "HBB", "HBA1", "ALAS2", "HBD", "PLIN2", "IER", "CTSB",
           "BIN1", "CD33", "MKI67", "XRCC5", "IFITM3", "IFITM1", "TFE3", "TFEC", "TFEB", "MITF", "EGR2", "ATF3",
           "PLAT", "CCR7", "PPARD", "ACY3", "IL4I1", "LY6E", "IFIT3", "XAF1", "ISG20")

if (TRUE){
  gene_list <- subset(genes, genes %in% rownames(scRNA[[AnalysisAssay]]$scale.data))
  
  if (!(dir.exists(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay))))){
    dir.create(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay)))
    dir.create(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "UMAP"))
    dir.create(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin"))
  }
  
  if (length(unique(scRNA@meta.data[[metaName]])) > 1){
    if (!(dir.exists(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("GroupBy_", metaName))))){
      dir.create(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("GroupBy_", metaName)))
    }
    if (!(dir.exists(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("SplitBy_", metaName))))){
      dir.create(file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("SplitBy_", metaName)))
    }
  }
  
  for (i in 1:length(gene_list)) {
    a <- NewFeaturePlot(object = scRNA, pt.size = 2,
                        features = gene_list[i],
                        cols = colfunc(4),
                        reduction = "umap.rpca",
                        order = T,
                        label = F)
    
    a <- a + theme(axis.line=element_line(linewidth=1),
                   axis.ticks.length = unit(0, 'cm'),
                   axis.text = element_text(family='Arial', face='bold', color ='black', size=0),
                   text = element_text(family='Arial', face='bold', color ='black', size=24),
                   panel.border = element_blank(),
                   plot.margin = unit(c(0,1,0,0), "cm"),
                   legend.position='none',
                   legend.background = element_rect(color = NA)) +
      xlab("UMAP 1") + ylab("UMAP 2") + labs(title=NULL) + theme_nothing()
    
    png(filename=file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "UMAP", paste0("UMAP_", Exp, "_", ID, "_", gene_list[[i]], ".png")),
         width=5, height=5, units="in", res=300)
    
    print(a)
    
    dev.off()
    
    # Create Violin plots for each gene
    VlnPlot(object = scRNA,
            pt.size = 0,
            cols = Tcols,
            features = gene_list[[i]])
    
    ggsave(filename=file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("Vln_", Exp, "_", ID, "_", gene_list[[i]], ".png")),
           width=10, height=5, units="in", dpi=300)
    
    if (length(unique(scRNA@meta.data[[metaName]])) > 1){
      # Create Violin plots for each gene grouped by metaName
      VlnPlot(object = scRNA,
              pt.size = 0,
              features = gene_list[[i]],
              group.by = metaName)
      
      ggsave(filename=file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("GroupBy_", metaName), paste0("Vln_", Exp, "_", ID, "_", gene_list[[i]], "_GroupBy_", metaName, ".png")),
             width=2.5*length(unique(scRNA@meta.data[[metaName]])),
             height=5,
             units="in",
             dpi=300)
      
      # Create Violin plots for each gene split by metaName
      VlnPlot(object = scRNA,
              pt.size = 0,
              features = gene_list[[i]],
              cols = treatCols,
              split.by = metaName)
      
      ggsave(filename=file.path(outDir, "Gene_Plots", paste0("Assay_", AnalysisAssay), "Violin", paste0("SplitBy_", metaName), paste0("Vln_", Exp, "_", ID, "_", gene_list[[i]], "_SplitBy_", metaName, ".png")),
             width=2*length(unique(Idents(scRNA))),
             height=5,
             units="in",
             dpi=300)
    }
  }
}


clust <- "IFN"

if (TRUE){
  
  vlnOut <- file.path(outDir, 
                      "Gene_Plots", 
                      paste0("Assay_", AnalysisAssay), 
                      "Violin", 
                      paste0("GroupBy_", metaName), paste0(clust, "Only"))
  
  if (!(dir.exists(vlnOut))){
    dir.create(vlnOut, recursive = T)
  }
  
  for (i in 1:length(gene_list)) {
    
    if (length(unique(scRNA@meta.data[[metaName]])) > 1){
      # Create Violin plots for each gene grouped by metaName
      b <- VlnPlot(object = scRNA, 
                   layer = "scale.data",
                   pt.size = 0, 
                   idents = clust,
                   cols = treatCols,
                   features = gene_list[[i]],
                   group.by = metaName)
      
      b <- b + theme_classic() +
        theme(axis.line = element_line(linewidth=1),
              axis.text = element_text(size = 0),
              axis.ticks = element_line(linewidth = 0),
              axis.title = element_text(face = 'bold', size = 24),
              plot.title = element_text(face = 'bold'),
              legend.position = "NA") +
        ylab("Expression") +
        xlab(clust) + ggtitle(NULL)
      
      ggsave(filename=file.path(vlnOut, paste0("Vln_", Exp, "_", ID, "_", gene_list[[i]], "_GroupBy_", metaName, ".png")),
             width=4,
             height=3.5,
             units="in",
             dpi=300)
    }
  }
}

# ============================================================================================================
# ================================= Save the barcodes for each cluster =======================================
# ============================================================================================================

for (j in levels(Idents(scRNA))){
  
  assign(paste0("Cells_", j), data.frame(Cells = WhichCells(scRNA, idents=j)))
  
  write.table(eval(parse(text=paste0("Cells_", j))),
              file.path(outDir, "Barcodes", paste0(Exp, "_", ID, "_Clust_", j, "_Barcodes.tsv")),
              quote=F, sep='\t', row.names = F)
}



