# Detaching all packages
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T,unload = T,force = T))

# Loading new packages
suppressMessages({
  library(Seurat)
  library(sctransform)
  library(dplyr)
  library(cowplot)
  library(ggplot2)
  library(jsonlite)
})
# ===================================== Updated Functions =========================================
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

LoadAlignedData <- function(base.path=NULL, ID=NULL, metaCol=NULL, metaVal=NULL)
{
  
  # Load the CellRanger data
  CR.path <- file.path(base.path, "Cell_Ranger", ID)
  CR.data <- Read10X(paste0(CR.path, "/outs/filtered_feature_bc_matrix/"))
  
  # If samples were aligned to a mixed genome, this will isolate the human gene names
  Hu <- grep("GRCh", rownames(CR.data), value=T)
  # If a combination of species genes were detected, subset to only retain human genes
  if (length(Hu) > 0){
    CR.data <- CR.data[which(rownames(CR.data) %in% Hu), ]
    CR.genes <- as.matrix(rownames(CR.data))
    rownames(CR.data) <- sapply(strsplit(CR.genes, split='_'), '[', 2)
  }
  
  # Create a Seurat object with the human gene data
  obj <- CreateSeuratObject(CR.data,
                            min.cells = 3,
                            min.features = 0,
                            project = ID)
  
  # Add the sample ID to the beginning of the cell barcodes
  obj <- RenameCells(obj, add.cell.id = gsub("_","", ID))
  
  # If values for both metadata column names (metaCol) and metadata values (metaVal)
  # exits and are the same length, metadata will be added to the Seurat object
  if (!(is.null(metaCol)) && !(is.null(metaVal)) && (length(metaCol) == length(metaVal))) {
    for (i in 1:length(metaCol)) {
      obj[[metaCol[[i]]]] <- factor(metaVal[[i]])
    }
  }
  
  return(obj)
  
}

CreateQCDirectories <- function(base.path=NULL, ID=NULL)
{
  # Create the Seurat output directory if it doesn't exists
  dir.create(file.path(base.path, "Seurat"), showWarnings = F)
  # Creating the necessary Contamination and Preliminary analysis output directories
  # Warnings are silenced to quiet the process in case directories already exist
  for(dir in c("Barcodes", "Plots")){
    # Creating the above directories at the con.path and pre.dir locations
    dir.create(file.path(con.path, ID, dir), recursive = T, showWarnings = F)
    dir.create(file.path(pre.dir, ID, dir), recursive = T, showWarnings = F)
  }
  # Create directories in the contamination directory for cluster vs contamination barcodes
  dir.create(file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes"), recursive = T, showWarnings = F)
  dir.create(file.path(con.path, ID, "Barcodes", "Cluster_Barcodes"), recursive = T, showWarnings = F)
}


# =================================================================================================
# Setting a random seed to get consistent results across analysis days
set.seed(030989)
# Setting the color palette for the output files. Currently has colors for up to 20 clusters
cols <- c("black", "red", "darkgreen", "blue", "purple", "orange", "yellow", "lightblue", "darkred", 
          "magenta", "green", "yellow4", "wheat4", "gray48", "plum1", "chocolate1", "bisque1", 
          "cornflowerblue", "burlywood4", "lightcoral")
# ===================================== Import Data from Cell Ranger =========================================
# User inputs the basic  the location and sample ID for the file being processed
storageDir <- "/path/to/main/storage/directory/"
study <- "2024_Chadarevian_NEP"
Exp <- "Evans_2023" # Change to match current analysis

# Current sample ID
ID <- "Met2" # Change to match current sample

# User can specify metadata to add to the dataset
# Multiple columns can be added, but only a single value for each column 
# will be added (i.e., all cells will have the same value for each column)
metaCol <- c("Condition") # Column name for the metadata to be added
metaVal <- c("Metastasis") # Value to be added to the dataset

# Generate a pointer to the base.path for the experimental input and output
base.path <- file.path(storageDir, study, Exp)

# Setting the output file directories
con.path <- file.path(base.path, "Seurat", "Contamination")
pre.dir <- file.path(base.path, "Seurat", "Preliminary")

# Loading the CellRanger dataset
scRNA <- LoadAlignedData(base.path=base.path, ID = ID, metaCol = metaCol, metaVal = metaVal)

# Creating the necessary output directories
CreateQCDirectories(base.path, ID)

# ================================================================================================================
# =========================================== Begin QC Filtering =================================================
# ================================================================================================================

if(TRUE){
  # Adding metadata columns for the percent of mitochondrial and ribosomal reads
  scRNA[["percent.mito"]] <- PercentageFeatureSet(object = scRNA, pattern = "^MT-")
  scRNA[["percent.ribo"]] <- PercentageFeatureSet(object = scRNA, pattern = "^RP[SL]")
  
  # Visualizing the basic QC metrics for the dataset
  cowplot::plot_grid(FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.mito"),
                     FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo"),
                     FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), 
                     nrow = 3, ncol=1)
  
}

# ===============================================================================================================
# ======================================= Scaling and Normalization =============================================
# ===============================================================================================================

if(TRUE){
  # Normalizing the data using default settings
  scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)
  
  # Scaling and centering the data using the full transcriptome
  scRNA <- ScaleData(scRNA, features = rownames(scRNA), block.size = 10000)
  
  # Performing PCA analysis on the full transcriptome
  scRNA <- RunPCA(object = scRNA, features = rownames(scRNA), do.print = TRUE, ndims.print = 1:5, 
                  nfeatures.print = 10)
  
  # Printing an elbow plot for PC selection
  ElbowPlot(object = scRNA, ndims = 50)
  # Save the elbow plot
  ggsave(filename= file.path(con.path, ID, paste0("ElbowPlot_", Exp, "_", ID, ".png")),
         width=8, height=5, units="in", dpi=300, bg="white")
  
  # Printing an elbow plot for PC selection
  ElbowPlot(object = scRNA, ndims = 50)
  
}

# ===============================================================================================================    
# ========================================== Clustering the cells ===============================================
# ===============================================================================================================

# Select the number of PCs to use for the nearest neighbor analysis and the
# resolution for the clustering. For QC clustering, 0.6-0.9 is typically sufficient
# to distinguish the relevant clusters.
PCs <- c(1:20) # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx
Res <- 0.8 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx

# Perform clustering of data
scRNA <- FindNeighbors(scRNA, dims = PCs, force.recalc = T)
scRNA <- FindClusters(object = scRNA, resolution = Res)

# ===================================================================================================================   
# =========================================== Generating UMAP plots =================================================
# ===================================================================================================================

# Running UMAP with the previously calculated parameters
scRNA <- RunUMAP(scRNA, dims = PCs)
DimPlot(scRNA, reduction = "umap", label = T)

# ============================================================================================================
# ======================================= Saving QC Analysis Data ============================================
# ============================================================================================================

# =========================================== Save UMAP plot ================================================= 

DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap') + theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.line=element_line(linewidth=1), 
        axis.ticks = element_line(color='black', linewidth=3), 
        axis.ticks.length = unit(0, 'cm'), 
        axis.text = element_blank(),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename= file.path(con.path, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Clusters.png")),
       width=5, height=5, units="in", dpi=300, bg='white')

# ==================================== Add Cluster Signature Scores ==========================================   

# Load dividing gene list
divGenes <- read.table("/path/to/QC_DividingGenes.tsv", header = T)

# Add dividing gene signature
scRNA <- AddModuleScore(scRNA, 
                        features = as.data.frame(divGenes), 
                        name = "Dividing",
                        assay = "RNA",
                        layer = "data")

# Create a list of signatures names to call in FeaturePlot
featureList <- c("nFeature_RNA", "nCount_RNA", "percent.ribo", 
                 "percent.mito", "Dividing1")
# Creating color palette for heatmap
colfunc <- colorRampPalette(c("midnightblue", "lightskyblue1", "orange", "red"))
# Save a layout with UMAPs for all signatures
NewFeaturePlot(scRNA, pt.size = 0.5,
            features = featureList, cols = colfunc(10), 
            reduction = "umap",
            order = T)

ggsave(filename= file.path(con.path, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Signatures.png")), 
       width=14, height=8, units="in", dpi=100, bg='white')

# ================================= Save the barcodes for each cluster =======================================

for (i in unique(eval(parse(text=paste0("scRNA$RNA_snn_res.", Res))))){
  
  assign(paste0("Cells", i), data.frame(Cells = WhichCells(scRNA, idents=i)))
  
  write.table(eval(parse(text=paste0("Cells", i))), 
              file.path(con.path, ID, "Barcodes", "Cluster_Barcodes", paste0(Exp, "_", ID, "_Clust", i, "_Barcodes.tsv")), 
              quote=F, sep='\t', row.names = F)
}

# Store all barcodes in the dataset for use later
Allcodes <- colnames(scRNA)

# ======================================== Save QC Scatter plots =============================================
cowplot::plot_grid(FeatureScatter(object = scRNA, 
                                  feature1 = "nCount_RNA", 
                                  feature2 = "percent.mito", 
                                  cols = cols) + NoLegend(),
                   FeatureScatter(object = scRNA, 
                                  feature1 = "nCount_RNA", 
                                  feature2 = "percent.ribo",
                                  cols = cols) + NoLegend(),
                   FeatureScatter(object = scRNA, 
                                  feature1 = "nCount_RNA", 
                                  feature2 = "nFeature_RNA",
                                  cols = cols) + NoLegend(), nrow = 3, ncol=1)

ggsave(filename = file.path(con.path, ID, "Plots", paste0("Scatter_", Exp, "_", ID, ".png")), 
       width=12, height=9, units="in", dpi=300, bg='white')

# ============================================================================================================
# ==================================== Removing Contaminating Cells ==========================================
# ================================ Identify and save contaminating barcodes ==================================
# ============================================================================================================

# Dividing cells cluster numbers
div <- c(7, 9) # Include cluster numbers that are high in the dividing signature
# Damaged cells cluster numbers
dam <- c(11, 12, 13) # Include cluster numbers that are low in genes/UMIs or high in percent.mito/percent.ribo
# Doublet cells cluster numbers
doub <- NULL # Include cluster numbers that are abnormally high in genes/UMIs (doublets)

# Subset the Seurat object to remove the cells that didn't pass QC
test <- subset(scRNA, cells = WhichCells(scRNA, idents = c(div, dam, doub), invert = T))

# Visualize a preview of the sample QC data after applying the above cutoffs
cowplot::plot_grid(FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = cols),
                   FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "percent.ribo", cols = cols),
                   FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = cols), 
                   nrow =3, ncol=1)

if(TRUE){
  # Write dividing cell barcodes to a file
  if (is.null(div)){
    dividing <- data.frame(Cells="NULL")
    write.table(dividing, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", 
                                    "Dividing_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  } else{
    dividing <- data.frame(Cells=WhichCells(scRNA, idents = div))
    write.table(dividing, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", 
                                    "Dividing_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  }
  # Write damaged cell barcodes to a file
  if (is.null(dam)){
    damaged <- data.frame(Cells=NULL)
    write.table(damaged, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes",
                                   "GenePoor_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  } else {
    damaged <- data.frame(Cells=WhichCells(scRNA, idents = dam))
    write.table(damaged, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", 
                                   "GenePoor_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  }
  # Write doublet cell barcodes to a file
  if (is.null(doub)){
    doublet <- data.frame(Cells="NULL")
    write.table(doublet, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", 
                                   "Doublet_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  }else{
    doublet <- data.frame(Cells=WhichCells(scRNA, idents = doub))
    write.table(doublet, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", 
                                   "Doublet_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  }
}

# ========================= Check to see if the cutoffs look right before continuing ===========================  
# After finalizing the cutoffs, overwrite previous Seurat object
scRNA <- test
# Remove the temporary object
rm(test)
# Subset Allcodes to remove the above barcodes
Allcodes <- subset(Allcodes, !(Allcodes %in% c(damaged$Cells, dividing$Cells, doublet$Cells)))

# ========================== Remove Contaminating barcodes and recalc cell totals =============================    
# Recalculate the dataset details
details <- CollectObjectDetails(scRNA, output=T)
# Writing the post-contamination details to the analysis dataframe
postCon <- data.frame(Analysis_Date = day, 
                      Time = emit, 
                      Seurat_Version = versions[["Seurat"]],
                      CellRanger_Version = versions[["CellRanger"]], 
                      R_Version = versions[["R"]],
                      Analysis_Stage = "Post_Contamination", 
                      Sample_ID = ID, 
                      Total_Cells = details[["Cell_Count"]], 
                      Average_UMI = details[["Mean_UMI"]], 
                      Median_UMI = details[["Median_UMI"]], 
                      Average_Gene = details[["Mean_Gene"]], 
                      Median_Gene = details[["Median_Gene"]])

# Merging the postCon datafrome with the preCon details
SampleInfo <- merge(preCon, postCon, all = T, sort=F)

# ============================ Additional filtering of remaining cells =======================================  
# Visualize the current sample QC data
cowplot::plot_grid(FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = cols),
                   FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", cols = cols),
                   FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = cols), 
                   nrow = 3, ncol=1)

# Set high and low gene cutoffs
geneLow <- 1000
geneHigh <- (2*details[["Median_Gene"]])
# Set high and low UMI cutoffs
UMILow <- 1000
UMIHigh <- (2*details[["Median_UMI"]])
# Set mitochondrial percentage cutoffs
MitoLow <- -Inf
MitoHigh <- 20 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx
# Set ribosomal percentage cutoffs
RiboLow <- -Inf
RiboHigh <- 30 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx

# Subset a test dataset to determine if the cutoffs are sufficient
test <- subset(scRNA, subset = nFeature_RNA > geneLow & nFeature_RNA < geneHigh)
test <- subset(test, subset = nCount_RNA > UMILow & nCount_RNA < UMIHigh)
test <- subset(test, subset = percent.mito > MitoLow & percent.mito < MitoHigh)
test <- subset(test, subset = percent.ribo > RiboLow & percent.ribo < RiboHigh)

# Visualize a preview of the sample QC data after applying the above cutoffs
cowplot::plot_grid(FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "percent.mito", cols = cols),
                   FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "percent.ribo", cols = cols),
                   FeatureScatter(object = test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = cols), 
                   nrow = 3, ncol=1)

# ========================= Check to see if the cutoffs look right before continuing ===========================  
# After finalizing the cutoffs, overwrite previous Seurat object
scRNA <- test
# Remove the temporary object
rm(test)
# Identify and write the barcodes that were removed
QCcodes <- data.frame(Cells = subset(Allcodes, !(Allcodes %in% colnames(scRNA))))
write.table(QCcodes, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", "Cutoff_Barcodes.tsv"), 
            sep='\t', quote=F, row.names = F)


# ====================================== Save the Seurat Object ==============================================  
# Save the Seurat object with the cells that passed QC
saveRDS(scRNA, file = file.path(pre.dir, ID, paste0(Exp, "_", ID, "_Seurat_Object.rds")))

# ================================= Save the barcodes for each cluster =======================================

for (i in unique(eval(parse(text=paste0("scRNA$RNA_snn_res.", Res))))){
  
  assign(paste0("Cells", i), data.frame(Cells = WhichCells(scRNA, idents=i)))
  
  write.table(eval(parse(text=paste0("Cells", i))), 
              file.path(pre.dir, ID, "Barcodes", paste0(Exp, "_", ID, "_Clust", i, "_Barcodes.tsv")), 
              quote=F, sep='\t', row.names = F)
}

