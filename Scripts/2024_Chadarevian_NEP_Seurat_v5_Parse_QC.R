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
  library(Matrix)
  library(jsonlite)
  library(scDblFinder)
})

# Setting a random seed to get consistent results across analysis days
set.seed(030989)
# Setting the color palette for the output files. Currently has colors for up to 20 clusters
Tcols <- c("black", "red", "darkgreen", "blue", "purple", "orange", "yellow", "lightblue", "darkred",
          "magenta", "green", "yellow4", "wheat4", "gray48", "plum1", "chocolate1", "bisque1",
          "cornflowerblue", "burlywood4", "lightcoral")
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

# ===================================== Import Data from Cell Ranger =========================================
# User inputs the basic  the location and sample ID for the file being processed
storageDir <- "/path/to/main/storage/directory"
study <- "2024_Chadarevian_NEP"
Exp <- "Cuprizone"
ID <- "Sublibrary1"

# User can specify metadata to add to the dataset
# Multiple columns can be added, but only a single value for each column
# will be added (i.e., all cells will have the same value for each column)
metaCol <- c("Sublibrary") # Column name for the metadata to be added
metaVal <- c("Sublibrary1") # Value to be added to the dataset

# Generate a pointer to the base.path for the experimental input and output
# Ensure that these align with how the data is organized on your system
base.path <- file.path(storageDir, study, Exp)
split.path <- file.path(base.path, paste0("SplitPipe_Output/", ID, "/all-sample/DGE_filtered"))

# Setting the output file directories
con.path <- file.path(base.path, "Seurat", "Contamination")
pre.dir <- file.path(base.path, "Seurat", "Preliminary")

# Creating the necessary output directories
CreateQCDirectories(base.path, ID)

# Loading the Parse combined dataset matrix values
mat <- ReadParseBio(data.dir = split.path)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(split.path, "/cell_metadata.csv"), row.names = 1)

# Add treatment and genotype information to the metadata
sampInfo <- read.table("/path/to/Cuprizone_SampleGrouping.tsv", header = T)
cellCodes <- rownames(cell_meta)
cell_meta <- merge(cell_meta, sampInfo, by = "sample")
rownames(cell_meta) <- cellCodes

# Clean up the space
rm(cellCodes, sampInfo)
gc()

# Create object
scRNA <- CreateSeuratObject(mat,
                            min.features = 100,
                            min.cells = 100,
                            names.field = 0,
                            meta.data = cell_meta)
# ===========================================================
# Update the object metadata
# ===========================================================

# Add current sublibrary column
scRNA@meta.data$Sublibrary <- ID

# Reset the object identities to something generic
scRNA@meta.data$orig.ident <- factor(paste0(ID, "_", scRNA$sample))
Idents(scRNA) <- scRNA$Sublibrary

# Collect and print the current cell, UMI, and gene details
details <- CollectObjectDetails(scRNA, output=T)

# Order the object metadata levels
scRNA$Treatment <- factor(scRNA$Treatment, levels = c("Normal_Chow",
                                                      "Cuprizone"))
# Remove matrix file to save space
rm(mat)

# Clear unused memory
gc()

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
                     FeatureScatter(object = scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"), nrow = 1, ncol=3)

}

# ===============================================================================================================
# ======================================= Scaling and Normalization =============================================
# ===============================================================================================================

if(TRUE){
  # Normalizing the data using default settings
  scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)

  # Identifying variable features to scale on
  scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(scRNA), 10)
  plot1 <- VariableFeaturePlot(scRNA)
  LabelPoints(plot = plot1, points=top10, repel=T)

  # Scaling and centering the data using the full transcriptome
  scRNA <- ScaleData(scRNA, 
                     features = VariableFeatures(scRNA),
                     block.size = 1000)

  # Performing PCA analysis on the full transcriptome
  scRNA <- RunPCA(object = scRNA, features = VariableFeatures(scRNA), do.print = TRUE, ndims.print = 1:5,
                  nfeatures.print = 10, npcs = 100)

  # Printing an elbow plot for PC selection
  ElbowPlot(object = scRNA, ndims = 100)
  # Save the elbow plot
  ggsave(filename= file.path(con.path, ID, paste0("ElbowPlot_", Exp, "_", ID, ".png")),
         width=8, height=5, units="in", dpi=300, bg="white")

  # Saving the Seurat object in case ScaleData crashes due to RAM limit
  saveRDS(scRNA, file = file.path(con.path, ID, paste0(Exp, "_", ID, "_Seurat_Object.rds")))

}

# Printing an elbow plot for PC selection
ElbowPlot(object = scRNA, ndims = 100)

# ===============================================================================================================
# ========================================== Clustering the cells ===============================================
# ===============================================================================================================

# Select the number of PCs to use for the nearest neighbor analysis and the
# resolution for the clustering. For QC clustering, 0.6-0.9 is typically sufficient
# to distinguish the relevant clusters.
PCs <- c(1:25) # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx
Res <- 0.8 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx

# Perform clustering of data
scRNA <- FindNeighbors(scRNA, dims = PCs, force.recalc = T)
scRNA <- FindClusters(object = scRNA, resolution = Res)

# ===================================================================================================================
# =========================================== Generating UMAP plots =================================================
# ===================================================================================================================

# Running UMAP with the previously calculated parameters
scRNA <- RunUMAP(scRNA, dims = PCs)
DimPlot(scRNA, reduction = "umap", label = T, raster = F)

# ================================= Save the barcodes for each cluster =======================================

for (i in unique(eval(parse(text=paste0("scRNA$RNA_snn_res.", Res))))){

  assign(paste0("Cells", i), data.frame(Cells = WhichCells(scRNA, idents=i)))

  write.table(eval(parse(text=paste0("Cells", i))),
              file.path(con.path, ID, "Barcodes", "Cluster_Barcodes", paste0(Exp, "_", ID, "_Clust", i, "_Barcodes.tsv")),
              quote=F, sep='\t', row.names = F)
}

# Store all barcodes in the dataset for use later
Allcodes <- colnames(scRNA)

# ============================================================================================================
# ========================================= Removing Doublets ================================================
# ================================ Identify and save contaminating barcodes ==================================
# ============================================================================================================

# Identify and remove doublets
sce <- as.SingleCellExperiment(scRNA)
sce <- scDblFinder(sce, clusters = 'seurat_clusters')

scRNA@meta.data$scDblFinder.score <- sce$scDblFinder.score
scRNA@meta.data$scDblFinder.class <- sce$scDblFinder.class
scRNA@meta.data$scDblFinder.weighted <- sce$scDblFinder.weighted

DimPlot(scRNA, group.by = "scDblFinder.class")

ggsave(filename= file.path(con.path, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Doublets.png")),
       width=7, height=5, units="in", dpi=300, bg='white')

doublet <- data.frame(Cells = rownames(scRNA@meta.data[scRNA$scDblFinder.class == "doublet",]))

# Write doublet cell barcodes to a file
write.table(doublet, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes",
                               "Doublet_Barcodes.tsv"), sep='\t', quote=F, row.names = F)

if (TRUE) {
# Make doublet object
doubs <- subset(scRNA, cells = doublet$Cells)

doubFile <- file(file.path(con.path, ID, "DoubletInformation.txt"))
writeLines(c("Original Count by Group",
             paste0('\t', names(table(scRNA@meta.data[[metaName]])), ": ", table(scRNA@meta.data[[metaName]])),
             "Doublet Count by Group",
             paste0('\t', names(table(doubs@meta.data[[metaName]])), ": ", table(doubs@meta.data[[metaName]])),
             "Percent Doublets by Group",
             paste0('\t', names(table(doubs@meta.data[[metaName]])), ": ", (table(doubs@meta.data[[metaName]])/table(scRNA@meta.data[[metaName]]))*100),
             "Percent Doublets by Sample",
             paste0('\t', names(table(doubs$sample)), ": ", (table(doubs$sample)/table(scRNA$sample))*100)), doubFile)
close(doubFile)
}

# Remove doublets
scRNA <- subset(scRNA, cells = doublet$Cells, invert = T)


# ============================================================================================================
# ======================================= Saving QC Analysis Data ============================================
# ============================================================================================================
# =========================================== Save UMAP plot =================================================

DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap', raster = F) + theme_classic() +
  scale_color_manual(values = Tcols) +
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


DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap', raster = F, split.by = "sample") +
  theme_classic() +
  scale_color_manual(values = Tcols) +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks = element_line(color='black', linewidth=3),
        axis.ticks.length = unit(0, 'cm'),
        axis.text = element_blank(),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename= file.path(con.path, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Clusters_bySample.png")),
       width=3*(length(unique(scRNA$sample))), height=5, units="in", dpi=300, bg='white')



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


# ======================================== Save QC Scatter plots =============================================
cowplot::plot_grid(FeatureScatter(object = scRNA,
                                  feature1 = "nCount_RNA",
                                  feature2 = "percent.mito",
                                  cols = Tcols,
                                  raster = F) + NoLegend(),
                   FeatureScatter(object = scRNA,
                                  feature1 = "nCount_RNA",
                                  feature2 = "percent.ribo",
                                  cols = Tcols,
                                  raster = F) + NoLegend(),
                   FeatureScatter(object = scRNA,
                                  feature1 = "nCount_RNA",
                                  feature2 = "nFeature_RNA",
                                  cols = Tcols,
                                  raster = F) + NoLegend(), nrow = 1, ncol=3)

ggsave(filename = file.path(con.path, ID, "Plots", paste0("Scatter_", Exp, "_", ID, ".png")),
       width=12, height=9, units="in", dpi=300, bg='white')

# ============================================================================================================
# =============================== Removing Remaining Contaminating Cells =====================================
# ================================ Identify and save contaminating barcodes ==================================
# ============================================================================================================
# Dividing cells cluster numbers
div <- c(14, 22, 15, 16) # Include cluster numbers that are high in the dividing signature

# Damaged cells cluster numbers
dam <- NULL # Include cluster numbers that are low in genes/UMIs or high in percent.mito/percent.ribo

# Subset the Seurat object to remove the cells that didn't pass QC
if (length(c(div, dam)) > 0){
  test <- subset(scRNA, cells = WhichCells(scRNA, idents = c(div, dam), invert = T))

  # Visualize a preview of the sample QC data after applying the above cutoffs
  cowplot::plot_grid(FeatureScatter(object = test,
                                    feature1 = "nCount_RNA", feature2 = "percent.mito",
                                    raster = F, cols = Tcols),
                     FeatureScatter(object = test,
                                    feature1 = "nCount_RNA", feature2 = "percent.ribo",
                                    raster = F, cols = Tcols),
                     FeatureScatter(object = test,
                                    feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                                    raster = F, cols = Tcols),
                     nrow=1, ncol=3)

  NewFeaturePlot(test, pt.size = 0.5, raster = F,
                 features = featureList,
                 cols = colfunc(20),
                 reduction = "umap",
                 order = T, label=T)
} else {

  print("No clusters are being removed...")
  test <- scRNA

  }

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
    damaged <- data.frame(Cells="NULL")
    write.table(damaged, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes",
                                   "GenePoor_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  } else {
    damaged <- data.frame(Cells=WhichCells(scRNA, idents = dam))
    write.table(damaged, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes",
                                   "GenePoor_Barcodes.tsv"), sep='\t', quote=F, row.names = F)
  }


  # ========================= Check to see if the cutoffs look right before continuing ===========================
  # After finalizing the cutoffs, overwrite previous Seurat object
  scRNA <- test
  # Remove the temporary object
  rm(test)
  # Subset Allcodes to remove the above barcodes
  Allcodes <- subset(Allcodes, !(Allcodes %in% c(damaged$Cells, dividing$Cells, doublet$Cells)))
}

# ============================ Additional filtering of remaining cells =======================================
# Visualize the current sample QC data
cowplot::plot_grid(FeatureScatter(object = scRNA,
                                  feature1 = "nCount_RNA", feature2 = "percent.mito", raster=F, ),
                   FeatureScatter(object = scRNA,
                                  feature1 = "nCount_RNA", feature2 = "percent.ribo", raster=F),
                   FeatureScatter(object = scRNA,
                                  feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=F),
                   nrow = 1, ncol=3)

# Set high and low gene cutoffs
geneLow <- 500 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx
geneHigh <- (2*median(scRNA[["nFeature_RNA"]][[1]]))

# Set high and low UMI cutoffs
UMILow <- 500 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx
UMIHigh <- (2*median(scRNA[["nCount_RNA"]][[1]]))

# Set mitochondrial percentage cutoffs
MitoLow <- -Inf
MitoHigh <- 10 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx

# Set ribosomal percentage cutoffs
RiboLow <- -Inf
RiboHigh <- 1.5 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx

# Subset a test dataset to determine if the cutoffs are sufficient
test <- subset(scRNA, subset = nFeature_RNA > geneLow & nFeature_RNA < geneHigh)
test <- subset(test, subset = nCount_RNA > UMILow & nCount_RNA < UMIHigh)
test <- subset(test, subset = percent.mito > MitoLow & percent.mito < MitoHigh)
test <- subset(test, subset = percent.ribo > RiboLow & percent.ribo < RiboHigh)

# Visualize a preview of the sample QC data after applying the above cutoffs
cowplot::plot_grid(FeatureScatter(object = test,
                                  feature1 = "nCount_RNA", feature2 = "percent.mito", raster=F),
                   FeatureScatter(object = test,
                                  feature1 = "nCount_RNA", feature2 = "percent.ribo", raster=F),
                   FeatureScatter(object = test,
                                  feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=F),
                   nrow = 1, ncol=3)

# ========================= Check to see if the cutoffs look right before continuing ===========================
# After finalizing the cutoffs, overwrite previous Seurat object
scRNA <- test
# Remove the temporary object
rm(test)
# Identify and write the barcodes that were removed
QCcodes <- data.frame(Cells = subset(Allcodes, !(Allcodes %in% colnames(scRNA))))
write.table(QCcodes, file.path(con.path, ID, "Barcodes", "Contaminating_Barcodes", "Cutoff_Barcodes.tsv"),
            sep='\t', quote=F, row.names = F)

# ===============================================================================================================
# ======================================= Scaling and Normalization =============================================
# ===============================================================================================================

if(TRUE){
  # Normalizing the data using default settings
  scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)

  # Scaling and centering the data using the full transcriptome
  scRNA <- ScaleData(scRNA, features = rownames(scRNA), block.size = 1000)

  # Performing PCA analysis on the full transcriptome
  scRNA <- RunPCA(object = scRNA, features = rownames(scRNA), do.print = TRUE, ndims.print = 1:5,
                  nfeatures.print = 10)

  # Printing an elbow plot for PC selection
  ElbowPlot(object = scRNA, ndims = 50)
  # Save the elbow plot
  ggsave(filename= file.path(pre.dir, ID, paste0("ElbowPlot_", Exp, "_", ID, ".png")),
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
PCs <- c(1:19) # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx
Res <- 1.0 # Change according to Chadarevian_2024_NEP_SampleQCMetadata.xlsx

# Perform clustering of data
scRNA <- FindNeighbors(scRNA, dims = PCs, force.recalc = T)
scRNA <- FindClusters(object = scRNA, resolution = Res)

# ===================================================================================================================
# =========================================== Generating UMAP plots =================================================
# ===================================================================================================================

# Running UMAP with the previously calculated parameters
scRNA <- RunUMAP(scRNA, dims = PCs)
DimPlot(scRNA, reduction = "umap", label = T)


# ====================================== Save the Seurat Object ==============================================
# Save the Seurat object with the cells that passed QC
saveRDS(scRNA, file = file.path(pre.dir, ID, paste0(Exp, "_", ID, "_Seurat_Object.rds")))

# ============================================================================================================
# ==================================== Saving PostFilt Analysis Data =========================================
# ============================================================================================================

# =========================================== Save UMAP plot =================================================

DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap') + theme_classic() +
  scale_color_manual(values = Tcols) +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks = element_line(color='black', linewidth=3),
        axis.ticks.length = unit(0, 'cm'),
        axis.text = element_blank(),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename= file.path(pre.dir, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Clusters.png")),
       width=5, height=5, units="in", dpi=300, bg='white')


DimPlot(object = scRNA, pt.size = 1, label = T, reduction = 'umap', split.by = "sample") +
  theme_classic() +
  scale_color_manual(values = Tcols) +
  theme(axis.line=element_line(linewidth=1),
        axis.ticks = element_line(color='black', linewidth=3),
        axis.ticks.length = unit(0, 'cm'),
        axis.text = element_blank(),
        text = element_text(face='bold', color ='black', size=24),
        panel.border = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        legend.position='none') +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename= file.path(pre.dir, ID, "Plots", paste0("UMAP_", Exp, "_", ID, "_Clusters_bySample.png")),
       width=3*(length(unique(scRNA$sample))), height=5, units="in", dpi=300, bg='white')

# ================================= Save the barcodes for each cluster =======================================

for (i in unique(eval(parse(text=paste0("scRNA$RNA_snn_res.", Res))))){

  assign(paste0("Cells", i), data.frame(Cells = WhichCells(scRNA, idents=i)))

  write.table(eval(parse(text=paste0("Cells", i))),
              file.path(pre.dir, ID, "Barcodes", paste0(Exp, "_", ID, "_Clust", i, "_Barcodes.tsv")),
              quote=F, sep='\t', row.names = F)
}

