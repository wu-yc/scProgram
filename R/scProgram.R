#' scProgram
#'
#' inMatrix1
#' @param countexp
#' @keywords scProgram
#' @examples
#' dist.JSD_2m_fast()
#' GetFeatures()
#' HeatFeatures()
#' GetProgram()
#' @export dist.JSD_2m_fast
#' @export GetFeatures
#' @export HeatFeatures
#' @export GetProgram
#'
#'
dist.JSD_2m_fast  <- function(inMatrix1, inMatrix2) {
  matrixColSize1 <- ncol(inMatrix1)
  matrixColSize2 <- ncol(inMatrix2)
  resultsMatrix <- matrix(0, matrixColSize1, matrixColSize2)

  pb = dplyr::progress_estimated(matrixColSize1, min_time = 3)

  resultsMatrix <- apply(inMatrix1, 2, function(x) {
    pb$tick()$print()
    apply(inMatrix2, 2, function(y) {
      philentropy::jensen_shannon(x, y, F, unit = "log2")
    })
  })
  resultsMatrix = (t(resultsMatrix))

  return(resultsMatrix)
}



GetFeatures <- function(obj,
                        group.by = "cluster",
                        assay='RNA',
                        slot='data',
                        mode = "fast",
                        pct_exp = 0.25,
                        downsample = T,
                        downsample_size = 100,
                        genenumber = 50){


  library(philentropy)
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyverse)
  # library(Matrix)

  #mode: fast, standard
  #obj: Seurat object
  #group.by: Seurat object中想差异分析的列
  #features: 想纳入分析的基因，默认是所有


  if (mode == "standard"){
    cat("Standard mode\n")
    Idents(obj) = obj@meta.data[,group.by]
    markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = pct_exp, logfc.threshold = 0.25, verbose = F)
    top_markers <- markers %>%
      group_by(cluster) %>%
      top_n(genenumber, avg_log2FC) %>%
      mutate(row_num = row_number())
    marker_df <- top_markers %>%
      pivot_wider(names_from = cluster, values_from = gene, id_cols = row_num)
    marker_df$row_num <- NULL
    return(marker_df)
  }


  if (mode == "fast"){
    cat("Fast mode\n")


    t1 = Sys.time()

    # Pre-allocate matrices
    # genemat <- matrix(0, nrow = nrow(obj), ncol = ncol(obj))

    # Subset directly instead of assigning to new object
    metadata <- data.table(obj@meta.data)


    # Extract data directly to pre-allocated matrix
    genemat <- GetAssayData(object = obj, slot = slot, assay = assay)
    # genemat <- genemat[rowMeans(genemat) > 0, ]


    # Filter genes
    # low_expressed_genes <- rowSums(genemat) / colSums(genemat) > pct_exp
    # genemat <- genemat[low_expressed_genes, ]

    # percent_expressed <- apply(genemat, 1, function(x) sum(x > 0) / ncol(genemat) * 100)
    # filtered_matrix <- genemat[percent_expressed >= 10,]

    row_sums <- Matrix::rowSums(genemat)
    percent_expressed <- row_sums / ncol(genemat)
    genemat <- genemat[percent_expressed > pct_exp,]


    # Downsample
    if(downsample){

      cell_cluster_df <- data.table(Cell = colnames(genemat),
                                    Cluster = metadata[,group.by, with = FALSE])

      keep_cells <- unlist(lapply(split(cell_cluster_df, cell_cluster_df$Cluster),
                                  function(x) sample(x$Cell, size = min(nrow(x), downsample_size))))

      genemat2 <- genemat[, keep_cells]
      metadata <- obj@meta.data[keep_cells, ]
    }

    # Log transform
    genemat2 <- log2(genemat2/Matrix::rowSums(genemat2) + 1)


    all.ident = unique(as.character(metadata[,group.by]))
    all.ident <- all.ident[order(all.ident, decreasing = F)]
    all.ident = all.ident[!is.na(all.ident)]

    pseudo.comparemat = data.table(name = row.names(metadata),
                                   ident = metadata[,group.by])

    pseudo.comparemat_simple = matrix(NA, nrow = nrow(pseudo.comparemat), ncol = length(all.ident))
    colnames(pseudo.comparemat_simple) = all.ident
    for(j in seq_along(all.ident)) {
      pseudo.comparemat_simple[,j] = as.numeric(pseudo.comparemat$ident == all.ident[j])
    }
    pseudo.comparemat_simple2 = log2(pseudo.comparemat_simple/rowSums(pseudo.comparemat_simple) + 1)


    sum.distance.JSD = dist.JSD_2m_fast(t(as.matrix(genemat2)), (pseudo.comparemat_simple2))

    row.names(sum.distance.JSD) = row.names(genemat2)
    colnames(sum.distance.JSD) = row.names(pseudo.comparemat_simple2)


    sum.distance.JSD_ordered = data.frame(matrix(NA, genenumber, ncol(pseudo.comparemat_simple)))

    #0.08 sec
    for (j in 1:ncol(sum.distance.JSD_ordered)){
      tmpname <- row.names(sum.distance.JSD[order(sum.distance.JSD[,j], decreasing = F),])[1:genenumber]
      sum.distance.JSD_ordered[,j] = tmpname
    }
    colnames(sum.distance.JSD_ordered) = all.ident
    cat(paste("\nElapsed time: ", as.numeric(difftime(Sys.time(), t1, units = "secs")), " seconds\n", sep=""))

    return(sum.distance.JSD_ordered)
  }
}



HeatFeatures <- function(obj,
                         features  = FeatureMatrix,
                         group.by = "cluster",
                         show_rownames = F, show_colnames = T,
                         cols = c("white","white", "white", "#52A85F"),
                         assay='RNA',
                         slot='data'){

  library(Seurat)
  library(pheatmap)
  library(RColorBrewer)

  avgexp = AverageExpression(obj, group.by = group.by, return.seurat = T)
  heatmat <- GetAssayData(object = avgexp, slot = slot, assay = assay)



  # cc <- colorRampPalette(c("white", "grey97", "#F47E5D", "#463873")) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"
  cc <- colorRampPalette(cols) #
  pheatmap(heatmat[unique(as.character(as.matrix(FeatureMatrix))), colnames(FeatureMatrix)],
           show_rownames = show_rownames, show_colnames = show_colnames,
           color = cc(100),
           cluster_rows = F, cluster_cols = F, scale = "row")

}



GetProgram <- function(features  = FeatureMatrix,
                       geneset = "KEGG", #"HALLMARK"
                       pvalue_cutoff = 0.05,
                       cols = c("#F47E5D", "#CA3D74", "#7F2880", "#463873"),
                       plot_term_number =3
){

  library(clusterProfiler)
  library(ggplot2)

  # gmtfile <- "/www/data/signature/h.all.v7.2.symbols.c2.cp.kegg.v7.2.symbols.gmt" #

  gmtfile <- system.file("data", "h.all.v7.2.symbols.c2.cp.kegg.v7.2.symbols.gmt", package = "scProgram")

  c5 <- read.gmt(gmtfile)

  if (geneset == "KEGG") c5 = subset(c5, term %like% "KEGG")
  if (geneset == "HALLMARK") c5 = subset(c5, term %like% "HALLMARK")

  ego_comb<-character(0)
  cluster.names<-colnames(FeatureMatrix)
  for (i in 1:(length(cluster.names))) {
    # print(cluster.names[i])
    gene_list0 <- FeatureMatrix[,cluster.names[i]]

    ego0 <- enricher(gene_list0, TERM2GENE=c5)



    ego0_result<-ego0@result
    ego0_result$GeneRatio.num<-ego0_result$Count/length(gene_list0)
    ego0_result <- ego0_result[order(ego0_result$GeneRatio.num, decreasing = T),]



    ego0_result<-subset(ego0_result, ego0_result$pvalue < pvalue_cutoff)
    ego0_result<-cbind(ego0_result, cluster.names[i])
    colnames(ego0_result)[ncol(ego0_result)]<-"group"
    ego_comb<-rbind(ego_comb, ego0_result)
  }


  ego_comb$Description <- factor(ego_comb$Description, levels=rev(unique(ego_comb$Description)))
  ego_comb[1:3,]


  ego_comb_sub <- ego_comb %>% group_by(group) %>% top_n(-plot_term_number, pvalue)

  pal <- colorRampPalette(cols)(100) #"white", "#F47E5D", "#CA3D74", "#7F2880", "#463873"

  ego_comb_sub$logp<- -log10(ego_comb_sub$pvalue)
  ego_comb_sub$Description <- factor(ego_comb_sub$Description, levels=rev(unique(ego_comb_sub$Description)))
  ego_comb_sub$group <- factor(ego_comb_sub$group, levels=rev(unique(ego_comb_sub$group)))

  ggplot(ego_comb_sub, aes(x = group, y = Description)) +
    geom_point(aes(size = Count, color = logp)) +
    theme_bw(base_size = 14) +
    #scale_colour_gradient(limits=c(0, 0.10), low="red") +
    scale_color_gradientn(colours = pal) +
    ylab(NULL) +
    theme(axis.text.x=element_text(angle=45,hjust=1))+ #aspect.ratio=1,
    ggtitle("")

  # return(ego_comb)
}




