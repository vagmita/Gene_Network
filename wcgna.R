#setwd("/home/bioinfoadmin/vagmita/lsbio_colon/counts_11_27_normal/onto")
library(edgeR)
library(RColorBrewer)
library(scatterplot3d)
library("tidyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("dplyr", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("DESeq2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library(DESeq)
library(gplots)
library(ggplot2)
library(plyr)
library(readxl)
library(factoextra)

####filter low count
#####filter low variance
#####take only DE genes
#####Log2CPM
####normalizaaton may or may not
#####SPEARMAN CORELATION, Bi Weight, EUCLIDEAN, MUTUAL INFORMATION ETC for non linear

####SIgn of correlation, if we have 2 set of gene + - , should we gp together where one is opp then other, take abs value 
####if sign keep then, (cor+1)/2
############how to eliminate corelation whic are only by chance 
####Sigmoid transformation high up and low low, 1/1+e-^x
#####Power transformation. preserve only strong corelation by riassing to higher power.

#counts_all_colon <- read.csv("pure_gene_colon.csv", header=TRUE, stringsAsFactors = TRUE, row.names=1)
counts_mapk_colon <- read.csv("../counts_11_27_normal/onto/mapk_data1.csv", header=TRUE, row.names=1)
rownames(counts_mapk_colon) <- counts_mapk_colon$Genedid
counts_mapk_colon <- counts_mapk_colon[,-1]
#counts_all_colon <- counts_all_colon[rowSums(counts_all_colon >20) >=1,]
countF_mapk_colon <- counts_mapk_colon[rowSums(counts_mapk_colon >20) >=1,]
head(countF_mapk_colon)


#colnames(countF_all_colon) <- c("Ca1_1","Ca1_2","Ca2_3","SAA_4","SAA_5","SAA_6","TAA_7","TAA_8","TAA_9","TAA_10","ile_11","Ca1_1","Ca1_2","Ca2_3","SAA_4","SAA_5","SAA_6","TAA_7","TAA_8","TAA_9","TAA_10","ile_11","ANA_12","ANA_12","ANB_13","ANB_13","ANC_14","ANC_14")

#counts_all_colon <- counts_all_colon[rowSums(counts_all_colon >20) >=1,]



#runinfo <- read.csv("matrix")
#######################Diferential Exp Anlysis############################
(condition <- factor(c("Ca1_1","Ca1_2","Ca2_3","SAA_4","SAA_5","SAA_6","TAA_7","TAA_8","TAA_9","TAA_10","ile_11","Ca1_1","Ca1_2","Ca2_3","SAA_4","SAA_5","SAA_6","TAA_7","TAA_8","TAA_9","TAA_10","ile_11","ANA_12","ANA_12","ANB_13","ANB_13","ANC_14","ANC_14")))
head(colnames(countF_mapk_colon))
#runinfo <- DataFrame(runinfo)

(coldata <- data.frame(row.names=colnames(countF_mapk_colon), condition))
#seIdx <- match(coldata, runinfo$Stage)
dds <- DESeqDataSetFromMatrix(countData = countF_mapk_colon, colData = coldata, design =~condition)
#ddsCollapsed <- collapseReplicates(dds, groupby = dds$)

dds <- DESeq(dds)

colon_mapk_datatry1 <- newCountDataSet(countF_mapk_colon, condition)

colon_mapk_datatry1 <- estimateSizeFactors(colon_mapk_datatry1)
cds_full_mapk_colon <- estimateDispersions(colon_mapk_datatry1)
#cds_full_blind_colon <- estimateDispersions(colon_all_datatry1, method="blind")
vdsFull_mapk_colon <- varianceStabilizingTransformation(cds_full_mapk_colon)
select_colon_mapk <- order(rowMeans(counts(colon_mapk_datatry1)), decreasing=TRUE)[1:50]

hmcol_colon_all <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vdsFull_mapk_colon) [select_colon_mapk,], col=hmcol_colon_all, trace="none", margin=c(10,6))
#heatmap(counts(countF_all_colon)[select_colon_All,], col=hmcol_colon_all, trace="none", margin=c(10,6))
heatmap.2(counts(colon_mapk_datatry1)[select_colon_mapk,], col=hmcol_colon_all, trace="none", margin=c(10,6))
#####heatmap of the samples to sample distance


####taking lof of the data
countF_mapk_colon_log <- log2(countF_mapk_colon +1)
x <- melt(as.matrix(countF_mapk_colon_log))
colnames(x) = c('Var1', 'Var2', 'value')
ggplot(x, aes(x=value, color=Var2)) + geom_density()


#####try using limma for DE to get only significant genes####

####first taking netire mapk data for coexpression####

cordist <- function(dat) {
  cor_matrix <- cor(t(dat))
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1- (dist_matrix / max(dist_matrix))
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix) /2)
}
sim_matrix_mapk  <- cordist(countF_mapk_colon_log) 
###########random sample to plot
heatmap_indices <- sample(nrow(sim_matrix_mapk), 20)
heatmap.2(t(sim_matrix_mapk[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          trace='non', dendrogram = 'row',
          xlab='Var1', ylab='Var1',
          main='Similarity matrix',
          density.info = 'none', revC=TRUE)
###############create adjaency matrix######for coexpression###
adj_matrix_mapk <- adjacency.fromSimilarity(sim_matrix_mapk, power=4, type='signed')
#convert to matrix#
gene_id <- rownames(adj_matrix_mapk)
adj_matrix_mapk <- matrix(adj_matrix_mapk, nrow=nrow(adj_matrix_mapk))
rownames(adj_matrix_mapk) <- gene_id
colnames(adj_matrix_mapk) <- gene_id

heatmap.2(t(adj_matrix_mapk[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          trace='non', dendrogram = 'row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info = 'none', revC=TRUE)

####hclust In R###reciprocal
gene_tree_mapk <- hclust(as.dist(1 - adj_matrix_mapk), method="average")

###use cuttree to break into smaller module###
module_labels_mapk <- cutreeDynamicTree(dendro=gene_tree_mapk, minModuleSize = 15, deepSplit = TRUE )
module_colors_mapk <- labels2colors(module_labels_mapk)

####include color version of module colors for better assignermnt in cytoscape
#gene_info <- cbind(gene_id, module=module_colors_mapk)

#gene_info$color_rgb <- col2hex(gene_info[)
###try graphml####
export_network_to_graphml <- function (adj_matrix_mapk, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE)
  {
    library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='network.graphml'
  }
  
  max_edges <- max_edge_ratio * nrow(adj_matrix_mapk)
  
  edge_to_total_ratio <- max_edges / length(adj_matrix_mapk)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_matrix_mapk), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_matrix_mapk), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_matrix_mapk[abs(adj_matrix_mapk) < threshold] <- 0
  
  # Drop any genes with no edges (TODO: Make optional)
  orphaned <- (colSums(adj_matrix_mapk) == 0)
  adj_matrix_mapk <- adj_matrix_mapk[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_matrix_mapk == 0
  is_negative <- adj_matrix_mapk < 0
  
  adj_matrix_mapk <- (abs(adj_matrix_mapk) - threshold) / (max(adj_matrix_mapk) - threshold)
  adj_matrix_mapk[is_zero] <- 0
  adj_matrix_mapk[is_negative] <- -adj_matrix_mapk[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_matrix_mapk), sum(adj_matrix_mapk > 0)))
  }
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_matrix_mapk, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_matrix_mapk[adj_matrix_mapk != 0] <- 1
    g <- graph.adjacency(adj_matrix_mapk, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list] 
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}

g <- export_network_to_graphml(adj_matrix_mapk, filename='~/network.graphml',
                               threshold=0.7)
g
source("https://bioconductor.org/biocLite.R")
biocLite("RCy3")
biocLite("RCytoscape")
