setwd('~/GITCLONES/2017-08-18-cytof-jingping')
source('~/GITCLONES/cytof/fcs_processing_pipeline.R')
load("/home/x189259/GITCLONES/2017-07-25-wgcna-clusterinag-jingping/stopsle.rsem_w_HCSLELN/z_further analysis/HCSLELN_teston_combined-immune-rsem/gene_logfc.RData")
#for all packs
load("/home/x189259/GITCLONES/combined-immune-rnaseq-rsem/annotation-metadata.rda")
install.all.packs<- function(){
install.packages("glmnet") 
install.packages("pamr") 
install.packages("ggplot2")
install.packages("survival")
library("glmnet")
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
biocLite("impute")
install.packages("samr")
install.packages("brew")
install.packages("shiny")
install.packages("devtools")
library("devtools")
install_github("nolanlab/Rclusterpp")
install_github('nolanlab/citrus')
library('citrus')
source("http://bioconductor.org/biocLite.R")
biocLite()
setRepositories()
biocLite("FlowSOM")
biocLite("flowMeans")
biocLite("SamSPECTRAL")
biocLite("cytofkit")

library(cytofkit)
library(SamSPECTRAL)
library(FlowSOM)
library(flowMeans)
  install.packages("XML")
  install_github("nolanlab/cytofCore")
  library('cytofCore')
  install.packages('NbClust')
  library(NbClust)
  
  }
#
gene_logfc <- t(gene_logfc)





###basic.cell
setwd("~/GITCLONES/combined-immune-rnaseq-rsem")
source("load.R")
postImmuneLoad()

#data parsing
bc_meta <- immune.metadata[immune.metadata$basic.cells != "",]
basic_cell_types <- split(bc_meta, f = bc_meta$basic.cells)
basic_name_list <- names(basic_cell_types)
basic.logfc <- gene_logfc[make.names(rownames(gene_logfc)) %in% make.names(bc_meta[,3]),]



#this function need a expression matrix with gene name as colnames, sample names as rownames
basic.tsne.result <- cytof_dimReduction(basic.logfc, method = c("tsne", "pca", "isomap", "diffusionmap", "NULL"),
                                distMethod = "euclidean", out_dim = 2, sampleSeed = 42, num.cores = 0,
                                isomap_k = 5, isomap_ndim = NULL, isomapFragmentOK = TRUE, resultDir=".",
                                tsne_function_py = sprintf("%s/cytof/python_scripts/multicore_tsne.py",
                                                           Sys.getenv("GITCLONES")))

cluster_cytof_flowMeans <- cytof_cluster(basic.logfc,basic.tsne.result,method = 'flowMeans',
                               sampleSeed=42, cluster_k = 43, autoCluster=FALSE, resultDir=".", num.cores=0,
                               phenograph_function_py = sprintf("%s/cytof/python_scripts/multicore_clustering.py",
                                                                Sys.getenv("GITCLONES")))
cluster_cytof_FlowSOM <- cytof_cluster(basic.logfc,basic.tsne.result,method = 'FlowSOM',
                                         sampleSeed=42, cluster_k = 43, autoCluster=FALSE, resultDir=".", num.cores=0,
                                         phenograph_function_py = sprintf("%s/cytof/python_scripts/multicore_clustering.py",
                                                                          Sys.getenv("GITCLONES")))
cluster_cytof_Rclustercpp <- cytof_cluster(basic.logfc,basic.tsne.result,method = 'Rclustercpp',
                                         sampleSeed=42, cluster_k = 43, autoCluster=FALSE, resultDir=".", num.cores=0,
                                         phenograph_function_py = sprintf("%s/cytof/python_scripts/multicore_clustering.py",
                                                                          Sys.getenv("GITCLONES")))



#plot function need a dataframe combined of expr data and tsne results.
basic.input.dataframe <- as.data.frame(cbind(basic.logfc, basic.tsne.result,
                     flowMeans = cluster_cytof_flowMeans, FlowSOM = cluster_cytof_FlowSOM,
                     Rclustercpp = cluster_cytof_Rclustercpp))
 cytof_clusterPlot(data=basic.input.dataframe, xlab="tsne_1", ylab="tsne_2",
                        cluster="FlowSOM", sampleLabel = FALSE)
 cytof_clusterPlot(data=basic.input.dataframe, xlab="tsne_1", ylab="tsne_2",
                        cluster="flowMeans", sampleLabel = FALSE)
 cytof_clusterPlot(data=basic.input.dataframe, xlab="tsne_1", ylab="tsne_2",
                        cluster="Rclustercpp", sampleLabel = FALSE)
 
 
 
 
 
 
 
 
#all cell 
 
 
 ac_meta <- immune.metadata[immune.metadata$cells != "",]
all_cell_types <- split(ac_meta, f = ac_meta$cells)
 all_name_list <- names(all_cell_types)
 all.logfc <- gene_logfc[make.names(rownames(gene_logfc)) %in% make.names(ac_meta[,3]),]
 
 
all.tsne.result <- cytof_dimReduction(all.logfc, method = c("tsne", "pca", "isomap", "diffusionmap", "NULL"),
                                         distMethod = "euclidean", out_dim = 2, sampleSeed = 42, num.cores = 0,
                                         isomap_k = 5, isomap_ndim = NULL, isomapFragmentOK = TRUE, resultDir=".",
                                         tsne_function_py = sprintf("%s/cytof/python_scripts/multicore_tsne.py",
                                                                    Sys.getenv("GITCLONES")))
 
 cluster_cytof_flowMeans <- cytof_cluster(all.logfc,all.tsne.result,method = 'flowMeans',
                                          sampleSeed=42, cluster_k = 43, autoCluster=FALSE, resultDir=".", num.cores=0,
                                          phenograph_function_py = sprintf("%s/cytof/python_scripts/multicore_clustering.py",
                                                                           Sys.getenv("GITCLONES")))
 cluster_cytof_FlowSOM <- cytof_cluster(all.logfc,all.tsne.result,method = 'FlowSOM',
                                        sampleSeed=42, cluster_k = 43, autoCluster=FALSE, resultDir=".", num.cores=0,
                                        phenograph_function_py = sprintf("%s/cytof/python_scripts/multicore_clustering.py",
                                                                         Sys.getenv("GITCLONES")))
 cluster_cytof_Rclustercpp <- cytof_cluster(all.logfc,all.tsne.result,method = 'Rclustercpp',
                                            sampleSeed=42, cluster_k = 43, autoCluster=FALSE, resultDir=".", num.cores=0,
                                            phenograph_function_py = sprintf("%s/cytof/python_scripts/multicore_clustering.py",
                                                                             Sys.getenv("GITCLONES")))
 all.input.dataframe <- as.data.frame(cbind(all.logfc, all.tsne.result,flowMeans = cluster_cytof_flowMeans, FlowSOM = cluster_cytof_FlowSOM,
                                            Rclustercpp = cluster_cytof_Rclustercpp))
 cytof_clusterPlot(data=all.input.dataframe, xlab="tsne_1", ylab="tsne_2",
                   cluster="FlowSOM", sampleLabel = FALSE)
 cytof_clusterPlot(data=all.input.dataframe, xlab="tsne_1", ylab="tsne_2",
                   cluster="flowMeans", sampleLabel = FALSE)
 cytof_clusterPlot(data=all.input.dataframe, xlab="tsne_1", ylab="tsne_2",
                   cluster="Rclustercpp", sampleLabel = FALSE)
 
 
 
 ##this is a pipe line to color just specific cell types on tsne graph
setwd('~/GITCLONES/2017-08-18-cytof-jingping')
#set a dir for image saving

cell_type_func <- function(datainput,exprmatrix,namelist,celltypelist){
for(i in 1:length(namelist)){

        a <- data.frame(celltypelist[i])[,3]
        a <- make.names(a)
        #c <- a %in% colnames(gene_logfc)
        l <- a
        #check if name is in database, if not, remove it.
        # if("FALSE" %in% c){
        #   a <- as.data.frame(a)
        #   a <- a[c,]
        #   print(1)
        # }
        datainput$a <- c(rep(0,length(rownames(exprmatrix))))
 for(j in 1:length(rownames(datainput))){
          if (make.names(rownames(datainput))[j] %in% as.character(l)){
            datainput$a[j] <- "1"
            
          }
        }
        colnames(datainput)[length(colnames(datainput))] <- make.names(paste('cells_',names(celltypelist[i])))
        cytof_clusterPlot(data=datainput, xlab="tsne_1", ylab="tsne_2",title = paste(colnames(datainput)[length(colnames(datainput))]),
                          cluster=colnames(datainput)[length(colnames(datainput))], sampleLabel = FALSE,clusterColor = c("grey","red"), point_size = 0.8,addLabel = F)

          ggsave(file=paste(colnames(datainput)[length(colnames(datainput))],".png",sep = ""), plot = last_plot(), device = NULL, path = NULL,
                 scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
                 dpi = 300, limitsize = TRUE)
        }
} 
      
cell_type_func(datainput = all.input.dataframe,exprmatrix = all.logfc,namelist = all_name_list,celltypelist = all_cell_types)  
      
