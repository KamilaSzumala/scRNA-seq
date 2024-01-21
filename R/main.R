# ---- wczytanie danych ----
data_file <- "data/BreastCancer/BreastCancer_filtered_normalized.csv"
meta_file <- "data/BreastCancer/BreastCancer_meta.csv"

data <- read.table(data_file, sep = ",", header = T, row.names = 1)
meta <- read.table(meta_file, sep = ",", header = T, row.name = 1)

tmp_env <- new.env()
load("data/All_path.RData", envir = tmp_env)
paths <- get("BC_pathways", tmp_env)

# ---- PCA i tSNE ----
library(ggplot2)
library(Rtsne)
library(ggpubr)

pca.data <- prcomp(t(data), scale = T, center = T)
pca.data.x <- as.data.frame(pca.data$x)

var_explained <- data.frame(x = pca.data$sdev^2 / sum(pca.data$sdev^2))
var_explained$comp <- 1:ncol(data)

# wykres łokciowy do określenia liczby najważniejszych komponent (ok 90% wariancji)
ggplot(var_explained[1:50,],aes(x = comp, y = x))+ geom_point() + theme_bw()

set.seed(123)
data.tSNE <- Rtsne(pca.data.x[,1:15], check_duplicates = F, perplexity = 20,
                   max_iter = 2500, theta = 0, pca = F, normalize = F)
tSNE <- as.data.frame(data.tSNE$Y)

# ---- obliczanie macierzy PAS ----
library(matrixStats)
library(dplyr)
library(GSEABase)
library(AUCell)
library(Biobase)
library(sparsepca)
library(GSVA)
source("R/cal_PAS.R")

# ---- filtrowanie ścieżek -----
pathway <- inactive_path(paths, data)

# ---- transformacje PAS ----
TopData <- topVar(data, 0.2)
dfMean <- calMean(pathway, data)
dfPlage <- calPlage(pathway, data)
dfZsc <- calZsc(pathway, data)
dfCerno <- calCerno(pathway, data)
dfAucell <- calAUCell(pathway, data)
dfDropRatio <- calDropRatio(pathway, data)
dfSpca <- calSpca(pathway, data)
dfJAS <- calJASMINE(pathway, data)

# GSVA
dfGsva <- gsva(as.matrix(data), pathway, kcdf = "Gaussian", mx.diff = F)
dfGsva <- as.data.frame(dfGsva)

# ssGSEA
dfssGSEA <- gsva(as.matrix(data), pathway, method = "ssgsea", kcdf = "Gaussian", mx.diff = F)
dfssGSEA <- as.data.frame(dfssGSEA)

# ---- grupowanie ----
library(factoextra)
library(cluster)
library(igraph)
library(scran)
source("grupowanie.R")

frames <- list(TopData, dfAucell, dfCerno, dfDropRatio, dfGsva, dfMean, dfPlage, 
               dfSpca, dfZsc, dfJAS, dfssGSEA)
seedsy <- c(87, 1, 4, 8, 7, 3, 2, 9, 6, 10, 11)

#* K-means ----
new_variable <- c("top_data", "AUCell", "Cerno", "DropRatio", "Gsva", "Mean",
                  "Plage", "Spca", "Zscore", "JASMINE", "ssGSEA")

for(i in 1:length(frames)){
  k_list <- cluster_kmeans(PASmtx = frames[[i]], k = 15, seed = seedsy[i])
  assign(new_variable[i], k_list$cluster)
}

#* Grupowanie hierarchiczne ----
new_variable <- c("TopModel", "AUCellModel", "CernoModel", "DropRatioModel", 
                  "GsvaModel", "MeanModel", "PlageModel", "SpcaModel", "ZscModel", 
                  "JASModel", "ssGSEAModel")

for(i in 1:length(frames)){
  h_list <- cluster_hclust(PASmtx = frames[[i]], k = 15)
  assign(new_variable[i], h_list)
}

#* Grupowanie Louvain ----
new_variable <- c("TopCom", "AUCellCom", "CernoCom", "DropRatioCom", 
                  "GsvaCom", "MeanCom", "PlageCom", "SpcaCom", "ZscoreCom", 
                  "JASCom", "ssGseaCom")

for(m in 1:length(frames)){
  l_list <- cluster_louvain(PASmtx = frames[[m]], k = 20, seed = seedsy[m])
  assign(new_variable[m], l_list)
}


# ---- Ocena jakości transformacji ----
library(clusterSim)
library(amap)

#NALEŻY WCZYTAĆ true_lab DLA ODPOWIEDNIEGO ZBIORU ZE SKRYPTU meta.R!

name <-  c("Data", "AUCell", "Cerno", "DropRatio", "GSVA", "Mean", "Plage", 
           "SparsePCA", "Z-score", "JASMINE", "ssGSEA")

DavBoul <- list()
sils<-list()
for (i in 1:length(frames)){
  dist_mtx <- Dist(t(frames[[i]]),"spearman")
  
  sil <- silhouette(true_lab, dist_mtx)
  sil<-group_by(as.data.frame.matrix(sil),cluster)
  
  db <- index.DB(t(frames[[i]]), true_lab, d = dist_mtx)
  
  sils[[i]] <- as.data.frame(summarise(sil, mean_sil=mean(sil_width)))
  DavBoul[[i]] <- db
}

names(DavBoul) <- name
names(sils) <- name

#---- ocena jakości grupowania ----
library(mclustcomp)

load("data/BreastCancer/PARC.RData")
#NALEŻY WCZYTAĆ true_lab DLA ODPOWIEDNIEGO ZBIORU ZE SKRYPTU meta.R!

JAS <- JASMINE
ssGsea <- ssGSEA
ZscoreModel <- ZscModel
ssGseaModel <- ssGSEAModel


pas_name_vec <- c("Top", "AUCell", "Cerno", "DropRatio", "Gsva", "Mean", "Plage", "Spca", "Zscore", "JAS", "ssGsea")

K_tool_opti <- lapply(c("top_data", pas_name_vec[-1]), function(x) get(x))
H_tool_opti <- lapply(paste0(pas_name_vec, "Model"), function(x) get(x))
L_tool <- lapply(paste0(pas_name_vec, "Com"), get)
Parc <- lapply(paste0(pas_name_vec, "Graf"), get)

Opti_BC_ARI <- clust_quality(K_tool_opti, H_tool_opti, L_tool, Parc, true_lab, "adjrand")
Opti_BC_NMI <- clust_quality(K_tool_opti, H_tool_opti, L_tool, Parc, true_lab, "nmi2")
Opti_BC_FMI <- clust_quality(K_tool_opti, H_tool_opti, L_tool, Parc, true_lab, "fmi")
