col_values_C = c("1" = "cornflowerblue",
                 "2" = "deeppink4",
                 "3" = "cyan3",
                 "4"= "darkgoldenrod",
                 "5" = "lightpink4",
                 "6" = "chocolate3",
                 "7" = "seagreen4",
                 "8" = "firebrick3",
                 "9" = "gray29",
                 "10" = "mediumpurple2",
                 "11" = "dodgerblue3",
                 "12" = "indianred2",
                 "13" = "lightblue4",
                 "14" = "orchid1",
                 "15" = "chartreuse",
                 "16" = "cadetblue4",
                 "17" = "khaki2",
                 "18" = "burlywood2",
                 "19" = "maroon1",
                 "20" = "navy")

#----Bone Marrow----
meta <- read.table("data/BoneMarrow/BoneMarrow_meta.csv", sep = ",", header = T, row.names = 1)

col_values_B <- c("dendritic cell" = "cornflowerblue",
                  "endothelial cell" = "deeppink4",
                  "erythroid progenitor cell" = "cyan3",
                  "granulocyte"= "darkgoldenrod",
                  "granulocyte monocyte progenitor cell" = "lightpink4",
                  "hematopoietic multipotent progenitor cell" = "chocolate3",
                  "mast cell" = "seagreen4",
                  "mature B cell" = "firebrick3",
                  "megakaryocyte" = "gray29",
                  "monocyte" = "mediumpurple2",
                  "natural killer cell" = "dodgerblue3",
                  "precursor B cell" = "indianred2",
                  "pro-B cell" = "lightblue4",
                  "progenitor cell" = "orchid1")

true_lab <- meta$Cell.type.Ontology
true_lab <- sapply(true_lab, switch, "dendritic cell" = 1, "endothelial cell" = 2, 
                   "erythroid progenitor cell" = 3, "granulocyte" = 4, 
                   "granulocyte monocyte progenitor cell" = 5, 
                   "hematopoietic multipotent progenitor cell" = 6,, "mast cell" = 7,
                   "mature B cell" = 8, "megakaryocyte" = 9, "monocyte" = 10,
                   "natural killer cell" = 11, "precursor B cell" = 12,
                   "pro-B cell" = 13, "progenitor cell" = 14)

true_lab <- as.integer(true_lab)
names(true_lab) <- meta$cell.ID

#----Breast Cancer ----
meta <- read.table("data/BreastCancer/BreastCancer_meta.csv", sep = ",", header = T, row.names = 1)

col_values_B = c("HER +" = "cornflowerblue",
                 "Luminal A" = "deeppink4",
                 "Luminal B" = "cyan3",
                 "Normal tissue"= "darkgoldenrod",
                 "TNBC" = "lightpink4")

true_lab <- meta$group
true_lab <- sapply(true_lab, switch, "HER +" = 1, "Luminal A" = 2, 
                   "Luminal B" = 3, "Normal tissue" = 4, "TNBC" = 5)

true_lab <- as.integer(true_lab)
names(true_lab) <- meta$cell.ID

# ---- COVID ---- 
meta <- read.table("data/COVID/COVID_meta.csv", sep = ",", header = T, row.names = 1)

col_values_B = c("B cell" = "cornflowerblue",
                 "erythroid lineage cell" = "deeppink4",
                 "monocyte" = "cyan3",
                 "neutrophil"= "darkgoldenrod",
                 "platelet" = "lightpink4",
                 "T cell" = "chocolate3")

true_lab <- meta$Cell.type.Ontology
true_lab <- sapply(true_lab, switch, "B cell" = 1, "erythroid lineage cell" = 2, 
                   "monocyte" = 3, "neutrophil" = 4, "platelet" = 5, "T cell" = 6)

true_lab <- as.integer(true_lab)
names(true_lab) <- meta$cell.ID

#---- Liver ----
meta <- read.table("data/Liver/Liver_meta.csv", sep = ",", header = T, row.names = 1)

col_values_B = c("B cell" = "cornflowerblue",
                 "cycling cell" = "deeppink4",
                 "endothelial cell of hepatic sinusoid" = "cyan3",
                 "hepatocyte"= "darkgoldenrod",
                 "macrophage" = "lightpink4",
                 "natural killer cell" = "chocolate3",
                 "T cell" = "seagreen4",
                 "cholangiocyte" = "firebrick3",
                 "endothelial cell of vascular tree" = "gray29",
                 "hematopoietic stem cell" = "mediumpurple2",
                 "Kupffer cell" = "dodgerblue3",
                 "plasma cell" = "lightblue4")

true_lab <- meta$Cell.type.org
true_lab <- sapply(true_lab, switch, "B cell" = 1, "cholangiocyte" = 2, 
                   "cycling cell" = 3, "endothelial cell of hepatic sinusoid" = 4,
                   "endothelial cell of vascular tree" = 5, 
                   "hematopoietic stem cell" = 6, "hepatocyte" = 7, "Kupffer cell" = 8,
                   "macrophage" = 9, "natural killer cell" = 10, "plasma cell" = 11,
                   "T cell" = 12)

true_lab <- as.integer(true_lab)
names(true_lab) <- meta$ID

# ---- PBMC ----
meta <- read.table("data/PBMC/PBMC_meta.csv", sep = ",", header = T, row.names = 1)

col_values_B = c("CD14+ monocyte" = "cornflowerblue",
                 "CD16+ monocyte" = "deeppink4",
                 "CD4+ T cell" = "cyan3",
                 "B cell"= "darkgoldenrod",
                 "Dendritic cell" = "lightpink4",
                 "Cytotoxic T cell" = "chocolate3",
                 "Natural killer cell" = "seagreen4",
                 "Plasmacytoid dendritic cell" = "firebrick3",
                 "Megakaryocyte" = "gray29")

true_lab <- meta$CellType
true_lab <- sapply(true_lab, switch, "B cell" = 1, "CD14+ monocyte" = 2, 
                   "CD16+ monocyte" = 3, "CD4+ T cell" = 4, "Cytotoxic T cell" = 5,
                   "Dendritic cell" = 6, "Megakaryocyte" = 7, "Natural killer cell" = 8,
                   "Plasmacytoid dendritic cell" = 9)

true_lab <- as.integer(true_lab)
names(true_lab) <- meta$NAME

# ---- Pancreas ----
meta <- read.table("data/Pancreas/Pancreas_meta.csv", sep = ",", header = T, row.names = 1)

col_values_B <- c("acinar" = "cornflowerblue",
                "activated_stellate" = "deeppink4",
                "alpha" = "cyan3",
                "beta"= "darkgoldenrod",
                "delta" = "lightpink4",
                "ductal" = "chocolate3",
                "endothelial" = "seagreen4",
                "gamma" = "firebrick3",
                "macrophage" = "dodgerblue3",
                "quiescent_stellate" = "gray29")

true_lab <- meta$CellType
true_lab <- sapply(true_lab, switch, "acinar" = 1, "activated_stellate" = 2, 
                   "alpha" = 3, "beta" = 4, "delta" = 5, "ductal" = 6, 
                   "endothelial" = 7, "gamma" = 8, "macrophage" = 9, 
                   "quiescent_stellate" = 10)

true_lab <- as.integer(true_lab)
names(true_lab) <- meta$index
