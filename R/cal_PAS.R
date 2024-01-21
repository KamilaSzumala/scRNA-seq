# ---- odfiltrowanie nieaktywnych ścieżek sygnałowych ----
#' Funkcja służy odfiltrowaniu ścieżek sygnałowych, które ykazują reprezentacje
#' mniej niż 35% składających się na nią genów, w danym zbiorze danych.
#' 
#' @param pathway lista ścieżek sygnałowych, gdzie każdy jej element zawiera 
#'                  wektor z nazwami genów do niej należących.
#'                    
#' @param cells \code{ramka danych} zawierająca poziom ekspresji genów (wiersze)
#'              we wszystkich komórkach (kolumny). Wiersze powiny być nazwane
#'              według genu, który reprezentują.
#' 
#' @return Funkcja zwraca liste tylko aktywnych ścieżek w danym zbiorze danych. 
#' 
#' @export
inactive_path <- function(pathway, cells){
  
  name <- c()
  for (i in 1:length(pathway)){
    df <- cells[rownames(cells) %in% pathway[[i]],]
    if(nrow(df) != 0){
      if(round(nrow(df)/length(pathway[[i]]), digits = 2) < .65){
        print(names(pathway[i]))
        name <- append(name, names(pathway[i]))
      }
    }else{
      name <- c(name, names(pathway[i]))
    }
  }
  pathway[name] <- NULL
  
  return(pathway)
}

# ---- top 20% var -----
#' Funkcja służąca odfiltrowaniu zadanego procenta genów o najwyższej wariancji.
#' Reprezentuje ona standardowe podejście do analizy danych scRNA-seq.
#' 
#' @param cells \code{ramka danych} zawierająca poziom ekspresji genów (wiersze)
#'              we wszystkich komórkach (kolumny). Wiersze powiny być nazwane
#'              według genu, który reprezentują.
#'              
#' @param percent zadany procent do odfiltrowania genów o najwyższej wariancji.
#'                Podany jako liczba dziesiętna. 
#'                
#' @return Funkcja zwraca macierz zawierającą geny o najwyższej wariancji.
#' 
#' @export            
topVar <- function(cells, percent){
  
  top_data <- cells
  top_data$var <- apply(top_data, 1, var)
  top_data <- top_data[order(top_data$var, decreasing = T),]
  top_data <- as.data.frame(top_data[1:(percent*nrow(top_data)),1:ncol(cells)])
  
  return(top_data)
}

# ---- Algorytmy PAS ----
#' Funkcje służące obliczeniu transformacji macierzy poziomu ekspresji genów 
#' do poziomu aktywacji ścieżek sygnałowych (PAS). 
#'
#' @param pathways lista ścieżek sygnałowych, gdzie każdy jej element zawiera 
#'                  wektor z nazwami genów do niej należących.
#'                    
#' @param cells \code{ramka danych} zawierająca poziom ekspresji genów (wiersze)
#'              we wszystkich komórkach (kolumny). Wiersze powiny być nazwane
#'              według genu, który reprezentują.
#' 
#' @return funkcje zwracają macierze po transformacji wybraną metodą.
#'
#' @example 
#' \dontrun{
#' data_file <- "data/BreastCancer/BreastCancer_filtered_normalized.csv"
#' data <- read.table(data_file, sep = ",", header = T, row.names = 1)
#' 
#' tmp_env <- new.env()
#' load("data/All_path.RData", envir = tmp_env)
#' paths <- get("BC_pathways", tmp_env)
#' 
#' pathway <- inactive_path(paths, data)
#' dfPlage <- calPlage(pathway, data)
#' }
#' 
#' @export
#* Mean ---- 
calMean <- function(pathways, cells){
  
  dfMean <- data.frame()
  for (i in 1:length(pathways)){
    df <- Path_df(pathways[[i]], cells)
    dfMean <- rbind(dfMean, colMeans(df))
  }
  rownames(dfMean) <- names(pathways)
  
  return(dfMean)
}

#* PLAGE ----
calPlage <- function(pathways, cells){
  
  dfPlage <- data.frame(); 
  for (i in 1:length(pathways)){
    df <- Path_df(pathways[[i]], cells)
    pca.path <- prcomp(t(df), scale = T, center = T)
    pca.path.x <- as.data.frame(pca.path$x)
    dfPlage <- rbind(dfPlage, pca.path.x$PC1)
  }
  rownames(dfPlage) <- names(pathways) 
  colnames(dfPlage) <- colnames(cells)
  
  return (dfPlage)
}

#* Z-score ----
calZsc <- function(pathways, cells){
  
  Nor_data <- matrix(NaN, nrow(cells), ncol(cells))
  for (j in 1:nrow(cells)){
    row <- as.numeric(as.matrix(cells[j,]))
    Nor_data[j,] <- sapply(row, Zscore, mean(row), sd(row))
  }
  Nor_data <- as.data.frame(Nor_data)
  colnames(Nor_data) <- colnames(cells)
  rownames(Nor_data) <- rownames(cells)
  
  dfZsc <- data.frame()
  for (i in 1:length(pathways)){
    df <- Path_df(pathways[[i]], Nor_data)
    s <- colSums(df)/sqrt(nrow(df))
    dfZsc <- rbind(dfZsc, s)
  }
  rownames(dfZsc) <- names(pathways)
  colnames(dfZsc) <- colnames(cells)
  
  return (dfZsc)
}

#* CERNO ----
calCerno<- function(pathways, cells){
  
  Rank_data <- matrix(1, nrow(cells), ncol(cells))
  for (i in 1:ncol(as.matrix(Rank_data))){
    colRank <- c()
    col <- as.numeric(as.vector(cells[,i]))
    colRank <- rank(dplyr::desc(col), ties.method = "average")
    Rank_data[,i] <- colRank
  }
  
  Rank_data <- as.data.frame(Rank_data)
  colnames(Rank_data) <- colnames(cells)
  rownames(Rank_data) <- rownames(cells)
  
  dfCerno <- data.frame()
  for (i in 1:length(pathways)){
    df <- Path_df(pathways[[i]], Rank_data)
    
    rowAUC <- c()
    for (j in 1:ncol(as.matrix(df))){
      cell <- as.numeric(as.vector(df[,j]))
      AUC <- (nrow(df)*(nrow(Rank_data)-nrow(df)) + (nrow(df)*(nrow(df)+1)/2) - sum(cell))/(nrow(df)*(nrow(Rank_data)-nrow(df)))
      rowAUC <- c(rowAUC, AUC)
    }
    
    dfCerno <-rbind(dfCerno, rowAUC)
  }
  rownames(dfCerno) <- names(pathways)
  colnames(dfCerno) <- colnames(cells)
  
  return(dfCerno)
}

#* AUCells----
calAUCell <- function(pathways, cells){
  
  cells_rankings <- AUCell_buildRankings(as.matrix(cells))
  geneSets <- pathways
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=0.25*nrow(cells), nCores=1)
  dfAucell <- as.data.frame(cells_AUC@assays@data@listData$AUC)
  
  return(dfAucell)
}
#----DropRatio----
calDropRatio <- function(pathways, cells){
  
  dfDropRatio <- data.frame()
  for (i in 1:length(pathways)){
    df <- Path_df(pathways[[i]], cells)
    row <- colSums(((df!=0)/nrow(df)))
    dfDropRatio <- rbind(dfDropRatio, row)
  } 
  rownames(dfDropRatio) <- names(pathways)
  colnames(dfDropRatio) <- colnames(cells)
  
  return (dfDropRatio)
}

#* Sparse PCA----
calSpca <- function(pathways, cells){
  
  dfSpca <- data.frame()
  for (i in 1:length(pathways)){
    df <- Path_df(pathways[[i]], cells)
    pca.path <- spca(t(df), scale = T, center = T, max_iter = 800)
    pca.path.x <- as.data.frame(pca.path$scores)
    dfSpca <- rbind(dfSpca, t(pca.path.x$V1))
  }
  rownames(dfSpca) <- names(pathways)
  colnames(dfSpca) <- colnames(cells)
  
  return (dfSpca)
}

#* JASMINE ----
calJASMINE <- function(pathways, cells){
  
  JAS_data <- matrix(0, length(pathways), ncol(cells))
  for (c in 1:ncol(cells)){
    tmp <- cells[,c]; names(tmp) <- rownames(cells)
    tmp <- tmp[!tmp == 0]
    tmp <- rank(tmp, ties.method = "average")
    
    Jas_vec <- matrix(0,length(pathways), 1)
    for(i in 1:length(pathways)){
      rk <- tmp[pathways[[i]]]; rk <- rk[!is.na(rk)]
      Jas_vec[i] <- mean(rk)/length(tmp)
    }
    JAS_data[,c] <- Jas_vec
  }
  
  JAS_data[is.na(JAS_data)] <- 0
  
  dfJAS <- JAS_data
  for (i in 1:nrow(JAS_data)){
    dfJAS[i,] <- sapply(JAS_data[i,], norm_jas, min(JAS_data[i,]), max(JAS_data[i,]))
  }
  
  dfJAS <- as.data.frame(dfJAS)
  rownames(dfJAS) <- names(pathways)
  colnames(dfJAS) <- colnames(cells)
  
  return(dfJAS)
}

# ---- funkcje wspierające ----
Path_df <- function(path, exp_mtx){
  
  x <- exp_mtx[rownames(exp_mtx) %in% path,]
  tmp <- rowVars(as.matrix(x))
  m <- x[!(is.na(tmp)==T | tmp == 0),]
  
  return(m)
}

Zscore <- function(i, sr, sd){
  z <- (i-sr)/sd
  return(z)
}

norm_jas <- function(jas, min, max){
  score <- (jas-min)/(max- min)
  return(score)
}




