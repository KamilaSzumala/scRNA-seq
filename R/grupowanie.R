# ---- K-means ----
#' Funkcja służy do wykonania grupowania metodą k-średnich dla optymalnej  
#' ilości klastrów wyliczanej na podstawie WSS oraz wykresu łokciowego.
#' Wymaga interakcji użytkownika w postani wpisaina odpowiedniej ilości 
#' grup po analizie wykresu.
#' 
#' @param PASmtx Macierz poziomu aktywacji ścieżek sygnałowych. Kolumny mają 
#'               reprezentować obserwacje (komórki), a wiersze ścieżki sygnałowe.
#' @param k Maksymalna liczba klastrów.
#' @param seed Dowolna wartość służąca możliwości odtworzenia wyników 
#'            przy randomizacji warunków początkoweych algorytmu.
#'            
#'            
#' @returns Funkcja zwraca wektor zawierający etykiety klas uzyskanych 
#'          W wyniku grupowania.
#'          
#' @export  
cluster_kmeans <- function(PASmtx, k, seed){
  
  dataframe <- apply(PASmtx, 1, scale)
  set.seed(seed)
  print(fviz_nbclust(dataframe, FUNcluster = kmeans, method = "wss", k.max = k, iter.max = 80))
  n1 <- scan()
  model2 <- kmeans(dataframe,n1)
  
  return(model2)
}

# ---- Grupowanie hierarchiczne ----
#' Funkcja służy do wykonania grupowania metodą grupowania hierarchicznego
#' dla optymalnej ilości klastrów wyliczanej na podstawie WSS oraz wykresu łokciowego.
#' Wymaga interakcji użytkownika w postani wpisaina odpowiedniej ilości 
#' grup po analizie wykresu.
#' 
#' @param PASmtx Macierz poziomu aktywacji ścieżek sygnałowych. Kolumny mają 
#'               reprezentować obserwacje (komórki), a wiersze ścieżki sygnałowe.
#' @param k Maksymalna liczba klastrów.
#'
#'
#' @returns Funkcja zwraca wektor zawierający etykiety klas uzyskanych 
#'          W wyniku grupowania.
#'          
#' @export  
cluster_hclust <- function(PASmtx, k){
  
  df <- apply(PASmtx, 1, scale)
  res.dis <- dist(df, method = "euclidean")
  print(fviz_nbclust(df, FUNcluster = hcut, method = "wss", diss = res.dis, k.max = k, verbose = T))
  n <- scan()
  model2 <- cutree(hclust(res.dis, method = "ward.D2"), n)

  return(model2)
}

# ---- Klastrowanie Louvain ----
#' Funkcja służy do wykonania grupowania metodą Louvain dla optymalnej  
#' ilości klastrów wyliczonej na podstawie modalnośći grafu.
#' 
#' @param PASmtx Macierz poziomu aktywacji ścieżek sygnałowych. Kolumny mają 
#'               reprezentować obserwacje (komórki), a wiersze ścieżki sygnałowe.
#' @param k Maksymalna liczba klastrów.
#' @param seed Dowolna wartość służąca możliwości odtworzenia wyników 
#'            przy randomizacji warunków początkoweych algorytmu.
#'            
#'            
#' @returns Funkcja zwraca wektor zawierający etykiety klas uzyskanych 
#'          W wyniku grupowania.
#'          
#' @export  
cluster_louvain <- function(PASmtx, k, seed){
  
  set.seed(seed)
  pca <- prcomp(t(PASmtx), scale. = T, center = T)
  pca.x <- as.data.frame(pca$x)
  
  graph_k10 <- scran::buildSNNGraph(t(pca.x), k = k)
  clust_k10_louvain <- igraph::cluster_louvain(graph_k10)$membership
  
  print(table(clust_k10_louvain))
  return(as.vector(clust_k10_louvain))
}

# ---- obliczanie jakości grupowania ----
#' Funkcja obliczenia takie metryk jak ARI, NMI czy FMI przy użyciu pakietu *mclustcomp*.
#' Metryki ty służą ocenie jakości klasteryzacji w porównaniu do etykiet referencyjych.
#' 
#' @param K_tool_opti Lista zawierająca wektory etykiet uzyskanych w wyniku
#'                    grupowania K-means, dla wszytskich metod transformacji PAS.
#' @param H_tool_opti Lista zawierająca wektory etykiet uzyskanych w wyniku
#'                    grupowania hierarchicznego, dla wszytskich metod 
#'                    transformacji PAS.
#' @param L_tool Lista zawierająca wektory etykiet uzyskanych w wyniku
#'               grupowania Louvain, dla wszytskich metod 
#'               transformacji PAS.
#' @param Parc Lista zawierająca wektory etykiet uzyskanych w wyniku
#'             grupowania PARC, dla wszytskich metod 
#'             transformacji PAS.
#' @param true_lab Wektor zawierający etykiety referencyjne, nadane przez eksperta.
#' @param metric Ciąg znaków wskazujący na metrykę która ma zostać obliczona 
#'              (zgodnie z pakietem *mclustcomp*).
#'              
#' 
#' @return Funkcja zwraca \code{ramkę danych} zawierającą wartości metryki
#'         dla każdej metody grupowania oraz każdej metody transformacji PAS.
#' 
clust_quality <- function(K_tool_opti, H_tool_opti, L_tool, Parc, true_lab, metric){
  
  metric_mtx <- data.frame(kmeans_opti = as.numeric(sapply(K_tool_opti, mclustcomp, true_lab, c(metric))[2,]),
                            hclust_opti = as.numeric(sapply(H_tool_opti, mclustcomp, true_lab, c(metric))[2,]),
                            louvain = as.numeric(sapply(L_tool, mclustcomp, true_lab, c(metric))[2,]),
                            PARC = as.numeric(sapply(Parc, mclustcomp, true_lab, c(metric))[2,]))

  rownames(metric_mtx) <- c("TopData", "AUCell", "Cerno", "DropRatio", 
                            "GSVA", "Mean", "Plage", "S-PCA", 
                            "Z-score", "JASMINE", "ssGSEA")
  
  return(metric_mtx)
}

