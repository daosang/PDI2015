#' Bi-clustering using Co-Similarity (version in 2015)
#'
#' @import igraph
#' @description  Bi-clustering Biological Dataset using Co-similarity by Sang & Saad.
#' @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
#' @param y response variable (1 or 2)
#' @param itr number of iterations.
#' @export
#' @return Return a list of objects.
#' @examples
#' library(PDI2015)
#' # Test with Colon Cancer dataset
#' res = XSim2015(ColonCancer[, -1], ColonCancer[,1], itr = 4)
#'
#' # Test News Groups dataset
#' res1 = XSim2015(NG20SMI$M2[[1]],  NG20SMI$M2[["class"]],  itr = 4)
#' res2 = XSim2015(NG20SMI$M2[[2]],  NG20SMI$M2[["class"]],  itr = 4)
#' res3 = XSim2015(NG20SMI$M5[[1]],  NG20SMI$M5[["class"]],  itr = 4)
#' res4 = XSim2015(NG20SMI$M10[[1]], NG20SMI$M10[["class"]], itr = 4)
#' res5 = XSim2015(NG20SMI$NG1[[1]], NG20SMI$NG1[["class"]], itr = 4)
#' res6 = XSim2015(NG20SMI$NG2[[1]], NG20SMI$NG2[["class"]], itr = 4)
#' res7 = XSim2015(NG20SMI$NG3[[1]], NG20SMI$NG3[["class"]], itr = 4)
XSim2015 <- function(x, y, itr = 2) {

  temps <- proc.time()
  this.call = match.call()

  x = as.data.frame(x)
  M = x
  r_clusts = max(y) # number of row (doc) clusters
  act = y
  dimM_1 = dim(x)[1]
  dimM_2 = dim(x)[2]

  # Normalisation M par ligne. Ensuite, sauvegarder dans MR
  rS = as.matrix(rowSums(M)) # sum by ligne
  rS[which(rS == 0)] <- 1 # to avoid divide by zero
  MR = as.matrix(M/rS)

  # Normalisation M par colonne. Sauvagarder dans MC
  cS = t(as.matrix(colSums(M))); # sum by colum
  cS[which(cS == 0)] <- 1 # to avoid divide by zero
  tMC = matrix(1, dimM_1, dimM_2)
  for (i in 1:dimM_2) {
    tMC[,i] = cS[,i]
  }
  MC = as.matrix(M/tMC)

  # initialize SR and SC
  SR = diag(dimM_1)
  SC = diag(dimM_2)

  # Iterate between SR and SC
  for (ii in 1:itr) {
    print(paste("Traitement en cours l'iteration No: ", toString(ii)));
    SR = MR %*% SC %*% t(MR)
    diag(SR) = 1
    SR = SR ^ 0.5

    SC = t(MC) %*% SR %*% MC
    diag(SC) = 1
    SC = SC ^ 0.5
  }

  #Calculate the clusters, par row only.
  if (r_clusts == 0) {
    print("You cannot create 0 row clusters. please specify a value for r_clusts (>0)");
  }

  gKernel <- igraph::graph.adjacency(SR, weighted = TRUE, mode = "undirected");
  gLaplacianKernel = igraph::graph.laplacian(gKernel, normalized = TRUE, sparse = FALSE)
  row.names(gLaplacianKernel) = NULL

  # The Moore-Penrose pseudoinverse of the Laplacian matrix of the graph
  gPlusKernel = corpcor::pseudoinverse(gLaplacianKernel)

  # Formulaire (3) dans l'article "Graph Nodes Clustering based on the Commute-Time Kernel", 2006 - Luh Yen, Francois Fouss
  # Les auteurs proposent alpha = 7, mais je teste et je reconnais que 4 est mieux que 7
  SRx = 1/(1 + exp(4*gPlusKernel/sd(gPlusKernel)))

  # Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
  # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
  disSR = as.matrix(1 - data.frame(SRx))

  # Hierarchical aggomerative clustering + Ward's method
  XSIM_v2008 = hclust(dist(disSR),method = "ward.D2")
  groups     = cutree(XSIM_v2008, k = r_clusts); # cut tree into "r_clusts" clusters

  # Confusion Matrix
  matrix_confusionT = table(data.matrix(act),data.matrix(groups))
  matrix_confusion  = matrix(data.frame(matrix_confusionT)[,3], ncol = r_clusts)
  matrix_confusion  = t(matrix_confusion[nrow(matrix_confusion):1,])

  library(lpSolve);
  maximum  = lp.assign(matrix_confusion,direction = "max")$objval
  accuracy =  maximum/(length(t(act))[1])
  print(accuracy)
  print("DONE.")

  time = round(((proc.time() - temps)[3][[1]]), 4)
  print(paste("Time:",time, "(s)"))
  return(list(call = this.call, accuracy = accuracy, matrix.confusion = matrix_confusion, time = time, SR = SR, SC = SC))

}
