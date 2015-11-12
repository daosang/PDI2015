#' Bi-clustering using Co-Similarity width modification
#'
#' @description   A modification of co-similarity by Hussain for Gene Expression Data ( See in the reference)
#' @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
#' @param y response variable (1 or 2)
#' @param itr number of iterations.
#' @export
#' @return Return a list of objects.
#' @references
#' Hussain S.F. \emph{Bi-Clustering Gene Expression Data Using Co-Similarity}, 7th International Conference on Advanced Data Mining and Applications (ADMA), 16-19th Dec. 2011, Beijing, China.
#' @examples
#' library(PDI2015)
#' XSim.mod(ColonCancer[, -1], ColonCancer[,1], itr = 4)
XSim.mod <- function(x, y, itr = 4){
  temps <- proc.time()
  this.call = match.call()
   act = y
  r_clusts = max(y);  # number of row (doc) clusters
  numRows = dim(x)[1]
  numCols = dim(x)[2]

  ###########################################
  # Standardization x par ligne. Formulaire (6), page 195
  MR1 = t(scale(t(x), center = TRUE, scale = TRUE))
  MR1 = as.data.frame(MR1)
  rS = as.matrix(rowSums(abs(MR1))) # sum by ligne
  rS[which(rS == 0)] <- 1       # to avoid divide by zero
  MR = as.matrix(MR1/rS)

  ################################################
  # Standardization x  par colonne. Formulaire (5), page 195
  MC1 = scale(x, center = TRUE, scale = TRUE)
  MC1 = as.data.frame(MC1)
  cS = t(as.matrix(colSums(abs(MC1)))) # sum by colum
  cS[which(cS == 0)] <- 1        # to avoid divide by zero
  tMC = matrix(1, numRows, numCols)
  for (i in 1:numCols) {
    tMC[,i] = cS[,i]
  }
  MC = as.matrix(MC1/tMC)

  ########################################################################
  # XSIM unchanged
  # initialize SR and SC
  SR = diag(numRows) # row
  SC = diag(numCols) # column

  # Iterate between SR and SC
  for (ii in 1:itr) {
    print(paste("Traitement en cours l'iteration No: ", toString(ii)))
    SR = MR %*% SC %*% t(MR)
    for (j in 1:numRows) {
      SR[j,j] = 1
    }

    SC = t(MC) %*% SR %*% MC
    for (k in 1:numCols) {
      SC[k,k] = 1
    }
  }

  #Calculate the clusters, par row only.
  if (r_clusts == 0) {
    print("You cannot create 0 row clusters. please specify a value for r_clusts (>0)");
  }

  # Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
  # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
  disSR = as.matrix(1 - data.frame(SR))

  # Hierarchical aggomerative clustering + Ward's method
  XSIM_v2008 = hclust(dist(disSR), method = "ward.D2")
  groups     = cutree(XSIM_v2008, k = r_clusts); # cut tree into "r_clusts" clusters

  # Confusion Matrix
  matrix_confusionT = table(data.matrix(act),data.matrix(groups))
  matrix_confusion  = matrix(data.frame(matrix_confusionT)[,3],ncol = r_clusts)
  matrix_confusion  = t(matrix_confusion[nrow(matrix_confusion):1,])

  maximum  = lpSolve::lp.assign(matrix_confusion,direction = "max")$objval
  accuracy =  maximum/(length(t(y))[1])
  print(accuracy)

  time = round(((proc.time() - temps)[3][[1]]), 4)
  print(paste("Time:",time, "(s)"))
  return(list(call = this.call, accuracy = accuracy, matrix.confusion = matrix_confusion, time = time, SR = SR, SC = SC))

}

############## END #########333
