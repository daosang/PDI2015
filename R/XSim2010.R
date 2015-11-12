#' Bi-clustering using Co-Similarity (version in 2010)
#'
#' @import lpSolve
#' @description  The χ-Sim algorithm is a co-similarity based approach which builds on the idea of simultaneously generating the similarity matrices SR (between documents) and SC (between words), each of them built on the basis of the other. This method is to implement the algorithm X-SIM in 2010 developped by Feward Sye Hussain and Gilles Bission (See in the reference). Here we propose a generalization of this approach by introducing a notion of pseudo-norm and a pruning algorithm.
#' @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
#' @param y response variable (1 or 2)
#' @param itr number of iterations.
#' @param k notion of pseudo-norm
#' @param p It's pruning. We set to 0 the p percent of the lowest similarity values in the similarity matrices SR (between documents) and SC (between words).
#' @export
#' @return Return a list of objects.
#' @references
#' Hussain F., Grimal C., Bisson G. \emph{An improved Co-Similarity Measure for Document Clustering}, 9th IEEE International Conference on Machine Learning and Applications (ICMLA), 12-14th Dec. 2010, Washington, United States.
#' @examples
#' library(PDI2015)
#' # Test with Colon Cancer dataset
#' res = XSim2010(ColonCancer[, -1], ColonCancer[,1], itr = 4, k = 0.8, p =0)
#'
#' # Test News Groups dataset
#' res1 = XSim2010(NG20SMI$M2[[1]],  NG20SMI$M2[["class"]],  itr = 4, k = 0.8, p =0)
#' res2 = XSim2010(NG20SMI$M2[[2]],  NG20SMI$M2[["class"]],  itr = 4, k = 0.8, p =0)
#' res3 = XSim2010(NG20SMI$M5[[1]],  NG20SMI$M5[["class"]],  itr = 4, k = 0.8, p =0)
#' res4 = XSim2010(NG20SMI$M10[[1]], NG20SMI$M10[["class"]], itr = 4, k = 0.8, p =0)
#' res5 = XSim2010(NG20SMI$NG1[[1]], NG20SMI$NG1[["class"]], itr = 4, k = 0.8, p =0)
#' res6 = XSim2010(NG20SMI$NG2[[1]], NG20SMI$NG2[["class"]], itr = 4, k = 0.8, p =0)
#' res7 = XSim2010(NG20SMI$NG3[[1]], NG20SMI$NG3[["class"]], itr = 4, k = 0.8, p =0)
XSim2010 <- function(x, y, itr = 4, k = 0.8, p = 0) {
  temps <- proc.time()
  this.call = match.call()

  x = as.data.frame(x)
  r_clusts = max(y) # number of row (doc) clusters
  dimM_1 = dim(x)[1]
  dimM_2 = dim(x)[2]
  Mk = as.matrix(x ^ k)


  # initialize SR and SC
  SR = diag(dimM_1)
  SC = diag(dimM_2)

  # Iterate between SR and SC
  for (ii in 1:itr) {
    print(paste("Traitement en cours l'iteration No: ", toString(ii)));
    SR = Mk %*% SC %*% (t(Mk))
    SR = SR ^ (1/k)
    norSR = (diag(SR) %*% t(diag(SR))) ^ 0.5
    SR = SR/norSR
    SR = prune(SR, p = p)

    SC = t(Mk) %*% SR %*% Mk
    SC = SC ^ (1/k)
    row.names(SC) = NULL
    norSC = (diag(SC) %*% t(diag(SC))) ^ 0.5
    SC = SC/norSC
    SC = prune(SC, p = p)
  }
  #Calculate the clusters, par row only.
  if (r_clusts == 0) {
    print("You cannot create 0 row clusters. please specify a value for r_clusts (>0)");
  }

  # Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it.
  # http://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
  disSR = as.matrix(1 - data.frame(SR))

  # Hierarchical aggomerative clustering + Ward's method
  XSIM_v2008 = hclust(dist(disSR),method = "ward.D2")
  groups     = cutree(XSIM_v2008, k = r_clusts) # cut tree into "r_clusts" clusters

  # Confusion Matrix
  matrix_confusionT = table(data.matrix(y),data.matrix(groups))
  matrix_confusion  = matrix(data.frame(matrix_confusionT)[,3],ncol = r_clusts)
  matrix_confusion  = t(matrix_confusion[nrow(matrix_confusion):1,])

  maximum  = lp.assign(matrix_confusion,direction = "max")$objval
  accuracy =  maximum/(length(t(y))[1])
  print(accuracy)
  print("DONE.")

  time = round(((proc.time() - temps)[3][[1]]), 4)
  print(paste("Time:",time, "(s)"))
  return(list(call = this.call, accuracy = accuracy, matrix.confusion = matrix_confusion, time = time, SR = SR, SC = SC))
}


# Cette fonction a pour assigner p% les valeurs au plus bas, vers zero.
prune <- function(A, p = 0){
  dimA = dim(A)[1]
  A1 = as.array(c(A))                # Conversation la matrice A en un tableau (array) A1
  S1 = sort.int(A1,index.return = T) # Arranger le tableau A1, en considerant des indexations
  Idx = as.array(S1$ix)              # Separation des indexations
  num = trunc(p*dim(Idx))            # p% est combien d'elements de la matrice A? c'est "num"

  for (i in 1:num) {A1[Idx[i]] = 0}     # On va mettre p% les elements les plus petits en 0

  CR = matrix(0,dimA,dimA)            # Rebuilding la matrice A, apres d'avoir supprim? p% les elements les plus petits
  for (i in 1:dimA) {
    for (j in 1:dimA) {
      CR[i,j] = A1[(j - 1)*dimA + i]
    }
  }
  return(CR)
}
