% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/XSim2010.R
\name{XSim2010}
\alias{XSim2010}
\title{Bi-clustering using Co-Similarity (version in 2010)}
\usage{
XSim2010(x, y, itr = 4, k = 0.8, p = 0)
}
\arguments{
\item{x}{matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.}

\item{y}{response variable (1 or 2)}

\item{itr}{number of iterations.}

\item{k}{notion of pseudo-norm}

\item{p}{It's pruning. We set to 0 the p percent of the lowest similarity values in the similarity matrices SR (between documents) and SC (between words).}
}
\value{
Return a list of objects.
}
\description{
The χ-Sim algorithm is a co-similarity based approach which builds on the idea of simultaneously generating the similarity matrices SR (between documents) and SC (between words), each of them built on the basis of the other. This method is to implement the algorithm X-SIM in 2010 developped by Feward Sye Hussain and Gilles Bission (See in the reference). Here we propose a generalization of this approach by introducing a notion of pseudo-norm and a pruning algorithm.
}
\examples{
library(PDI2015)
# Test with Colon Cancer dataset
res = XSim2010(ColonCancer[, -1], ColonCancer[,1], itr = 4, k = 0.8, p =0)

# Test News Groups dataset
res1 = XSim2010(NG20SMI$M2[[1]],  NG20SMI$M2[["class"]],  itr = 4, k = 0.8, p =0)
res2 = XSim2010(NG20SMI$M2[[2]],  NG20SMI$M2[["class"]],  itr = 4, k = 0.8, p =0)
res3 = XSim2010(NG20SMI$M5[[1]],  NG20SMI$M5[["class"]],  itr = 4, k = 0.8, p =0)
res4 = XSim2010(NG20SMI$M10[[1]], NG20SMI$M10[["class"]], itr = 4, k = 0.8, p =0)
res5 = XSim2010(NG20SMI$NG1[[1]], NG20SMI$NG1[["class"]], itr = 4, k = 0.8, p =0)
res6 = XSim2010(NG20SMI$NG2[[1]], NG20SMI$NG2[["class"]], itr = 4, k = 0.8, p =0)
res7 = XSim2010(NG20SMI$NG3[[1]], NG20SMI$NG3[["class"]], itr = 4, k = 0.8, p =0)
}
\references{
Hussain F., Grimal C., Bisson G. \emph{An improved Co-Similarity Measure for Document Clustering}, 9th IEEE International Conference on Machine Learning and Applications (ICMLA), 12-14th Dec. 2010, Washington, United States.
}

