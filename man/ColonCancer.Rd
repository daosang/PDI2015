\name{ColonCancer}
\alias{ColonCancer}
\docType{data}
\title{
Colon Cancer Dataset
}
\description{
Data pertaining to the article ‘Broad patterns of gene expression revealed by clustering of tumor and normal colon tissues probed by oligonucleotide arrays’.
}
\usage{data("ColonCancer")}
\format{
  The format is:
 num [1:62, 1:2001] 1 2 1 2 1 2 1 2 1 2 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:62] "V1" "V2" "V3" "V4" ...
  ..$ : chr [1:2001] "cls" "" "" "" ...
}
\details{
The matrix I2000 contains the expression of the 2000 genes with highest minimal intensity across the 62 tissues. The genes are placed in order of descending minimal intensity. Each entry in I2000 is a gene intensity derived from the ~20 feature pairs that correspond to the gene on the chip, derived using the filtering process described in the ‘materials and methods’ section. The data is otherwise unprocessed (for example it has not been normalized by the mean intensity of each experiment).

The file ‘names’ contains the EST number and description of each of the 2000 genes, in an order that corresponds to the order in I2000. Note that some ESTs are repeated which means that they are tiled a number of times on the chip, with different choices of feature sequences. The descriptions UMGAP, HSAC07 and I correspond to control RNAs spiked with each experiment.

The identity of the 62 tissues is given in file tissues. The numbers correspond to patients, a positive sign to a normal tissue, and a negative sign to a tumor tissue.

The clustering algorithm (matlab version 5.1) is given in file cluster. This program accepts a matrix of input data, where each object to be clustered is represented by a column. It outputs three variables: Ord, Num and BetaVal. Ord contains the post-clustering order of objects along the binary tree. Additional information about the binary tree is found in Num and BetaVal. Num contains the sizes of the clusters at each splitting, so that in the second row of Num are two nonzero entries corresponding to the sizes of the two clusters in the first division in the tree, the third row has 4 entries, etc. The matrix BetaVal contains the b values of the cluster splits (see ‘materials and methods’).
}
\source{
http://genomics-pubs.princeton.edu/oncology/affydata/index.html
}
\references{
‘Broad patterns of gene expression revealed by clustering of tumor and normal colon tissues probed by oligonucleotide arrays’

U. Alon, N. Barkai, D. A. Notterman, K. Gish, S. Ybarra, D. Mack, and A. J. Levine,

Proc. Natl. Acad. Sci. USA, Vol. 96, Issue 12, 6745-6750, June 8, 1999
}
\examples{
data(ColonCancer)
## maybe str(ColonCancer) ; plot(ColonCancer) ...
}
\keyword{datasets}
