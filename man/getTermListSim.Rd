\name{getTermListSim}
\alias{getTermListSim}
\title{ Similarity matrix for two sets of HPO terms }
\description{
  Given two sets of HPO terms, the function calculates similarity between two sets.
}
\usage{
  getTermListSim(anno1, anno2, combinemethod="funSimMax",
                    method="Resnik", IC, verbose=FALSE)
}
\arguments{
  \item{anno1}{ one set of HPO terms }
  \item{anno2}{ another set of HPO terms }
  \item{combinemethod}{ method to calculate the similarity between diseases based on , one of "max", "mean", "funSimMax", "funSimAvg" and ""BMA }
  \item{method}{ method to compute the similarity of HPO terms, one of "Resnik", "JiangConrath", "Lin", "simIC", "relevance", "GIC" and "Wang" }
  \item{IC}{ IC used to calculate similarity between terms }
  \item{verbose}{ print out some information }
}
\details{
  The method to combine the similarities between HPO terms can either be: 
  \describe{
    \item{"max"}{the maximum similarity between any two HPO terms}
    \item{"mean"}{the average similarity between any two HPO terms1}
    \item{"funSimMax"}{ the average of best matching HPO term similarities. Take the maximum of the scores achieved by assignments of HPO terms from disease 1 to disease 2 and vice versa.}
    \item{"funSimAvg"}{ the average of best matching HPO term similarities. Take the average of the scores achieved by assignments of HPO terms from disease 1 to disease 2 and vice versa.}
    \item{"BMA"}{best match average approach} 
  }
}
\value{
  the value of the similarity of two sets of HPO terms
}
\author{
  Yue Deng<anfdeng@163.com>
}
\references{ 
[1] P. W. Lord, et al., "Investigating semantic similarity measures across the Gene Ontology: the relationship between sequence and annotation," Bioinformatics, vol. 19, pp. 1275-83, Jul 1 2003.

[2] A. Schlicker, F. Domingues, J. Rahnenfuehrer, T. Lengauer, A new measure for functional similarity of gene products based on Gene Ontology, BMC Bioinformatics, 7, 302, 2006.

[3] James Z.Wang,Zhidian Du, et al. A new method to measure the semantic similarity of GO terms.Bioinformatics 2007,Vol 23,1274-1281.

}
\seealso{  
  \code{\link{calcTermSim}}
}
\examples{
  .initialize()
  IC<-get("termIC",envir=HPOSimEnv)
  getTermListSim(c("HP:0000118", "HP:0000152", "HP:0000234", "HP:0000271"),
                 c("HP:0000284", "HP:0000478", "HP:0000479", "HP:0000488"),
                 combinemethod="funSimMax", method="Resnik", IC=IC, verbose=FALSE)
}
\keyword{ manip }
