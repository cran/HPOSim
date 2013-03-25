\name{getGeneListSim}
\alias{getGeneListSim}
\title{ Pairwise Similarity for a List of Genes}
\description{
  Given a list of genes, the function calculates the pairwise similarities for any two genes in the list using different strategies.
}
\usage{
  getGeneListSim(genelist,combinemethod="funSimMax",method="Resnik",ontology="PA",normalization=FALSE,normalizationmethod="Lin",verbose=FALSE)
}
\arguments{
  \item{genelist}{ character vector of Entrez gene IDs }
  \item{combinemethod}{ method to calculate the similarity between genes based on , one of "max", "mean", "funSimMax", "funSimAvg" and ""BMA }
  \item{method}{ method to compute the similarity of HPO terms, one of "Resnik", "JiangConrath", "Lin", "simIC", "relevance" and "Wang" }
  \item{ontology}{ the ontology used for similarity calculation }
  \item{normalization}{ normalize similarities yes/no}
  \item{normalizationmethod}{ one of "sqrt", "Lin" and "Tanimoto"}
  \item{verbose}{ print out some information }
}
\details{
  The method to calculate the pairwise disease similarity between gene products can either be: 
  \describe{
    \item{"max"}{the maximum similarity between any two DO terms}
    \item{"mean"}{the average similarity between any two DO terms1[1]}
    \item{funSimMax}{ the average of best matching DO term similarities. Take the maximum of the scores achieved by assignments of DO terms from gene 1 to gene 2 and vice versa. [2]}
    \item{funSimAvg}{ the average of best matching DO term similarities. Take the average of the scores achieved by assignments of DO terms from gene 1 to gene 2 and vice versa. [2]}
    \item{"BMA"}{best match average approach [3]} 
  }
}
\value{
  n*n similarity matrix (n = number of genes)
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{  
  \code{\link{getTermListSim}}
}
\examples{
  getGeneListSim(c("2","9","25"),combinemethod="mean")
}
\keyword{ manip }


