\name{getDiseaseListSim}
\alias{getDiseaseListSim}
\title{ Pairwise Similarity for a List of Diseases}
\description{
  Given a list of diseases, the function calculates the pairwise similarities for any two diseases in the list using different strategies.
}
\usage{
  getDiseaseListSim(diseaselist,combinemethod="funSimMax",method="Resnik",ontology="PA",normalization=FALSE,normalizationmethod="Lin",verbose=FALSE)
}
\arguments{
  \item{diseaselist}{ character vector of disease OMIM IDs }
  \item{combinemethod}{ method to calculate the similarity between diseases based on , one of "max", "mean", "funSimMax", "funSimAvg" and ""BMA }
  \item{method}{ method to compute the similarity of HPO terms, one of "Resnik", "JiangConrath", "Lin", "simIC", "relevance" and "Wang" }
  \item{ontology}{ the ontology used for similarity calculation }
  \item{normalization}{ normalize similarities yes/no}
  \item{normalizationmethod}{ one of "sqrt", "Lin" and "Tanimoto"}
  \item{verbose}{ print out some information }
}
\details{
  The combine method to calculate the pairwise disease similarity between disease can either be: 
  \describe{
    \item{"max"}{the maximum similarity between any two HPO terms}
    \item{"mean"}{the average similarity between any two HPO terms1[1]}
    \item{funSimMax}{ the average of best matching HPO term similarities. Take the maximum of the scores achieved by assignments of HPO terms from disease 1 to disease 2 and vice versa. [2]}
    \item{funSimAvg}{ the average of best matching HPO term similarities. Take the average of the scores achieved by assignments of HPO terms from disease 1 to disease 2 and vice versa. [2]}
    \item{"BMA"}{best match average approach [3]} 
  }
}
\value{
  n*n similarity matrix (n = number of diseases)
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{  
  \code{\link{getDiseaseSim}}
}
\examples{
  list<-c("OMIM:101900","OMIM:102000","OMIM:143470")
  getDiseaseListSim(list,combinemethod="funSimAvg")
}
\keyword{ manip }


