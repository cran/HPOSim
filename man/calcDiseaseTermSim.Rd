\name{calcDiseaseTermSim}
\alias{calcDiseaseTermSim}
\title{ Semantic Similarity of Two HPO Terms Based on Diseases }
\description{
  Given two HPO terms, return semantic similarity of the terms based on the map of HPO terms and OMIM IDs.
}
\usage{
  calcDiseaseTermSim(term1, term2, method = "Resnik", verbose = FALSE)
}
\arguments{
  \item{term1}{ one HPO term }
  \item{term2}{ another HPO term }
  \item{method}{ one of "Resnik", "JiangConrath", "Lin", "simIC", "relevance" and "Wang" }
  \item{verbose}{ print some information }
}
\value{
  Semantic similarity of the two terms.
}
\author{
  Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{getDiseaseListSim}}
  \code{\link{getDiseaseSim}}
}
\examples{
  print(calcDiseaseTermSim("HP:0003796","HP:0003797"))
}

\keyword{ manip }