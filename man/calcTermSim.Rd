\name{calcTermSim}
\alias{calcTermSim}
\title{ Semantic Similarity of Two HPO Terms }
\description{
  Given two HPO terms, return semantic similarity of the terms based on the map of HPO terms and genes.
}
\usage{
  calcTermSim(term1, term2, method = "Resnik", IC, verbose = FALSE)
}
\arguments{
  \item{term1}{ one HPO term }
  \item{term2}{ another HPO term }
  \item{method}{ one of "Resnik", "JiangConrath", "Lin", "simIC", "relevance" , "GIC" and "Wang" }
  \item{IC}{ the IC used in calculating similarity}
  \item{verbose}{ print some information }
}
\value{
  Semantic similarity of the two terms.
}
\author{
  Yue Deng <anfdeng@163.com>
}
\seealso{
  \code{\link{getGeneSim}}
}
\examples{
  .initialize()
  IC<-get("termIC", envir=HPOSimEnv)
  calcTermSim("HP:0000028","HP:0000033", method = "Wang", IC, verbose = FALSE)
  IC<-get("DiseasetermIC", envir=HPOSimEnv)
  calcTermSim("HP:0000028","HP:0000033", method = "Resnik", IC, verbose = FALSE)

}

\keyword{ manip }