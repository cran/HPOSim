\name{getSimWang}
\alias{getSimWang}
\title{ Semantic Similarity Between Two HPO Terms by Wang's Method }
\description{
  Given two HPO terms, this function will calculate the Wang's Semantic Similarity between them
}
\usage{
  getSimWang(term1, term2)
}

\arguments{
  \item{term1}{ one HPO term }
  \item{term2}{ another HPO term }
}
\value{ Semantic similarity.}
\references{
[1] J. Z. Wang, Z. Du, R. Payattakool, P. S. Yu, and C.-F. Chen, "A new method to measure the semantic similarity of GO terms", Bioinformatics, vol. 23, no. 10, pp. 1274-1281, May. 2007.
}
\author{
  Yue Deng <anfdeng@163.com>
}
\examples{
  getSimWang("HP:0000028","HP:0000033")
}
\keyword{ manip }
