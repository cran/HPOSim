\name{getTermOffsprings}
\alias{getTermOffsprings}
\title{Get Offspring of HPO terms}
\description{
  Given a list of HPO terms, the function returns the offsprings of each term.
}
\usage{
  getTermOffsprings(hpolist, verbose = FALSE)
}
\arguments{
  \item{hpolist}{ character vector of HPO IDs }
  \item{verbose}{ print out some information }
}
\value{
  List with entry names for each HPO ID. Each entry contains a character vector with the offspring HPO IDs.
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{ 
  \code{\link{getTermChildren}}
  \code{\link{getTermParents}}
  \code{\link{getTermAncestors}} 
}
\examples{
  terms<-c("HP:0000002","HP:0002719")
  result<-getTermOffsprings(terms)
  print(result)
}
\keyword{ manip }

