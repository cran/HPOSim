\name{getTermParents}
\alias{getTermParents}
\title{ Get Direct Parents of HPO Terms }
\description{
  Given a list of HPO terms, the function returns the direct parents of each term.
}
\usage{
  getTermParents(hpolist, verbose = FALSE)
}
\arguments{
  \item{hpolist}{ character vector of HPO IDs }
  \item{verbose}{ print out some information }
}
\value{
  List with entry names for each HPO ID. Each entry contains a character vector with the direct parent of HPO IDs. 
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{ 
  \code{\link{getTermOffsprings}}
  \code{\link{getTermChildren}}
  \code{\link{getTermAncestors}} 
}
\examples{
  terms<-c("HP:0007408","HP:0006538")
  result<-getTermParents(terms)
  print(result)
}
\keyword{ manip }
