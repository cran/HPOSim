\name{getTermAncestors}
\alias{getTermAncestors}
\title{ Get Ancestors of HPO Terms }
\description{
  Given a list of HPO terms, the function returns the ancestors of each term.
}
\usage{
  getTermAncestors(hpolist, verbose = FALSE)
}
\arguments{
  \item{hpolist}{ character vector of HPO IDs }
  \item{verbose}{ print out some information }
}
\value{
  List with entry names for each HPO ID. Each entry contains a charcter vector with the ancestor HPO IDs
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{ 
  \code{\link{getTermOffsprings}} 
  \code{\link{getTermChildren}} 
  \code{\link{getTermParents}} 
}
\examples{
  terms<-c("HP:0000002","HP:0002719")
  result<-getTermAncestors(terms)
  print(result)
}
\keyword{ manip }
