\name{getTermChildren}
\alias{getTermChildren}
\title{ Get Direct Children of HPO terms }
\description{
  Given a list of HPO terms, the function returns the direct children of each term.
}
\usage{
  getTermChildren(hpolist, verbose = FALSE)
}
\arguments{
  \item{hpolist}{ character vector of HPO IDs }
  \item{verbose}{ print out some information }
}
\value{
  List with entry names for each HPO ID. Each entry contains a character vector with the direct children of HPO IDs. 
}
\author{
Yue Deng<anfdeng@163.com>
}
\seealso{ 
  \code{\link{getTermOffsprings}}
  \code{\link{getTermParents}}
  \code{\link{getTermAncestors}} 
}
\examples{
  terms<-c("HP:0000002","HP:0002719")
  result<-getTermChildren(terms)
  print(result)
}
\keyword{ manip }

