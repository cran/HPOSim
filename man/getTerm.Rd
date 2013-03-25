\name{getTerm}
\alias{getTerm}
\title{ Get information of HPO terms }
\description{ Given a list of IDs of HPO term's, the function returns their names and synonyms.
}
\usage{
  getTerm(hpolist)
}
\arguments{
  \item{hpolist}{ character vector of HPO IDs }
}
\value{
  List of names and synonyms for each HPO ID. 
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{
  \code{\link{calcTermSim}}
}
\examples{
  terms<-c("HP:0000002","HP:0000012")
  result<-getTerm(terms)
  print(result)
}
\keyword{ manip }
