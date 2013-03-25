\name{getMinimumSubsumer}
\alias{getMinimumSubsumer}
\title{ Get minimum subsumer (MICA) of Two HPO Terms }
\description{
  Given two HPO terms, the function calculates the minimum subsumer(the common ancestor having the maximal information content) of two terms. 
}
\usage{
  getMinimumSubsumer(term1, term2)
}
\arguments{
  \item{term1}{ HPO term 1 }
  \item{term2}{ HPO term 2 }
}
\value{
  HPO term representing the minmum subsumber. If there is no minumum subsumer,the result is the string "NA".
}
\author{
  Yue Deng<anfdeng@163.com>
}
\examples{
  term1="HP:0000010"
  term2="HP:0000012"
  print(getMinimumSubsumer(term1,term2))
}
\keyword{ manip }
