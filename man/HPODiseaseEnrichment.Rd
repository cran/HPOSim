\name{HPODiseaseEnrichment}
\alias{HPODiseaseEnrichment}
\title{HPO enrichment analysis of Disease List}
\description{
	Given a list of disease, the function returns the enriched HPO terms.
}
\usage{
	HPODiseaseEnrichment(diseaselist,filter=5,cutoff=0.05,
                            background=getDiseaseDefaultBackground())
}

\arguments{
  \item{diseaselist}{ character vector of disease OMIM IDs }
  \item{filter}{ indicates that HPO terms must have at least 'filter' diseases annotated }
  \item{cutoff}{ significant cutoff for HPO enrichment analysis }
  \item{background}{	A character vector of diseases used as backgroud. Default is to use all the diseases annotated in HPO }
}
\value{
	Return a data.frame object with 8 columns.Details are below:
	
\item{"HPOID"}{enriched HPO IDs}
\item{"annDiseaseNumber"}{Disease number annotated to this HPO term in the inputed disease list}
\item{"annBgNumber"}{Disease number in the inputed disease list}
\item{"diseaseNumber"}{Disease number annotated to this HPO term in the backgroud list}
\item{"bgNumber"}{Disease number in the backgroud list}
\item{"odds"}{Calculated by \eqn{\frac{annDiseaseNumber/annBgNumber}{diseaseNumber/bgNumber}}{(annDiseaseNumber/annBgNumber)/(geneNumber/bgNumber)} }
\item{"pvalue"}{siginicance of the hypergeometric test for this HPO term}
\item{"qvale"}{multiple test correction value for pvalue using FDR}
	
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{  
  \code{\link{HPOGeneEnrichment}}
}
\examples{
  diseaselist<-c("OMIM:605685","OMIM:167100","OMIM:304200","OMIM:219300",cutoff=0.1)
  res<-HPODiseaseEnrichment(diseaselist)
}
\keyword{ manip }
