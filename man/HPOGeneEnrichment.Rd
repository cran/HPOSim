\name{HPOGeneEnrichment}
\alias{HPOGeneEnrichment}
\title{HPO enrichment analysis of Gene List}
\description{
	Given a list of genes, the function returns the enriched HPO terms.
}
\usage{
	HPOGeneEnrichment(genelist,filter=5,cutoff=0.05,background=getGeneDefaultBackground())
}

\arguments{
  \item{genelist}{ character vector of gene names }
  \item{filter}{ indicates that HPO terms must have at least 'filter' genes annotated }
  \item{cutoff}{ significant cutoff for HPO enrichment analysis }
  \item{background}{	A character vector of genes used as backgroud. Default is to use all the genes annotated in HPO }
}
\value{
	Return a data.frame object with 8 columns.Details are below:
	
\item{"HPOID"}{enriched HPO IDs}
\item{"annGeneNumber"}{Gene number annotated to this HPO term in the inputed gene list}
\item{"annBgNumber"}{Gene number in the inputed gene list}
\item{"geneNumber"}{Gene number annotated to this HPO term in the backgroud list}
\item{"bgNumber"}{Gene number in the backgroud list}
\item{"odds"}{Calculated by \eqn{\frac{annGeneNumber/annBgNumber}{geneNumber/bgNumber}}{(annGeneNumber/annBgNumber)/(geneNumber/bgNumber)} }
\item{"pvalue"}{siginicance of the hypergeometric test for this HPO term}
\item{"qvale"}{multiple test correction value for pvalue using FDR}
	
}
\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{  
  \code{\link{HPODiseaseEnrichment}}
}
\examples{
  genelist<-c("6299","9241","1717","54880","85358")
  res<-HPOGeneEnrichment(genelist)
  print(res)
}
\keyword{ manip }
