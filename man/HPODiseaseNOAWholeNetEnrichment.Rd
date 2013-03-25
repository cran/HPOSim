\name{HPODiseaseNOAWholeNetEnrichment}
\alias{HPODiseaseNOAWholeNetEnrichment}
\title{HPO enrichment analysis of disease network using NOA}
\description{
	Given a disease network, the function returns the enriched HPO terms using whole-net method of NOA.
}
\usage{
	HPODiseaseNOAWholeNetEnrichment(file,filter=5,cutoff=0.05)
}

\arguments{
  \item{file}{ CSV format disease network, each line of which is OMIM IDs of two diseases }
  \item{filter}{ indicates that HPO terms must have at least 'filter' edges annotated }
  \item{cutoff}{ significant cutoff for HPO enrichment analysis }
}
\value{
	Return a data.frame object with 8 columns.Details are below:
	
\item{"HPOID"}{enriched HPO IDs}
\item{"annEdgeNumber"}{Edge number annotated to this HPO term in the inputed disease network}
\item{"annBgNumber"}{Edge number in the inputed disease network}
\item{"edgeNumber"}{Edge number annotated to this HPO term in the backgroud network which is chosen as all possible links in the test network}
\item{"bgNumber"}{Edge number in the backgroud network}
\item{"odds"}{Calculated by \eqn{\frac{annEdgeNumber/annBgNumber}{edgeNumber/bgNumber}}{(annEdgeNumber/annBgNumber)/(edgeNumber/bgNumber)} }
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
  edges <- data.frame(node1=c("OMIM:275350","OMIM:212750","OMIM:212750","OMIM:275350","OMIM:212750","OMIM:300835","OMIM:260400","OMIM:275350"),
                      node2=c("OMIM:100050","OMIM:105600","OMIM:249270","OMIM:300751","OMIM:250250","OMIM:300751","OMIM:250250","OMIM:250250"))
  graph<-graph.data.frame(edges,directed=FALSE)
  write.csv(edges,file="HPODiseaseNOAWholeNetEnrichment-testwork.csv",row.names=FALSE) 
  HPODiseaseNOAWholeNetEnrichment("HPODiseaseNOAWholeNetEnrichment-testwork.csv",cutoff=0.8)
}
\keyword{ manip }
