\name{HPODiseaseNOASubNetEnrichment}
\alias{HPODiseaseNOASubNetEnrichment}
\title{HPO enrichment analysis of disease network using NOA}
\description{
	Given a disease network, the function returns the enriched HPO terms using sub-net method of NOA.
}
\usage{
	HPODiseaseNOASubNetEnrichment(testfile,backgroundfile,filter=5,cutoff=0.05)
}
\arguments{
  \item{testfile}{ CSV format disease network, each line of which is OMIM IDs of two diseases }
  \item{backgroundfile}{ CSV format disease network used as the reference set, each line of which is OMIM IDs of two diseases }
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
  \code{\link{HPODiseaseNOAWholeNetEnrichment}}
}
\examples{
  edges <- data.frame(node1=c("OMIM:275350","OMIM:212750","OMIM:212750",
                             "OMIM:275350","OMIM:212750","OMIM:260400"),
                      node2=c("OMIM:100050","OMIM:105600","OMIM:249270",
                              "OMIM:300751","OMIM:250250","OMIM:105600"))
  graph<-igraph::graph.data.frame(edges,directed=FALSE)
  write.csv(edges,file="HPODiseaseNOASubNetEnrichment-testnetwork.csv",row.names=FALSE) 
 
  edges <- data.frame(node1=c("OMIM:275350","OMIM:212750","OMIM:212750","OMIM:275350",
                              "OMIM:212750","OMIM:300835","OMIM:260400","OMIM:275350",
                              "OMIM:260400","OMIM:300835","OMIM:106230","OMIM:106240",
                              "OMIM:106250","OMIM:106260"),
                      node2=c("OMIM:100050","OMIM:105600","OMIM:249270","OMIM:300751",
                              "OMIM:250250","OMIM:300751","OMIM:250250","OMIM:250250",
                              "OMIM:105600","OMIM:105600","OMIM:106260","OMIM:106230",
                              "OMIM:106240","OMIM:106250"))
  graph<-igraph::graph.data.frame(edges,directed=FALSE)
  write.csv(edges,file="HPODiseaseNOASubNetEnrichment-backgroundnetwork.csv",row.names=FALSE) 

  HPODiseaseNOASubNetEnrichment("HPODiseaseNOASubNetEnrichment-testnetwork.csv",
                                "HPODiseaseNOASubNetEnrichment-backgroundnetwork.csv",cutoff=0.1)
}
\keyword{ manip }
