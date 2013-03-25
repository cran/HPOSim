\name{HPOGeneNOASubNetEnrichment}
\alias{HPOGeneNOASubNetEnrichment}
\title{HPO enrichment analysis of gene network using NOA}
\description{
	Given a gene network, the function returns the enriched HPO terms using sub-net method of NOA.
}
\usage{
	HPOGeneNOASubNetEnrichment(testfile,backgroundfile,filter=5,cutoff=0.05)
}
\arguments{
  \item{testfile}{ CSV format gene network, each line of which is Entrez IDs of two genes }
  \item{backgroundfile}{ CSV format gene network used as the reference set, each line of which is Entrez IDs of two genes }
  \item{filter}{ indicates that HPO terms must have at least 'filter' edges annotated }
  \item{cutoff}{ significant cutoff for HPO enrichment analysis }
}
\value{
	Return a data.frame object with 8 columns.Details are below:
	
\item{"HPOID"}{enriched HPO IDs}
\item{"annEdgeNumber"}{Edge number annotated to this HPO term in the inputed gene network}
\item{"annBgNumber"}{Edge number in the inputed gene network}
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
  \code{\link{HPOGeneEnrichment}}
  \code{\link{HPOGeneNOAWholeNetEnrichment}}
}
\examples{
  edges <- data.frame(node1=c("65018","65019","65055","65055","65125","65266"),
                      node2=c("65055","65055","65125","65266","65266","79001"))
  graph<-graph.data.frame(edges,directed=FALSE)
  write.csv(edges,file="HPOGeneNOASubNetEnrichment-testnetwork.csv",row.names=FALSE) 
 
  edges <- data.frame(node1=c("65018","65018","65018","65018","65018","65055","65055","65055","65055","65109","65109","65109","65125","65125","65266"),
                      node2=c("65055","65109","65125","65266","79001","65019","65125","65266","79001","65125","65266","79001","65266","79001","79001"))
  graph<-graph.data.frame(edges,directed=FALSE)
  write.csv(edges,file="HPOGeneNOASubNetEnrichment-backgroundnetwork.csv",row.names=FALSE) 

  HPOGeneNOASubNetEnrichment("HPOGeneNOASubNetEnrichment-testnetwork.csv","HPOGeneNOASubNetEnrichment-backgroundnetwork.csv",cutoff=0.1)
}
\keyword{ manip }
