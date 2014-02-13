\name{HPOGeneNOAWholeNetEnrichment}
\alias{HPOGeneNOAWholeNetEnrichment}
\title{HPO enrichment analysis of gene network using NOA}
\description{
	Given a gene interaction network, the function returns the enriched HPO terms using whole-net method of NOA.
}
\usage{
	HPOGeneNOAWholeNetEnrichment(file,filter=5,cutoff=0.05)
}

\arguments{
  \item{file}{ CSV format gene interaction network, each line of which is Entrez IDs of two interacting genes/proteins }
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
}
\examples{
 
  edges <- data.frame(node1=c("581", "581", "581", "1499", "1499", "1630", 
                              "60", "257", "999", "1282", "1723", "1760"),
                      node2=c("1499", "1630", "2033", "1630", "2033", 
                              "2033", "581", "1499", "2033", "1630", "1760", "581"))
  graph<-igraph::graph.data.frame(edges,directed=FALSE)
  write.csv(edges,file="HPOGeneNOAWholeNetEnrichment-testwork.csv",row.names=FALSE) 
  HPOGeneNOAWholeNetEnrichment("HPOGeneNOAWholeNetEnrichment-testwork.csv",cutoff=0.1)
}
\keyword{ manip }
