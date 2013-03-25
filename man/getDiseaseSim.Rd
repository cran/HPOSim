\name{getDiseaseSim}
\alias{getDiseaseSim}
\title{ Phenotype-based similarity for two diseases base on HPO }
\description{
  Given two diseases, the function calculates similarity of the diseases based on HPO.
}
\usage{
  getDiseaseSim(disease1,disease2,combinemethod="funSimMax",method="Resnik",ontology="PA",normalization=FALSE,normalizationmethod="Lin", verbose=FALSE)
}
\arguments{
  \item{disease1}{ one disease OMIM ID }
  \item{disease2}{ another disease OMIM ID }
  \item{combinemethod}{ method to calculate the similarity between diseases based on , one of "max", "mean", "funSimMax", "funSimAvg" and ""BMA }
  \item{method}{ method to compute the similarity of HPO terms, one of "Resnik", "JiangConrath", "Lin", "simIC", "relevance" and "Wang" }
  \item{ontology}{ the ontology used to calculate the similarity, one of "PA"(Phenotypic Abnormality), "MI"(Mode of Inheritance) and "OC"(Onset and Clinical Course)}
  \item{normalization}{ whether do the normalization }
  \item{normalizationmethod}{ one of "sqrt", "Lin" and "Tanimoto"}
  \item{verbose}{ print out some information }
}
\details{
  The method to combine the similarities between HPO terms can either be: 
  \describe{
    \item{"max"}{the maximum similarity between any two HPO terms}
    \item{"mean"}{the average similarity between any two HPO terms1}
    \item{"funSimMax"}{ the average of best matching HPO term similarities. Take the maximum of the scores achieved by assignments of HPO terms from disease 1 to disease 2 and vice versa.}
    \item{"funSimAvg"}{ the average of best matching HPO term similarities. Take the average of the scores achieved by assignments of HPO terms from disease 1 to disease 2 and vice versa.}
    \item{"BMA"}{best match average approach} 
  }
}
\value{
  the value of the similarity of two diseases
}
\author{
  Yue Deng<anfdeng@163.com>
}
\references{ 
[1] P. W. Lord, et al., "Investigating semantic similarity measures across the Gene Ontology: the relationship between sequence and annotation," Bioinformatics, vol. 19, pp. 1275-83, Jul 1 2003.

[2] A. Schlicker, F. Domingues, J. Rahnenfuehrer, T. Lengauer, A new measure for functional similarity of gene products based on Gene Ontology, BMC Bioinformatics, 7, 302, 2006.

[3] James Z.Wang,Zhidian Du, et al. A new method to measure the semantic similarity of GO terms.Bioinformatics 2007,Vol 23,1274-1281.

}
\seealso{  
  \code{\link{getDiseaseListSim}}
}
\examples{
  getDiseaseSim("OMIM:101900","OMIM:102000")
}
\keyword{ manip }
