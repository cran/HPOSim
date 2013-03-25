\name{getTermListSim}
\alias{getTermListSim}
\title{ Get pairwise similarities of a list of HPO terms }
\description{
  Given a list of HPO terms, the function calculates the pairwise similarities between the terms.}
\usage{
  getTermListSim(termlist, method = "Resnik", verbose = FALSE)
}
\arguments{
  \item{termlist}{ some HPO terms }
  \item{method}{ one of "Resnik", "JiangConrath", "Lin", "CoutoResnik", "CoutoJiangConrath", "CoutoLin", "relevance", "GIC", "simIC" and "Wang" }
  \item{verbose}{ print some information }
}
\details{
  Currently the following methods for computing DO term similarities are implemented: 
  \describe{
    \item{"Resnik"}{information content of minimum subsumer (ICms) [1]}
    \item{"JiangConrath"}{\eqn{1-\min(1,IC(term1)-2ICms+IC(term2))}{1-min(1,IC(term1)-2ICms+IC(term2))}[2]}           
    \item{"Lin"}{\eqn{\frac{2ICms}{(IC(term1) + IC(term2))}}{2ICms/(IC(term1)+IC(term2))} [3]}
    \item{"CoutoResnik"}{average information content of common disjunctive ancestors of \code{term1} and \code{term2} (ICshare) [4]}
    \item{"CoutoJiangConrath"}{\eqn{1 - \min(1, IC(term1) -2ICshare + IC(term2))}{1-min(1,IC(term1)-2ICshare+IC(term2))} [4]}
    \item{"CoutoLin"}{\eqn{\frac{2ICshare}{(IC(term1) + IC(term2))}}{2ICshare/(IC(term1) + IC(term2))} [4]}
    \item{"relevance"}{sim_Lin * (1 - exp(-ICms)) [5]}
    \item{"GIC"}{summed information content of common ancestors divided by summed information content of all ancestors of term1 and term2 [7]}
    \item{"simIC"}{sim_Lin * (1 - 1/(1 + ICms)) [7]}
    \item{"Wang"}{\eqn{Sim(term1,term2)=\frac{\sum_{t\in T_{term1}\cap T_{term2}} (S_{term1}(t)+S_{term2}(t))}{SV(term1)+SV(term2)}} [8]}
  }
}
\value{
  n * n matrix with similarities between HPO terms (n = number of HPO terms in the list) . 
}

\references{
[1] P. Resnik, Using Information Content to evaluate semantic similarity in a taxonomy, Proc. 14th Int. Conf. Artificial Intel., 1995 

[2] J. Jiang, D. Conrath, Semantic Similarity based on Corpus Statistics and Lexical Taxonomy, Proc. Int. Conf. Research in Comp. Ling., 1998 

[3] D. Lin, An Information-Theoretic Definition of Similarity, Proc. 15th Int. Conf. Machine Learning, 1998 

[4] Couto, F.; Silva, M. & Coutinho, P., Semantic Similarity over the Gene Ontology: Family Correlation and Selecting Disjunctive Ancestors, Conference in Information and Knowledge Management, 2005  

[5] A. Schlicker, F. Domingues, J. Rahnenfuehrer, T. Lengauer, A new measure for functional similarity of gene products based on Gene Ontology, BMC Bioinformatics, 7, 302, 2006.

[6] C. Pesquita, D. Faria, H. Bastos, A. Falcao, F. Couto, Evaluating GO-based Semantic Similarity Measures, In: Proc. 10th Annual Bio-Ontologies Meeting 2007, 37 - 40, 2007

[7] B. Li, J. Wang, A. Feltus, J. Zhou, F. Luo, Effectively Integrating Information Content and Structural Relationship to Improve the GO-based Similarity Measure Between Proteins, BMC Bioinformatics, 2009. 

[8] James Z.Wang,Zhidian Du, et al. A new method to measure the semantic similarity of GO terms.Bioinformatics 2007,Vol 23,1274-1281.
}

\author{
  Yue Deng<anfdeng@163.com>
}
\seealso{ 
  \code{\link{getMinimumSubsumer}}
  \code{\link{calcTermSim}}
}

\examples{
  getTermListSim(c("HP:0000028","HP:0000033","HP:0000060"))
}
\keyword{ manip }


