getGeneListSim<-
function(genelist,combinemethod="funSimMax",method="Resnik",ontology="PA",normalization=FALSE,normalizationmethod="Lin",verbose=FALSE){
	.initialize()
	if(length(genelist)<2)
  		stop("Genelist should contain at least two genes!")
	ker<-matrix(0,nrow=length(genelist),ncol=length(genelist))	
	for(i in 1:length(genelist))	
		for(j in 1:i){
			ker[i,j]<-getGeneSim(genelist[i],genelist[j],combinemethod,method,ontology,normalization,normalizationmethod,verbose)
			ker[j,i]<-ker[i,j]
		}
	ker
}
getGeneSim <-
function(gene1,gene2, combinemethod="funSimMax", method="Resnik", ontology="PA", normalization=FALSE, normalizationmethod="Lin", verbose=FALSE){
#   gene1<-"410"
#   gene2<-"650"
#   combinemethod="funSimMax"
#   method="Resnik"
#   ontology="PA"
#   normalization=FALSE
#   normalizationmethod="Lin"
#   verbose=F
  if(verbose)
       message(paste("Calculating Gene Sim: ", gene1, "and", gene2))
  .initialize()
	if(length(gene1) > 1||length(gene2)>1)
		stop("Gene1 and Gene2  should contain only one element!")
	IC<-get("termIC",envir=HPOSimEnv)
	gene2hpo<-get("gene2hpo",envir=HPOSimEnv)

	Terms1<-as.character(unlist(gene2hpo[gene1]))
	Terms2<-as.character(unlist(gene2hpo[gene2]))
  if(length(Terms1)==0) {		
    warning(paste(gene1,"has No effective annotation!"))
    return(0)
  }
  if(length(Terms2)==0) {  	
    warning(paste(gene2,"has No effective annotation!"))
    return(0)
  }

  Terms1<-RemoveTermsWithoutIC(Terms1,ontology,IC)
  Terms2<-RemoveTermsWithoutIC(Terms2,ontology,IC)
  

	if(length(Terms1) >= 1 && length(Terms2)>=1){		
		Ker<-getTermListSim(Terms1, Terms2, combinemethod, method, IC, verbose)

		if(normalization){		
		  Ker1<-getTermListSim(Terms1, Terms1, combinemethod, method, IC, verbose)
		  Ker2<-getTermListSim(Terms2, Terms2, combinemethod, method, IC, verbose)
			Ker = normalize.kernel(Ker,normalizationmethod,Ker1,Ker2)
		}			
	}
	else{
		if(length(Terms1) == 0)		
		  warning(paste(gene1,"has No effective annotation!"))
		if(length(Terms2) == 0)
		  warning(paste(gene2,"has No effective annotation!"))
		Ker<-0					
	}
	return(Ker)
}
