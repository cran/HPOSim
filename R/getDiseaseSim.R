getDiseaseListSim<-
  function(diseaselist,combinemethod="funSimMax",method="Resnik",ontology="PA",normalization=FALSE,normalizationmethod="Lin",verbose=FALSE){
    initialize()
    if(length(diseaselist)<2)
      stop("Diseaselist should contain at least two diseases!")
    ker<-matrix(0,nrow=length(diseaselist),ncol=length(diseaselist))	
    for(i in 1:length(diseaselist))	
      for(j in 1:i){
        ker[i,j]<-getDiseaseSim(diseaselist[i],diseaselist[j],combinemethod,method,ontology,normalization,normalizationmethod,verbose)
        ker[j,i]<-ker[i,j]
      }
    ker
  }

getDiseaseSim <-
  function(disease1,disease2, combinemethod="funSimMax", method="Resnik", ontology="PA", normalization=FALSE, normalizationmethod="Lin", verbose=FALSE){
    if(verbose)
      print(paste("Calculating similarity between disease", disease1, "and", disease2, "with term similarity measure",method, "and combined with", combinemethod))
    initialize()
    if(length(disease1) > 1||length(disease2)>1)
      stop("Disease1 and Disease2  should contain only one element!")

    IC<-get("DiseasetermIC",envir=HPOSimEnv)
    disease2hpo<-get("disease2hpo",envir=HPOSimEnv)

    Terms1<-disease2hpo[disease1]
    Terms2<-disease2hpo[disease2]

    if(length(Terms1[[1]])==0) {  	
      warning(paste(disease1,"has No effective annotation!"))
      return(0)
    }
    if(length(Terms2[[1]])==0) {  	
      warning(paste(disease2,"has No effective annotation!"))
      return(0)
    }
    for(i in 1:length(Terms1[[1]])) 
    {
      if(is.na(IC[IC[1]==Terms1[[1]][i]][1]) || IC[IC[1]==Terms1[[1]][i]][2]!=ontology || IC[IC[1]==Terms1[[1]][i]][3]==0)
        Terms1[[1]][i]<-NA
    }
    for(i in 1:length(Terms2[[1]]))  
    {
      if(is.na(IC[IC[1]==Terms2[[1]][i]][1]) || IC[IC[1]==Terms2[[1]][i]][2]!=ontology || IC[IC[1]==Terms2[[1]][i]][3]==0)
        Terms2[[1]][i]<-NA
    }
    Terms1<-na.omit(unlist(Terms1))
    Terms2<-na.omit(unlist(Terms2))
    
    if(length(Terms1) >= 1 && length(Terms2)>=1){		
      Ker<-getDSim(Terms1,Terms2,combinemethod,method,verbose)
      
      if(normalization){		
        Ker1<-getDSim(Terms1,Terms1,combinemethod,method,verbose)
        Ker2<-getDSim(Terms2,Terms2,combinemethod,method,verbose)
        Ker = normalize.kernel(Ker,normalizationmethod,Ker1,Ker2)
      }			
    }
    else{
      if(length(Terms1) == 0)
        warning(paste(disease1,"has No effective annotation!"))					
      if(length(Terms2) == 0)
        warning(paste(disease2,"has No effective annotation!"))
      Ker<-0					
    }
    return(Ker)
  }
