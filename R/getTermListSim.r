fmax<-function(a,anc,IC){
  fff<-function(x){
    if(!exists("Ancestors",envir=HPOSimEnv)) getAncestors()
    ancestor<-get("Ancestors",envir=HPOSimEnv)
    anca<-c(ancestor[x]$HP,x); anca<-anca[!is.na(anca)];
    anall<-intersect(anca,anc)
    if(length(anall)==0){
      return (0)
    }
    b<-as.double(IC[IC[,1] %in% anall,][,3])
    return (max(b))
  }
  res<-apply(a,1,fff)
  return (mean(res))
}
dp<-function(a,anc,IC){
  dfs<-function(t){
    if(!exists("Parents",envir=HPOSimEnv)) getParents()
    parent<-get("Parents",envir=HPOSimEnv)
    if(t %in% anc) return (as.double(IC[IC[,1] %in% t,][,3]))
    parn<-parent[t]$HP; parn<-parn[is.na(parn)];
    if(length(parn)==0) return (0)
    tem<-apply(parn,1,dfs)
    return (max(tem))
  }
  res<-apply(a,1,dfs)
  return (mean(res))
}

getTermListSim <- function(anno1, anno2, combinemethod="funSimMax", method="Resnik", IC, verbose=FALSE){	  
##test1
#       anno1<-c("HP:0000118", "HP:0000152", "HP:0000234", "HP:0000271")
#       anno2<-c("HP:0000284", "HP:0000478", "HP:0000479", "HP:0000488")
#       combinemethod="funSimMax"
#       method="Resnik"
#       IC<-get("termIC",envir=HPOSimEnv)
#       verbose=FALSE
  ##test2
#       anno2<-c("HP:0000118","HP:0000951","HP:0001000","HP:0001574","HP:0200015")
#       anno1<-c("HP:0000118","HP:0000951","HP:0000962","HP:0001072","HP:0001574","HP:0001597","HP:0001807","HP:0002164") 
#       combinemethod="funSimMax"
#       method="Resnik"
#       IC<-get("DiseasetermIC",envir=HPOSimEnv)
#       verbose=FALSE

  
  if(length(anno1)*length(anno2) == 0) {
    warning(paste("No HPO information for",anno1,anno2,". Similarity set to NaN."))  	
    return(NaN)
  }
  
  if(!exists("Parents",envir=HPOSimEnv)) getParents()
  parent<-get("Parents",envir=HPOSimEnv)
  if(!exists("Ancestors",envir=HPOSimEnv)) getAncestors()
  ancestor<-get("Ancestors",envir=HPOSimEnv)
  
  #rowScore=mean(max value of each row)
  #colScore=mean(max value of each column)
  #funSimMax=max(rowScore,colScore)
  #funSimAvg=mean(rowScore,colScore)
  
  if(method=="Resnik") { # special optimize for Resnik 
    
    if(combinemethod=="funSimMax"){
      if(length(anno1)==length(anno2) && all(anno1==anno2)){
        b<-as.double(IC[IC[,1] %in% anno1,][,3])
        return (mean(b))
      }
      
      dim(anno1)<-length(anno1); anc1<-unlist(apply(anno1,1,function(x) ancestor[x]$HP)); anc1<-unique(union(anc1,anno1)); anc1<-anc1[!is.na(anc1)];
      dim(anno2)<-length(anno2); anc2<-unlist(apply(anno2,1,function(x) ancestor[x]$HP)); anc2<-unique(union(anc2,anno2)); anc2<-anc2[!is.na(anc2)];
      
      colMax<-fmax(anno1,anc2,IC) 
      rowMax<-fmax(anno2,anc1,IC)
      return (max(colMax,rowMax))
    }
    # 
    if(combinemethod=="funSimAvg"){
      if(length(anno1)==length(anno2) && all(anno1==anno2)){
        b<-as.double(IC[IC[,1] %in% anno1,][,3])
        return (mean(b))
      }
      
      dim(anno1)<-length(anno1); anc1<-unlist(apply(anno1,1,function(x) ancestor[x]$HP)); anc1<-unique(union(anc1,anno1)); anc1<-anc1[!is.na(anc1)];
      dim(anno2)<-length(anno2); anc2<-unlist(apply(anno2,1,function(x) ancestor[x]$HP)); anc2<-unique(union(anc2,anno2)); anc2<-anc2[!is.na(anc2)];
      
      colMax<-fmax(anno1,anc2,IC)
      rowMax<-fmax(anno2,anc1,IC)
      return(0.5*(rowMax + colMax))
    }
    
  }
  ker<-matrix(0,nrow=length(anno1),ncol=length(anno2),dimnames=list(anno1,anno2))	 ##
  kerzero<-matrix(0,nrow=length(anno1),ncol=length(anno2),dimnames=list(anno1,anno2))   ##
  #     
  for(i in 1:length(anno1)){
    for(j in 1:length(anno2)) {
      if(kerzero[i,j]==0) {
        ker[i,j]<-calcTermSim(anno1[i],anno2[j], method, IC)  
        kerzero[i,j]<-1
        if( (anno1[i] %in% anno2) && (anno2[j] %in% anno1)){
          ker[anno2[j],anno1[i]]<-ker[i,j]
          kerzero[anno2[j],anno1[i]]<-1
        }
      }
    }
  }
  
  if(combinemethod == "max"){				
    return(max(ker))
  }
  else if(combinemethod == "mean"){				
    return(mean(ker))
  }  
  else if(combinemethod == "funSimAvg"){
    rowMax = mean(apply(ker,1,max))
    colMax = mean(apply(ker,2,max))
    return(0.5*(rowMax + colMax))
  }
  else if(combinemethod == "funSimMax"){
    rowMax = mean(apply(ker,1,max))
    colMax = mean(apply(ker,2,max))
    return(max(rowMax, colMax))
  }
  else if(combinemethod =="BMA"){
    m=nrow(ker)
    n=ncol(ker)
    return((sum(apply(ker,1,max))+sum(apply(ker,2,max)))/(m+n))
  }	
  else
    stop(paste("getTermListSim: Unknown gene combinemethod",combinemethod,"!"))
}