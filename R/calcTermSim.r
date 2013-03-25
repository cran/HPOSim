calcTermSim <-
  function(term1, term2, method="Resnik", verbose=FALSE){
    initialize()
    IC<-get("termIC", envir=HPOSimEnv)
    if(method== "Resnik"){
      an=getMinimumSubsumer(term1,term2)
      if(an == "NA")
        return(0)
      else {
        return(as.double(IC[IC[1]==an][3]))
      }  
    }
    else if(method == "JiangConrath"){
      an=getMinimumSubsumer(term1,term2)
      if(an == "NA"||an=="HP:0000001"||an=="HP:0000004"||an=="HP:0000005"||an=="HP:0000118")
        return(0)
      if(term1==term2)
        return(1)  
      else{
        res= - 1/ (  1 + 2*as.double(IC[IC[1]==an][3]) - as.double(IC[IC[1]==term1][3]) - as.double(IC[IC[1]==term2][3]) )     
        return(res)
      }
    }
    
    else if(method == "Lin"){
      an=getMinimumSubsumer(term1,term2)
      if(an == "NA"||an=="HP:0000001"||an=="HP:0000004"||an=="HP:0000005"||an=="HP:0000118")
        return(0)
      else{
        res = 2*as.double(IC[IC[1]==an][3])/(as.double(IC[IC[1]==term1][3])+as.double(IC[IC[1]==term2][3]))
        return(ifelse(is.na(res), 1, res)) 
      }  
    }
    
    else if(method == "simIC"){ # Li et al.
      an = getMinimumSubsumer(term1,term2)
      if(an == "NA"||an=="HP:0000001"||an=="HP:0000004"||an=="HP:0000005"||an=="HP:0000118")
        return(0)
      else{
        res = 2*as.double(IC[IC[1]==an][3])/(as.double(IC[IC[1]==term1][3])+as.double(IC[IC[1]==term2][3])) * (1 - 1/(1 + as.double(IC[IC[1]==an][3])))
        return(ifelse(is.na(res), 1, res))
      }
    }
    else if(method == "relevance"){ # Schlicker et al.
      an = getMinimumSubsumer(term1,term2)
      if(an == "NA"||an=="HP:0000001"||an=="HP:0000004"||an=="HP:0000005"||an=="HP:0000118"){
        return(0)}
      
      else{
        res = (2*as.double(IC[IC[1]==an][3])/(as.double(IC[IC[1]==term1][3])+as.double(IC[IC[1]==term2][3])))*(1 - exp(-as.double(IC[IC[1]==an][3])))
        return(ifelse(is.na(res), 1, res))
      }  
    }
    
    else if(method == "GIC") # graph information content
      return(getGIC(term1, term2))
    
    else if(method == "Wang"){
      res=getSimWang(term1,term2)
      return(res)
    }
    else
      stop(paste("calcTermSim: Unknown term similarity",method))
  }


getSimWang<-function(term1,term2){
  initialize()
  
  if(term1 == term2){
    return(1);
  }
  if(!exists("Ancestors",envir=HPOSimEnv)) getAncestors()
  ancestor<-get("Ancestors",envir=HPOSimEnv)
  we=0.7
  
  an1<-ancestor[names(ancestor)==term1]$HP
  an2<-ancestor[names(ancestor)==term2]$HP
  
  an1<-c(an1,term1)
  an2<-c(an2,term2)
  
  an1<-unique(an1)
  an2<-unique(an2)
  common<-intersect(an1,an2)
  if(length(common)==0){
    return(0);
  }
  
  SA<-list()
  SB<-list()
  SA[term1]=1;
  SB[term2]=1;
  
  done <- FALSE
  while (!done) {
    if(all(an1 %in% names(SA))){
      done=TRUE;
    }else{
      parents<-unique(unlist(getTermParents(names(SA),verbose=FALSE)))
      parents<-parents[!is.na(parents)]
      v<-parents[! parents %in% names(SA)]
      if(length(v) == 0){
        done=TRUE;
      }else{
        vv<-v[1]
        vvchildren<-unique(unlist(getTermChildren(vv,verbose=FALSE)))
        s=max(sapply(SA[vvchildren[vvchildren %in% names(SA)]],function(x){we*x}))
        SA[vv]<-s	
      }
    }
  }
  
  done<-FALSE
  while (!done) {
    if(all(an2 %in% names(SB))){
      done=TRUE;
    }else{
      parents<-unique(unlist(getTermParents(names(SB),verbose=FALSE)))
      parents<-parents[!is.na(parents)]
      v<-parents[! parents %in% names(SB)]
      if(length(v) == 0){
        done=TRUE;
      }else{
        vv<-v[1]
        vvchildren<-unique(unlist(getTermChildren(vv,verbose=FALSE)))
        s=max(sapply(SB[vvchildren[vvchildren %in% names(SB)]],function(x){we*x}))
        SB[vv]<-s	
      }			
    }
  }
  
  SA["HP:0000001"]<-0
  SB["HP:0000001"]<-0
  
  SVA= sum(unlist(SA))
  SVB= sum(unlist(SB))
  res<-(sum(unlist(SA[common]))+sum(unlist(SB[common])))/(SVA+SVB)
  
  res	
  
}


getGIC <-
  function(term1, term2){
    initialize()	
    if(term1 == term2){
      return(1)
    }
    IC<-get("termIC", envir=HPOSimEnv)
    ancestor<-get("Ancestors",envir=HPOSimEnv)
    an1<-ancestor[names(ancestor) == term1]$HP
    an2<-ancestor[names(ancestor) == term2]$HP
    
    ancommon <- intersect(an1, an2)
    anunion <- union(an1, an2)
    a<-0
    for(i in 1:length(ancommon))
    {
      a <- a+as.double(IC[IC[1]==ancommon[i]][3])
    }		
    
    b<-0
    for(i in 1:length(anunion))
    {
      b <- b+as.double(IC[IC[1]==anunion[i]][3])
    }		
    
    res=a/b
    return(res)
  }

