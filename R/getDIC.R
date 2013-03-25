getDIC <-
  function(term1, term2){
    initialize()	
    if(term1 == term2){
      return(1)
    }
    IC<-get("DiseasetermIC", envir=HPOSimEnv)
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

