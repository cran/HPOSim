getLCA <-
  function(term1, term2, IC){  
    .initialize()
    if(!exists("Ancestors",envir=HPOSimEnv)) getAncestors()
    ancestor<-get("Ancestors",envir=HPOSimEnv)
    an1<-ancestor[[term1]]
    an2<-ancestor[[term2]]
    an1[length(an1)+1]<-term1 #add the term itself to the list of its ancestors
    an2[length(an2)+1]<-term2 #
    # find common ancestor with maximal information content
    anall<-intersect(an1, an2) #所有公共祖先
    #IC<-get("termIC", envir=HPOSimEnv)
    info<-IC[IC[,1] %in% anall,] # IC of all the common ancestors
    info<-info[order(info[,3],decreasing=T),] # IC in decreasing order
    lca<-as.character(info[1,1])
    return(lca)
  }