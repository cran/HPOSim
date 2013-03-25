calcDiseaseTermSim <-
  function(term1, term2, method="Resnik", verbose=FALSE){
    initialize()
    IC<-get("DiseasetermIC", envir=HPOSimEnv)
    if(verbose)
      print(paste("Terms:",term1,",",term2,"( method:",method,")"))  
    
    if(method== "Resnik"){
      an=getMinimumSubsumer(term1,term2)
      if(an == "NA")
        return(0)
      else
        return(as.double(IC[IC[1]==an][3]))
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
      return(getDIC(term1, term2))
    
    else if(method == "Wang"){
      res=getSimWang(term1,term2)
      return(res)
    }
    else
      stop(paste("calcTermSim: Unknown term similarity",method))
  }

