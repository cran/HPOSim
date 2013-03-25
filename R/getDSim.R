getDSim <-
  function(anno1, anno2, combinemethod="max", method="Resnik", verbose=FALSE){    
    if(length(anno1) < length(anno2)){
      a1<-anno1
      a2<-anno2
      swap<-FALSE
    }
    else{
      a1<-anno2
      a2<-anno1
      swap<-TRUE
    }
    ker<-matrix(0,nrow=length(a1),ncol=length(a2))	
    for(i in 1:length(a1)){
      for(j in 1:length(a2))
        ker[i,j]<-calcDiseaseTermSim(a1[i],a2[j], method, verbose)	
    }
    
    
    if(length(a1)*length(a2) > 0){
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
        stop(paste("getGSim: Unknown gene combinemethod",combinemethod,"!"))
    }
    else{	
      warning(paste("No HPO information for",a1,a2,". Similarity set to NaN."))		
      return(NaN)
    }
  }
