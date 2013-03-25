getMinimumSubsumer <-
function(term1, term2){	
	initialize()
	if(!exists("Ancestors",envir=HPOSimEnv)) getAncestors()
	
	ancestor<-get("Ancestors",envir=HPOSimEnv)
	if(term1 == term2){
		ms<-term1	
		return(ms)
	}
	an1<-ancestor[names(ancestor) == term1]$HP
	an2<-ancestor[names(ancestor) == term2]$HP
	case1<-which(an2 == term1)  # term1 is the ms of term2
	case2<-which(an1 == term2) # term2 is the ms of term1	  
	if(length(case1) > 0){
		ms<-term1	
	} else if(length(case2) > 0) {
		ms<-term2	
	} else {
		# find common ancestor with maximal information content
		anall<-intersect(an1, an2) 
		IC<-get("termIC", envir=HPOSimEnv)
		a<-0
		for(i in 1:length(anall))
		{
			a[i]<-IC[IC[1]==anall[i]][3]
		}		
		ms<-anall[as.double(which.max(a))]
	}	
	if(is.null(ms) | length(ms) == 0)
		ms <- "NA"	
	ms
}

