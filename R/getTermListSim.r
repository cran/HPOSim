getTermListSim <-
function(termlist, method="Resnik", verbose=FALSE){
	S<-matrix(0,nrow=length(termlist),ncol=length(termlist))
	colnames(S)<-termlist
	rownames(S)<-termlist
	for(i in 1:length(termlist)){
		S[i,i] <- calcTermSim(termlist[i],termlist[i], method, verbose)				
		if(i > 1){
			for(j in 1:(i-1)){				
				S[i,j]<- calcTermSim(termlist[i],termlist[j], method, verbose)				
				S[j,i]<-S[i,j]
			}
		}
	}
	S
}

