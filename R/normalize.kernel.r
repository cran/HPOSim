normalize.kernel <-
function(Ker, method="none",Ker1,Ker2){
	if(method != "none"){
		if(method == "sqrt"){ # result between -1 and 1
			Kd<-sqrt(Ker1*Ker2)
			Ker<-Ker/Kd			
		}
		else if(method == "Lin"){ # result: diagonal = 1		
			Ker = 2*Ker / (Ker1+Ker2)
		}
		else if(method == "Tanimoto"){ 
			Ker = Ker / (Ker1 + Ker2 - Ker)
		}								
		else
			stop(paste("Unknown normalization method", method))
	}
	Ker
}
