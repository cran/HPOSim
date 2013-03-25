initialize<-function(pos = 1,envir = as.environment(pos)){
  if(!exists("HPOSimEnv") || length(HPOSimEnv)<1) {
    print("initializing HPOSim package ...")		
    #assign("HPOSimEnv",new.env(),envir=.GlobalEnv)  
    assign("HPOSimEnv",new.env(),envir=envir)  
    #HPOSimEnv<-new.env()

    data("GeneTermIC",envir=HPOSimEnv)	
    data("DiseaseTermIC",envir=HPOSimEnv)	
    data("GeneToPhenotype",envir=HPOSimEnv)
    data("PhenotypeToGene",envir=HPOSimEnv)
    data("DiseaseToPhenotype",envir=HPOSimEnv)
    data("PhenotypeToDisease",envir=HPOSimEnv)
    print("finished.")
  }
}

getOffsprings<-function(){
  require(HPO.db)
  initialize()
  res<-AnnotationDbi::as.list(HPOOFFSPRING)
  assign("Offsprings",res,envir=HPOSimEnv)
}

getAncestors<-function(){
  require(HPO.db)
  initialize()
  res<-AnnotationDbi::as.list(HPOANCESTOR)
  assign("Ancestors",res,envir=HPOSimEnv)
}

getParents<-function(){
  require(HPO.db)
  initialize()
  res<-AnnotationDbi::as.list(HPOPARENTS)
  assign("Parents",res,envir=HPOSimEnv)
}

getChildren<-function(){
  require(HPO.db)
  initialize()
  res<-AnnotationDbi::as.list(HPOCHILDREN)
  assign("Children",res,envir=HPOSimEnv)
}

getTermOffsprings<-
  function(hpolist,verbose=FALSE){
    if(is.list(hpolist)){
      hpolist<-unique(unlist(hpolist))
    }else{
      hpolist<-unique(hpolist)
    }
    if(verbose){
      print("Start to fetch the offsprings")
    }
    if(!exists("HPOSimEnv")||!exists("Offsprings",envir=HPOSimEnv)) getOffsprings()
    offspring<-get("Offsprings",envir=HPOSimEnv)
    res<-offspring[hpolist[hpolist %in% names(offspring)]]
    notmatch<-hpolist[! hpolist %in% names(offspring)]
    if(length(notmatch)>0){
      for(i in 1:length(notmatch)){
        res[[notmatch[i]]]<-NA
      }
    }
    res
  }

getTermParents<-
  function(hpolist,verbose=FALSE){
    if(is.list(hpolist)){
      hpolist<-unique(unlist(hpolist))
    }else{
      hpolist<-unique(hpolist)
    }
    if(verbose){
      print("Start to fetch the parents")
    }
    if(!exists("HPOSimEnv")||!exists("Parents",envir=HPOSimEnv)) getParents()
    parent<-get("Parents",envir=HPOSimEnv)
    res<-parent[hpolist[hpolist %in% names(parent)]]
    notmatch<-hpolist[! hpolist %in% names(parent)]
    if(length(notmatch)>0){
      for(i in 1:length(notmatch)){
        res[[notmatch[i]]]<-NA
      }
    }
    res
  }

getTermAncestors<-
  function(hpolist,verbose=FALSE){
    if(is.list(hpolist)){
      hpolist<-unique(unlist(hpolist))
    }else{
      hpolist<-unique(hpolist)
    }
    if(verbose){
      print("Start to fetch the ancestors")
    }
    if(!exists("HPOSimEnv")||!exists("Ancestors",envir=HPOSimEnv)) getAncestors()
    ancestor<-get("Ancestors",envir=HPOSimEnv)
    res<-ancestor[hpolist[hpolist %in% names(ancestor)]]
    notmatch<-hpolist[! hpolist %in% names(ancestor)]
    if(length(notmatch)>0){
      for(i in 1:length(notmatch)){
        res[[notmatch[i]]]<-NA
      }
    }
    res
  }

getTermChildren<-
  function(hpolist,verbose=FALSE){
    if(is.list(hpolist)){
      hpolist<-unique(unlist(hpolist))
    }else{
      hpolist<-unique(hpolist)
    }
    if(verbose){
      print("Start to fetch the children")
    }
    if(!exists("HPOSimEnv")||!exists("Children",envir=HPOSimEnv)) getChildren()
    child<-get("Children",envir=HPOSimEnv)
    res<-child[hpolist[hpolist %in% names(child)]]
    notmatch<-hpolist[! hpolist %in% names(child)]
    if(length(notmatch)>0){
      for(i in 1:length(notmatch)){
        res[[notmatch[i]]]<-NA
      }
    }
    res
  }

getTerms<-function(){
  require(HPO.db)
  initialize()
  res<-AnnotationDbi::as.list(HPOTERM)
  assign("HPOTerms",res,envir=HPOSimEnv)
}


getTerm<-
  function(hpolist){
    if(is.list(hpolist)){
      hpolist<-unique(unlist(hpolist))
    }else{
      hpolist<-unique(hpolist)
    }
    if(!exists("HPOSimEnv")||!exists("HPOTerm",envir=HPOSimEnv)) getTerms()
    term<-get("HPOTerms",envir=HPOSimEnv)
    res<-term[hpolist[hpolist %in% names(term)]]
    notmatch<-hpolist[! hpolist %in% names(term)]
    if(length(notmatch)>0){
      warning(paste("===>",length(notmatch),"of",length(hpolist),"HPIDs not mapped to current disease ontology\n"))
    }
    res
  }


calculateICS<-function(){
  require(HPO.db)
  initialize()
  PhenotypeToGenes<-read.csv(file="PhenotypeToGenes.csv",header=FALSE,sep=",")
  nrows<-length(PhenotypeToGenes[,1])
  if(!exists("HPOSimEnv")||!exists("HPOTerm",envir=HPOSimEnv)) 
    getTerms()
  term<-get("HPOTerms",envir=HPOSimEnv)
  termlist<-as.factor(names(term))
  termfrequency<-data.frame(termlist,0)
  
  ##frequencies of special terms£¬root nodes of 3 ontologies
  termfrequency[termfrequency[1]=="HP:0000001",][2]<-1826  ##all
  termfrequency[termfrequency[1]=="HP:0000004",][2]<-669   ##root nodes of 3 ontologies
  termfrequency[termfrequency[1]=="HP:0000005",][2]<-1764
  termfrequency[termfrequency[1]=="HP:0000118",][2]<-1812
  
  ##calculate frequencies of all terms
  i<-1
  for(i in 5400:nrows)
  {
    xx<-PhenotypeToGenes[i,]
    xx<-xx[!is.na(xx)]
    ##termfrequency[i,2]<-length(xx)-1
    termfrequency[termfrequency[1]==xx[1],][2]<-length(xx)-1
  }
  
  OC<-getTermOffsprings("HP:0000004",verbose=FALSE)
  OC<-OC$HP
  MI<-getTermOffsprings("HP:0000005",verbose=FALSE)
  MI<-MI$HP
  PA<-getTermOffsprings("HP:0000118",verbose=FALSE)
  PA<-PA$HP
  
  ##calculate IC
  termIC<-data.frame(termlist,ontology="",0) 
  termIC$ontology<-as.character(termIC$ontology)
  termIC[1,]$ontology="All"
  termIC[termIC[1]=="HP:0000004",]$ontology<-"OC"
  termIC[termIC[1]=="HP:0000005",]$ontology<-"MI"
  termIC[termIC[1]=="HP:0000118",]$ontology<-"PA"
  for(i in 2:length(termlist))
  {
    xx<-as.character(termIC[i,]$termlist)
    if(xx %in% PA)
    {termIC[termIC[1]==xx,]$ontology<-"PA"
     termIC[termIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/1812)
    }
    if(xx %in% OC)
    {termIC[termIC[1]==xx,]$ontology<-"OC"
     termIC[termIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/669)
    }
    if(xx %in% MI)
    {termIC[termIC[1]==xx,]$ontology<-"MI"
     termIC[termIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/1764)
    }	
  }
  save(termIC,file="termIC.rda")
}

