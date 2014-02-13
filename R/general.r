.initialize<-function(pos = 1,envir = as.environment(pos)){
  if(!exists("HPOSimEnv") || length(HPOSimEnv)<1) {
    message("initializing HPOSim package ...")		
    assign("HPOSimEnv",new.env(),envir=envir)  

    data("GeneTermIC",envir=HPOSimEnv)
    data("DiseaseTermIC",package="HPOSim", envir=HPOSimEnv)	
    data("GeneToPhenotype",envir=HPOSimEnv)
    data("PhenotypeToGene",envir=HPOSimEnv)
    data("DiseaseToPhenotype",envir=HPOSimEnv)
    data("PhenotypeToDisease",envir=HPOSimEnv)
    #tryCatch(utils::data(list=c("GeneTermIC","DiseaseTermIC","GeneToPhenotype","PhenotypeToGene","DiseaseToPhenotype","PhenotypeToDisease"), envir=HPOSimEnv))

    message("done.")
  }
}

.initial<-function(pos = 1,envir = as.environment(pos)){
  if(!exists("HPOSimEnv") || length(HPOSimEnv)<1) {
    packageStartupMessage("initializing HPOSim package ...")		
    assign("HPOSimEnv",new.env(),envir=envir)  
    #tryCatch(utils::data(list=c("GeneTermIC","DiseaseTermIC","GeneToPhenotype","PhenotypeToGene","DiseaseToPhenotype","PhenotypeToDisease"), envir=HPOSimEnv))

    packageStartupMessage("done.")
  }
}

getOffsprings<-function(){
  #require(HPO.db)
  .initialize()
  res<-AnnotationDbi::as.list(HPOOFFSPRING)
  assign("Offsprings",res,envir=HPOSimEnv)
}

getAncestors<-function(){
  #require(HPO.db)
  .initialize()
  res<-AnnotationDbi::as.list(HPOANCESTOR)
  assign("Ancestors",res,envir=HPOSimEnv)
}

getParents<-function(){
  #require(HPO.db)
  .initialize()
  res<-AnnotationDbi::as.list(HPOPARENTS)
  assign("Parents",res,envir=HPOSimEnv)
}

getChildren<-function(){
  #require(HPO.db)
  .initialize()
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
      message("Start to fetch the offsprings")
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
      message("Start to fetch the parents")
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
      message("Start to fetch the ancestors")
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
      message("Start to fetch the children")
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
  #require(HPO.db)
  .initialize()
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
  
getTermOffspringDiseases<-function(t){
    hpo2disease<-get("hpo2disease",envir=HPOSimEnv)
    diseases<-unique(unlist(lapply(getTermOffsprings(t),function(x){hpo2disease[x]})))
    return(diseases)
}
  
calculateGeneTermIC<-function(){
  #require(HPO.db)
  .initialize()
  getTerms()
  term<-get("HPOTerms",envir=HPOSimEnv)
  termlist<-as.factor(names(term))
  hpo2gene<-get("hpo2gene",envir=HPOSimEnv)
  
  OCTerms<-unlist(getTermOffsprings("HP:0000004",verbose=FALSE))
  OCgenenumber<-length(hpo2gene["HP:0000004"]$HP)
  MITerms<-unlist(getTermOffsprings("HP:0000005",verbose=FALSE))
  MIgenenumber<-length(hpo2gene["HP:0000005"]$HP)
  PATerms<-unlist(getTermOffsprings("HP:0000118",verbose=FALSE))
  PAgenenumber<-length(hpo2gene["HP:0000118"]$HP)

  termfreq<-lapply(termlist,function(x){length(hpo2gene[as.character(x)]$HP)})
  termfrequency<-data.frame(termlist,unlist(termfreq))

  ##calculate IC
  termIC<-data.frame(termlist,ontology="",0) 
  termIC$ontology<-as.character(termIC$ontology)
  termIC[1,]$ontology="All"
  termIC[termIC[1]=="HP:0000004",]$ontology<-"OC"
  termIC[termIC[1]=="HP:0000005",]$ontology<-"MI"
  termIC[termIC[1]=="HP:0000118",]$ontology<-"PA"
  termsnumber<-length(termlist)
  for(i in 2:termsnumber)
  {
    message(paste(i,"of",termsnumber))
    xx<-as.character(termIC[i,]$termlist)
    if(xx %in% PATerms){
      termIC[termIC[1]==xx,]$ontology<-"PA"
      termIC[termIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/PAgenenumber)
    }
    if(xx %in% OCTerms){
      termIC[termIC[1]==xx,]$ontology<-"OC"
      termIC[termIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/OCgenenumber)
    }
    if(xx %in% MITerms){
      termIC[termIC[1]==xx,]$ontology<-"MI"
      termIC[termIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/MIgenenumber)
    }
  }
  save(termIC,file="GeneTermIC.rda",compress="xz")
}

calculateDiseaseTermIC<-function(){
  #require(HPO.db)
  .initialize()
  getTerms()
  term<-get("HPOTerms",envir=HPOSimEnv)
  termlist<-as.factor(names(term))
  hpo2disease<-get("hpo2disease",envir=HPOSimEnv)
  
  OCTerms<-unlist(getTermOffsprings("HP:0000004"))
  OCdiseasenumber<-length(hpo2disease[["HP:0000004"]])
  MITerms<-unlist(getTermOffsprings("HP:0000005"))
  MIdiseasenumber<-length(hpo2disease[["HP:0000005"]])
  PATerms<-unlist(getTermOffsprings("HP:0000118"))
  PAdiseasenumber<-length(hpo2disease[["HP:0000118"]])
  
  #OCdiseasenumber<-length(c(getTermOffspringDiseases("HP:0000004"),hpo2disease["HP:0000004"]))
  #MIdiseasenumber<-length(c(getTermOffspringDiseases("HP:0000005"),hpo2disease["HP:0000005"]))
  #PAdiseasenumber<-length(c(getTermOffspringDiseases("HP:0000118"),hpo2disease["HP:0000118"]))
  
  #termfreq<-lapply(termlist,function(x){length(unique(c(getTermOffspringDiseases(x),hpo2disease[x])))})

  termfreq<-lapply(termlist,function(x){length(hpo2disease[as.character(x)]$"HP")})
  termfrequency<-data.frame(termlist,unlist(termfreq))
  
  ##calculate Disease Term IC
  DiseasetermIC<-data.frame(termlist,ontology="",0) 
  DiseasetermIC$ontology<-as.character(DiseasetermIC$ontology)
  DiseasetermIC[1,]$ontology="All"
  DiseasetermIC[DiseasetermIC[1]=="HP:0000004",]$ontology<-"OC"
  DiseasetermIC[DiseasetermIC[1]=="HP:0000005",]$ontology<-"MI"
  DiseasetermIC[DiseasetermIC[1]=="HP:0000118",]$ontology<-"PA"
  termsnumber<-length(termlist)
  for(i in 2:termsnumber)
  {
    message(paste(i,"of",termsnumber))
    xx<-as.character(DiseasetermIC[i,]$termlist)
    if(xx %in% PATerms){
      DiseasetermIC[DiseasetermIC[1]==xx,]$ontology<-"PA"
      DiseasetermIC[DiseasetermIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/PAdiseasenumber)
    }
    if(xx %in% OCTerms){
      DiseasetermIC[DiseasetermIC[1]==xx,]$ontology<-"OC"
      DiseasetermIC[DiseasetermIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/OCdiseasenumber)
    }
    if(xx %in% MITerms){
      DiseasetermIC[DiseasetermIC[1]==xx,]$ontology<-"MI"
      DiseasetermIC[DiseasetermIC[1]==xx,][3]<--log(termfrequency[termfrequency[1]==xx,][2]/MIdiseasenumber)
    }
  }
  save(DiseasetermIC,file="DiseaseTermIC.rda",compress="xz")
}

RemoveTermsWithoutIC <- function (terms,ontology,IC){
  info<-IC[IC[,1] %in% terms,] #所有条目对应的IC信息
  info<-info[info[,2]==ontology,] #筛选本体
  info<-info[info[,3]!=Inf,] #筛选IC的值
  return (as.character(info[,1]))
}
