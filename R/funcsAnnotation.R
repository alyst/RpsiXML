getAbstractByPMID <- function(pmID){

pmXML <- pubmed(pmID)
pmRoot <- xmlRoot(pmXML)
numAbst <- length(pmRoot)
abstracts <- vector("list", length=numAbst)
for(i in 1:numAbst){
  abstracts[[i]] <- buildPubMedAbst(pmRoot[[i]])
}
names(abstracts) <- pmID
return(abstracts)

}
