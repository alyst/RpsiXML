# misc functions
genBPGraph <- function(bpMat, directed=TRUE, bp=TRUE){
  bpMat1 <- bpMat
  b <- rownames(bpMat)
  p <- colnames(bpMat)
  
  if(!bp){
    if(sum(b != p) != 0){
      stop("The rownames and the colnames must be identical.")
    }
  } else{
    baits <- union(rownames(bpMat), colnames(bpMat))
    preys <- baits
    
    bpMat1 <- matrix(0, length(baits), length(preys))
    dimnames(bpMat1) <- list(baits, preys)
    
    bpMat1[b,p] <- bpMat
    if(!directed) {
      bpMat1 <- bpMat1 + t(bpMat1)
      mode(bpMat1) <- "logical"
      mode(bpMat1) <- "numeric"
    }
  }
  
  bpGraph <- as(bpMat1, "graphNEL")
  if(!directed){
    bpGraph <- ugraph(removeSelfLoops(bpGraph))
  }
  bpGraph  
}

statusDisplay <- function(...) {
  cat(...)
}

statusIndicator <- function(x, length, N=40) {
  stages <- round(length/N + 0.5)
  if (x > length) {warning("Indicator received wrong message!\n"); x <- length}
  if (x %% stages == 0 | x == length) {
    per <- round(x/length,2)
    statusDisplay("\r  ",per*100, "% ", rep("=",round(N*per)), ">",sep="")
  }
}

isLengthOneAndNA <- function(x) {
  return(length(x) == 1 && is.na(x))
}

XMLvalueByPath <- function(doc, path, namespaces) {
  x <- unlist(xpathApply(doc=doc,path=path,fun=xmlValue, namespaces=namespaces))
  return(x)
}

nonNullXMLvalueByPath <- function(doc, path, namespaces) {
  x <- XMLvalueByPath(doc=doc,path=path, namespaces=namespaces)
  return(null2na(x))
}

XMLattributeValueByPath <- function(doc, path, name, namespaces) {
  x <- unlist(xpathApply(doc=doc,path=path,fun=xmlGetAttr, name=name, namespaces=namespaces))
  return(x)
}

nonNullXMLattributeValueByPath <- function(doc, path, name, namespaces) {
  x <- XMLattributeValueByPath(doc=doc, path=path, name=name, namespaces=namespaces)
  return(null2na(x))
}

getValueByMatchingMatrixColumn <- function(x, matrix, nameCol, valueCol) {
  names <- matrix[,nameCol]
  ind <- match(x, names)
  if (length(ind) == 1 && is.na(ind)) {
    value <- NA_character_
  } else {
    value <- unlist(matrix[ind, valueCol])
  }
  return(value)
}
