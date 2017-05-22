##------------------------------------------------------------##
## PSI-MI 2.5 XML entry parsers
## Private low-level parsers
##------------------------------------------------------------##

##----------------------------------------##
## Experiment Parser
##----------------------------------------##
parseXmlExperimentNode <- function(root, namespaces, sourceDb) {
  subDoc <- xmlDoc(root)
  interactionDetectionMethod <- nonNullXMLvalueByPath(doc = subDoc,
                                           path = "/ns:experimentDescription/ns:interactionDetectionMethod/ns:names/ns:shortLabel",
                                           namespaces = namespaces)
  
  if(isLengthOneAndNA(interactionDetectionMethod)) {
    interactionDetectionMethod <- nonNullXMLvalueByPath(doc = subDoc,
                                             path = "/ns:experimentDescription/ns:interactionDetectionMethod/ns:names/ns:fullName",
                                             namespaces = namespaces)
  }
  
  ## it seems that xpathApply treats XPath in a case-insensitive way, to be confirmed
  expPubMed <- nonNullXMLattributeValueByPath(doc = subDoc,
                                              path = "/ns:experimentDescription/ns:bibref/ns:xref/ns:primaryRef[@db='pubmed']", 
                                              name = "id",
                                              namespaces = namespaces)
  
  ## experiment source Id
  sourceId <- nonNullXMLattributeValueByPath(doc = subDoc,
                                             path = paste("/ns:experimentDescription/ns:xref/ns:primaryRef[@db='",sourceDb,"']",sep=""),
                                             name="id", namespaces=namespaces)
  
  ## if sourceId not found, try alternatives
  if(isLengthOneAndNA(sourceId)) {
    sourceId <- nonNullXMLattributeValueByPath(doc = subDoc,
                                               path = paste("ns:experimentList/ns:experimentDescription/ns:xref/ns:secondaryRef[@db='",
                                                 sourceDb,
                                                 "']|/ns:experimentDescription",sep=""), 
                                               name = "id", namespaces = namespaces)
  }
  
  free(subDoc)
  experiment <- new("psimi25Experiment",
                    sourceDb = sourceDb,
                    interactionDetectionMethod = interactionDetectionMethod,
                    sourceId = sourceId,
                    expPubMed = expPubMed)
  return(experiment)
}

getExperimentNodeSetPath <- function(basePath) {
  path <- paste(basePath, "/ns:experimentList/ns:experimentDescription", sep = "", collapse = "")
  return(path)
}

getXmlExperimentNodeSet <- function(doc, basePath, namespaces) {
  experimentPath <- getExperimentNodeSetPath(basePath)
  experimentNodes <- getNodeSet(doc, experimentPath, namespaces)
}
parseXmlExperimentNodeSet <- function(nodes, psimi25source, namespaces, verbose) {
  if (verbose) {statusDisplay("  Parsing experiments:\n")}
  experimentEnv <- new.env(parent=emptyenv(), hash=TRUE)
  for (exp_ix in seq_along(nodes)) {
    if (verbose) {statusIndicator(exp_ix, length(nodes))}
    exper <- parseXmlExperimentNode(nodes[[exp_ix]], namespaces=namespaces, sourceDb=sourceDb(psimi25source))
    assign(xmlGetAttr(nodes[[exp_ix]], name="id"), exper, envir = experimentEnv)
  }
  if (verbose) statusDisplay("\n")
  return(experimentEnv)
}

##----------------------------------------##
## Interactor Parser
##----------------------------------------##
parseXmlInteractorNode <- function(root, namespaces, sourceDb, uniprotsymbol) {
  subDoc <- xmlDoc(root)
  sourceIds <- XMLattributeValueByPath(doc = subDoc, path = "/ns:interactor", 
                                              name = "id", namespaces = namespaces)
  shortLabels <- nonNullXMLvalueByPath(doc = subDoc,
                                       path="/ns:interactor/ns:names/ns:shortLabel",
                                       namespaces = namespaces)
  if(isLengthOneAndNA(shortLabels)) {
    shortLabels <- nonNullXMLvalueByPath(doc = subDoc,
                                          path="/ns:interactor/ns:names/ns:fullName",
                                          namespaces = namespaces)
  }
  
  organismNames <- nonNullXMLvalueByPath(doc = subDoc,
                                         path = "/ns:interactor/ns:organism/ns:names/ns:fullName",
                                         namespaces = namespaces)
  taxIds <- nonNullXMLattributeValueByPath(doc = subDoc,
                                           path = "/ns:interactor/ns:organism",
                                           name = "ncbiTaxId",
                                           namespaces = namespaces)
  xrefNodeSet <- getNodeSet(doc = subDoc,
                            path = "/ns:interactor/ns:xref/ns:primaryRef|/ns:interactor/ns:xref/ns:secondaryRef",
                            namespaces = namespaces)
  xrefs_df <- data.frame(
    kind = sapply(xrefNodeSet, xmlName),
    db = sapply(sapply(xrefNodeSet, xmlGetAttr, "db"), null2na),
    refType = sapply(sapply(xrefNodeSet, xmlGetAttr, "refType"), null2na),
    id = sapply(sapply(xrefNodeSet, xmlGetAttr, "id"), null2na),
    stringsAsFactors = FALSE
  )

  sourceRef <- subset(xrefs_df, db == sourceDb & refType == "identity")
  if (nrow(sourceRef) > 1L && "primaryRef" %in% sourceRef$kind) {
    sourceRef <- subset(sourceRef, kind=="primaryRef")
  }
  if (nrow(sourceRef) > 1L) { sourceRef <- sourceRef[1:1, ] }

  uniprotRef <- subset(xrefs_df, db == uniprotsymbol & refType == "identity")
  if (nrow(uniprotRef) > 1L && "primaryRef" %in% uniprotRef$kind) {
    uniprotRef <- subset(uniprotRef, kind=="primaryRef")
  }
  if (nrow(uniprotRef) > 1L) { uniprotRef <- uniprotRef[1:1, ] }

  free(subDoc)
  interactor <- new("psimi25Interactor",
                    sourceDb = sourceDb,
                    sourceId = sourceIds,
                    shortLabel = shortLabels,
                    uniprotId = uniprotRef$id,
                    organismName = organismNames,
                    taxId = taxIds,
                    xrefs = xrefs_df
                    )
  return(interactor)
}

##----------------------------------------##
## Interaction Parser
##----------------------------------------##
getInteractionNodeSetPath <- function(basePath) {
  path <- paste(basePath, "/ns:interactionList/ns:interaction", 
                sep = "", collapse = "")
  return(path)
}

getInteractionNodeSet <- function(doc, basePath, namespaces) {
  interactionPath <- getInteractionNodeSetPath(basePath)

  interactionNodes <- getNodeSet(doc=doc,
                                 path=interactionPath, 
                                 namespaces = namespaces)
  return(interactionNodes)
}

getInteractionPubmedPath <- function(sourceDb) {
  path <- paste("/ns:interaction/ns:xref/ns:primaryRef[@db='",sourceDb,"']",sep="")
  return(path)
}

isEmptyNodeSet <- function(nodeset) {
  return(length(nodeset) == 1 && is.na(nodeset))
}

parseXmlInteractionNode <- function(node,
                                    psimi25source,
                                    expEnv,
                                    interactorInfo,
                                    namespaces,
                                    verbose) {
  subDoc <- xmlDoc(node)
  sourceDb <- sourceDb(psimi25source)
  psimi25Id <- nonNullXMLattributeValueByPath(doc = subDoc,
                                              path = getInteractionPubmedPath(sourceDb),
                                              name = "id", namespaces = namespaces)[[1]]
  expRef <- XMLvalueByPath(doc = subDoc,
                           path = "/ns:interaction/ns:experimentList/ns:experimentRef",
                           namespaces = namespaces)
  
  if ((!is.null(expRef)) && exists(expRef, envir = expEnv)) {
    expData <- get(expRef, envir = expEnv)
    interactionDetectionMethod <- expData@interactionDetectionMethod
    expPsimi25 <- expData@sourceId
    expPubMed <- expData@expPubMed
  }
  else {
    interactionDetectionMethod <- nonNullXMLvalueByPath(doc = subDoc,
                                             path = "/ns:interaction/ns:experimentList/ns:experimentDescription/ns:interactionDetectionMethod/ns:names/ns:shortLabel",
                                             namespaces = namespaces)[[1]]
    if(isLengthOneAndNA(interactionDetectionMethod)) {
      interactionDetectionMethod <- nonNullXMLvalueByPath(doc = subDoc,
                                               path = "/ns:interaction/ns:experimentList/ns:experimentDescription/ns:interactionDetectionMethod/ns:names/ns:fullName",
                                               namespaces = namespaces)[[1]]
    }
    expPsimi25 <- nonNullXMLattributeValueByPath(doc = subDoc,
                                                 path = sprintf("/ns:interaction/ns:experimentList/ns:experimentDescription/ns:xref/ns:primaryRef[@db='%s']",sourceDb),
                                                 name = "id",
                                                 namespaces = namespaces)[[1]]
    expPubMed <- nonNullXMLattributeValueByPath(doc = subDoc,
                                                path = "/ns:interaction/ns:experimentList/ns:experimentDescription/ns:bibref/ns:xref/ns:primaryRef[@db='pubmed']",
                                                name = "id",
                                                namespaces = namespaces)[[1]]
  }
  interactionType <- nonNullXMLvalueByPath(doc = subDoc,
                                           path = "/ns:interaction/ns:interactionType/ns:names/ns:shortLabel",
                                           namespaces = namespaces)[[1]]
  ## misc attributes
  ## confidence value
  confidenceValue <- nonNullXMLvalueByPath(doc=subDoc, namespaces=namespaces,
                                           path = "/ns:interaction/ns:confidenceList/ns:confidence/ns:value")

  isNegative <- XMLvalueByPath(doc=subDoc, namespaces=namespaces,
                               path = "/ns:interaction/ns:negative")
  isNegative <- !is.null(isNegative) && length(isNegative) == 1 && switch(isNegative, true = TRUE, false = FALSE, stop("Unknown <negative> value: ", isNegative))
  isModeled <- XMLvalueByPath(doc=subDoc, namespaces=namespaces,
                               path = "/ns:interaction/ns:modeled")
  isModeled <- !is.null(isModeled) && length(isModeled) == 1 && switch(isModeled, true = TRUE, false = FALSE, stop("Unknown <modeled> value: ", isModeled))
  isIntraMolecular <- XMLvalueByPath(doc=subDoc, namespaces=namespaces,
                               path = "/ns:interaction/ns:intraMolecular")
  isIntraMolecular <- !is.null(isIntraMolecular) && length(isIntraMolecular) == 1 && switch(isIntraMolecular, true = TRUE, false = FALSE, stop("Unknown <intraMolecular> value: ", isIntraMolecular))

  ## participant
  rolePath <- "/ns:interaction/ns:participantList/ns:participant/ns:interactorRef"
  participantRefs <- nonNullXMLvalueByPath(doc=subDoc, path=rolePath, namespaces=namespaces)
  if(isEmptyNodeSet(participantRefs)) {
    rolePath <- "/ns:interaction/ns:participantList/ns:participant/ns:interactor"
    participantRefs <- nonNullXMLattributeValueByPath(doc=subDoc,
                                                      path=rolePath,
                                                      name="id", namespaces=namespaces)
  }
  
  ## bait
  rolePath <- "/ns:interaction/ns:participantList/ns:participant[ns:experimentalRoleList/ns:experimentalRole/ns:names/ns:fullName='bait']/ns:interactorRef"
  baitRefs <- nonNullXMLvalueByPath(doc=subDoc, path = rolePath, namespaces=namespaces)
  if(isEmptyNodeSet(baitRefs)) {
    rolePath <- "/ns:interaction/ns:participantList/ns:participant[ns:experimentalRoleList/ns:experimentalRole/ns:names/ns:fullName='bait']/ns:interactor"
    baitRefs <- nonNullXMLattributeValueByPath(doc=subDoc,
                                               path=rolePath,
                                               name="id", namespaces=namespaces)
  }
  
  ##prey
  rolePath <- "/ns:interaction/ns:participantList/ns:participant[ns:experimentalRoleList/ns:experimentalRole/ns:names/ns:fullName='prey']/ns:interactorRef"
  preyRefs <- nonNullXMLvalueByPath(doc=subDoc, path = rolePath, namespaces=namespaces)
  if(isEmptyNodeSet(preyRefs)) {
    rolePath <- "/ns:interaction/ns:participantList/ns:participant[ns:experimentalRoleList/ns:experimentalRole/ns:names/ns:fullName='prey']/ns:interactor"
    preyRefs <- nonNullXMLattributeValueByPath(doc=subDoc,
                                               path=rolePath,
                                               name="id", namespaces=namespaces)
  }
  
  ## inhibitor
  rolePath <- "/ns:interaction/ns:participantList/ns:participant[ns:experimentalRoleList/ns:experimentalRole/ns:names/ns:fullName='inhibitor']/ns:interactorRef"
  inhibitorRefs <- nonNullXMLvalueByPath(doc=subDoc, path=rolePath, namespaces=namespaces)
  
  ## neutral component
  rolePath <- "/ns:interaction/ns:participantList/ns:participant[ns:experimentalRoleList/ns:experimentalRole/ns:names/ns:fullName='neutral component']/ns:interactorRef"
  neutralComponentRefs <- nonNullXMLvalueByPath(doc=subDoc, path = rolePath, namespaces=namespaces)
  free(subDoc)
  
  srcLabel <- "sourceId"; uniLabel <- "uniprotId"
  participantUniprot <- getValueByMatchingMatrixColumn(participantRefs, interactorInfo, srcLabel, uniLabel)
  preyUniprot <- getValueByMatchingMatrixColumn(preyRefs, interactorInfo, srcLabel, uniLabel)
  baitUniprot <- getValueByMatchingMatrixColumn(baitRefs, interactorInfo, srcLabel, uniLabel)
  inhibitorUniprot <- getValueByMatchingMatrixColumn(inhibitorRefs,  interactorInfo, srcLabel, uniLabel)
  neutralComponentUniprot <- getValueByMatchingMatrixColumn(neutralComponentRefs,  interactorInfo, srcLabel, uniLabel)
  
  interaction <- new("psimi25Interaction",
                     sourceDb = sourceDb,
                     sourceId = as.character(psimi25Id), ## FIXME: can we do it nullable?
                     interactionType = interactionType,
                     interactionDetectionMethod = interactionDetectionMethod,
                     expPubMed = expPubMed,
                     ##expSourceId = expPsimi25, 
                     confidenceValue = confidenceValue,
                     participant = participantUniprot,
                     bait = baitRefs,
                     baitUniProt = baitUniprot, 
                     prey = preyRefs,
                     preyUniProt = preyUniprot,
                     inhibitor = inhibitorUniprot, 
                     neutralComponent = neutralComponentUniprot,
                     isNegative = isNegative,
                     isModeled = isModeled,
                     isIntraMolecular = isIntraMolecular
                     )
  return(interaction)
  
}

parseXmlInteractionNodeSet <- function(nodes,
                                       psimi25source,
                                       expEnv,
                                       interactorInfo,
                                       namespaces,
                                       verbose) {

  if (verbose) {statusDisplay("  Parsing interactions:\n")}

  interactions <- lapply(seq_along(nodes), function(node_ix) {
    if (verbose) {statusIndicator(node_ix, length(nodes))}
    parseXmlInteractionNode(nodes[[node_ix]],
                         psimi25source=psimi25source,
                         expEnv=expEnv,
                         interactorInfo=interactorInfo,
                         namespaces=namespaces,
                         verbose=verbose)
  })
  if (verbose) {statusDisplay("\n")}

  return(interactions)
}

##----------------------------------------##
## complex parser
## TODO: needs to know why 'as.character' is needed
##----------------------------------------##

parseXmlComplexNode <- function(node,
                                namespaces,
                                psimi25source) {

  subDoc <- xmlDoc(node)
  sourcedb <- sourceDb(psimi25source)
  sourceId <- nonNullXMLattributeValueByPath(doc=subDoc,
                                             path=paste("/ns:interaction/ns:xref/ns:primaryRef[@db='",sourcedb,"']",sep=""),
                                             name="id",
                                             namespaces=namespaces)

  shortLabel <- nonNullXMLvalueByPath(doc=subDoc,
                                      path="/ns:interaction/ns:names/ns:shortLabel",
                                      namespaces=namespaces)

  fullName <- nonNullXMLvalueByPath(doc=subDoc,
                                    path="/ns:interaction/ns:names/ns:fullName",
                                    namespaces=namespaces)

  interactorRef <- nonNullXMLvalueByPath(doc=subDoc,
                                         path="/ns:interaction/ns:participantList/ns:participant/ns:interactorRef",
                                         namespaces=namespaces)
  if(isEmptyNodeSet(interactorRef)) {
    interactorRef <- nonNullXMLattributeValueByPath(doc=subDoc,
                                                      path="/ns:interaction/ns:participantList/ns:participant/ns:interactor",
                                                      name="id", namespaces=namespaces)
  }

  attributesList <- parseXmlAttributesListByPath(doc=subDoc,
                                                 path="/ns:interaction/ns:attributeList/ns:attribute[@name]",
                                                 namespaces=namespaces)

  free(subDoc)
  complex <- new("psimi25Complex",
                 sourceDb=sourcedb,
                 sourceId=sourceId,
                 shortLabel=shortLabel,
                 fullName=fullName,
                 interactorRef=interactorRef,
                 attributesList=attributesList)
  return(complex)

}

annotateComplexesWithInteractors <- function(complexes,
                                              interactors) {
  interactorIds <- sapply(interactors, sourceId)
  interactorInfo <- interactorInfo(interactors)

  for(i in seq(along=complexes)) {
    thisComplex <- complexes[[i]]
    thisNodeInteractorIds <- interactorRef(thisComplex)
    participantRef <- split(thisNodeInteractorIds, 
                            as.factor(thisNodeInteractorIds))
    multiplicity <- sapply(participantRef, length)
    if ( !all(names(participantRef) %in% interactorIds) ) {
      msg <- paste("parse complex ", shortLabel(thisComplex), 
                   "can't resolve all interactor references.")
      stop(msg)
    }	
    participantIndex <- interactorIds %in% names(participantRef)
    participants <- interactorInfo[participantIndex, "sourceId"];
    participantsUniProt <- interactorInfo[participantIndex, "uniprotId"];
    names(participants) <- as.character(participants)
    participants <- data.frame(sourceId=participants, uniprotId = participantsUniProt, multiplicity=multiplicity)
    thisOrganism <- unique(as.character(interactorInfo[participantIndex,"organismName"]))
    thisTaxid <- unique(as.character(interactorInfo[participantIndex, "taxId"]))
    
    taxId(thisComplex) <- thisTaxid
    organismName(thisComplex) <- thisOrganism
    members(thisComplex) <- participants
    complexes[[i]] <- thisComplex
  }

  return(complexes)
}
  
##----------------------------------------##
## Entry Parser
##----------------------------------------##
getEntryBasePath <- function(index) {
  basePath <- paste("/ns:entrySet/ns:entry[", index, "]", sep = "", collapse = "")
  return(basePath)
}

getReleaseDatePath <- function(basePath) {
  releaseDatePath <- paste(basePath, "/ns:source", sep = "", collapse = "")
  return(releaseDatePath)
}

parseReleaseDate <- function(doc, basePath, namespaces) {
  releaseDatePath <- getReleaseDatePath(basePath)
  releaseDate <- nonNullXMLattributeValueByPath(doc = doc,
                                                path = releaseDatePath, 
                                                name = "releaseDate",
                                                namespaces = namespaces)
  return(releaseDate)
}
getPrimaryInteractorPath <- function(basePath) {
  interactorPath <- paste(basePath, "/ns:interactorList/ns:interactor", 
                          sep = "", collapse = "")
  return(interactorPath)
}


getSecondaryInteractorPath <- function(basePath) {
  interactorPath <-  paste(basePath, "/ns:interactionList/ns:interaction/ns:participantList/ns:participant/ns:interactor",
                           sep="", collapse="")
  return(interactorPath)
}

getInteractorNodeSet <- function(doc, basePath, namespaces) {
  interactorPath <- getPrimaryInteractorPath(basePath)
  interactorNodes <- getNodeSet(doc, interactorPath, namespaces)
  if(length(interactorNodes) == 0) { ## in case no interactorList was provided
    interactorPath <- getSecondaryInteractorPath(basePath)
    interactorNodes <- getNodeSet(doc, interactorPath, 
                                  namespaces)
  }
  return(interactorNodes)
}

parseXmlInteractorNodeSet <- function(nodes, psimi25source,
                                      namespaces, verbose) {
  interactorIds <- sapply(nodes, xmlGetAttr, name = "id")
  interactorCount <- length(nodes)
  if (verbose) {statusDisplay(paste(interactorCount, "interactor(s) found\n"))}
  if (verbose) {statusDisplay("  Parsing interactors:\n")}

  interactors <- vector("list",length=interactorCount)
  if (interactorCount > 0) {
    for (p in seq(interactorCount)) {
      if(verbose)
        statusIndicator(p, interactorCount)
      theRes <- parseXmlInteractorNode(root=nodes[[p]],
                                       namespaces=namespaces,
                                       sourceDb=sourceDb(psimi25source),
                                       uniprotsymbol=uniqueIdentifierSymbol(psimi25source))
      interactors[[p]] <- theRes
    }
  }
  if(verbose)
    statusDisplay("\n")
  
  interactorInfMat <- interactorInfo(interactors)
  names(interactors) <- interactorInfMat[,"uniprotId"]
  return(interactors)
}

## shortcut of getting interactor node set and parse them
parseXmlEntryInteractors <- function(doc,
                                     basePath,
                                     psimi25source,
                                     namespaces,
                                     verbose=TRUE) {
  interactorNodes <- getInteractorNodeSet(doc=doc,
                                          basePath=basePath,
                                          namespaces=namespaces)
  interactors <- parseXmlInteractorNodeSet(nodes=interactorNodes,
                                           psimi25source=psimi25source,
                                           namespaces=namespaces,
                                           verbose=verbose)
  return(interactors)
}

parseXmlEntryNode <- function(doc, index, namespaces, psimi25source, verbose=TRUE) {
  if(verbose) {statusDisplay(paste0("Parsing entry #",index,"\n"))}

  basePath <- getEntryBasePath(index)
  thisEntry <- new("psimi25InteractionEntry")

  ## experiment
  experimentNodes <- getXmlExperimentNodeSet(doc=doc, basePath=basePath,
                                             namespaces=namespaces)
  experimentEnv <- parseXmlExperimentNodeSet(nodes=experimentNodes, psimi25source=psimi25source, 
                                             namespaces=namespaces, verbose=verbose)
  
  ## misc information
  releaseDate(thisEntry) <- parseReleaseDate(doc=doc,
                                             basePath=basePath,
                                             namespaces=namespaces)
  
  ## interactor
  interactors <- parseXmlEntryInteractors(doc=doc,
                                          basePath=basePath,
                                          psimi25source=psimi25source,
                                          namespaces=namespaces,
                                          verbose=verbose)
  interactorInfMat <- interactorInfo(interactors)

  
  organismName <- unique(unlist(interactorInfMat[, "organismName"]))
  organismName(thisEntry) <- organismName
  taxId <- unique(unlist(interactorInfMat[, "taxId"]))
  taxId(thisEntry) <- taxId

  ## interaction
  interactionNodes <- getInteractionNodeSet(doc=doc,
                                            basePath=basePath,
                                            namespaces=namespaces)
  if (verbose) {statusDisplay(paste0(length(interactionNodes), " interaction(s) found\n"))}
  sourcedb <- sourceDb(psimi25source)
  interactions <- parseXmlInteractionNodeSet(nodes=interactionNodes,
                                             psimi25source = psimi25source,
                                             expEnv = experimentEnv,
                                             interactorInfo = interactorInfMat,
                                             namespaces=namespaces,
                                             verbose=verbose)
  
  interactions(thisEntry) <- interactions
  rownames(interactorInfMat) <- interactorInfMat[, "uniprotId"]
  interactors(thisEntry) <- interactors

  return(thisEntry)
}


parseXmlEntryNodeSet <- function(psimi25file, psimi25source, verbose=TRUE) {

  psimi25Doc <- xmlTreeParse(psimi25file, useInternalNodes = TRUE)
  
  psimi25NS <- getDefaultNamespace(psimi25Doc, simplify=TRUE)
  namespaces <- c(ns = psimi25NS)
  entry <- getNodeSet(psimi25Doc, "/ns:entrySet/ns:entry", namespaces)

  if (verbose) {statusDisplay(paste0(length(entry)," entries found\n"))}

  entryList <- lapply(seq_along(entry), function(i) {
    parseXmlEntryNode(doc=psimi25Doc, index=i,
                      namespaces=namespaces,
                      psimi25source=psimi25source,
                      verbose=verbose)
  })

  free(psimi25Doc)
  if (length(entryList) > 1) {
    el <- new("psimi25InteractionEntry")
    organismName(el) <- unique(unlist(sapply(entryList, organismName, simplify=FALSE)))
    taxId(el) <-  unique(unlist(sapply(entryList, taxId, simplify=FALSE)))
    releaseDate(el) <- unique(unlist(sapply(entryList, releaseDate,simplify=FALSE)))
    interactors(el) <- unique(unlist(sapply(entryList, interactors)))
    interactions(el) <- unique(unlist(sapply(entryList, interactions)))
    return(el)
  } else {
    return(entryList[[1]])
  }
}

##------------------------------------------------------------##
## High-level public parsers
##------------------------------------------------------------##

## File parser: parsing file into interaction entries
parsePsimi25Interaction <- function (psimi25file, psimi25source, verbose=TRUE) {
  if (verbose) message("Parsing PSI-MI 2.5 interactions file \"", psimi25file, "\"")
  parsedEntry <- parseXmlEntryNodeSet(psimi25file, psimi25source, verbose=verbose)
  return(parsedEntry)
}


## File parser: parsing file into complex
parsePsimi25Complex <- function(psimi25file, psimi25source, verbose=FALSE) {
  if (verbose) message("Parsing PSI-MI 2.5 complexes file \"", psimi25file, "\"")
  psiDoc <- xmlTreeParse(psimi25file, useInternalNodes=TRUE)
  psiNS <- xmlNamespaceDefinitions(psiDoc)
  namespaces <- c(ns=psiNS[[1]]$uri)

  entry <- getNodeSet(psiDoc, "/ns:entrySet/ns:entry", namespaces)
  if (verbose) {statusDisplay(paste0(length(entry)," entries found\n"))}
  if (length(entry) != 1L)
    stop("Internal RpsiXML Error: parsePsimi25Complex() does not support ", length(entry), " entries")

  basePath <- getEntryBasePath(1)

  releaseDate <- parseReleaseDate(doc=psiDoc,
                                  basePath=basePath,
                                  namespaces=namespaces)
  ## interactor
  interactors <- parseXmlEntryInteractors(doc=psiDoc,
                                          basePath=basePath,
                                          psimi25source=psimi25source,
                                          namespaces=namespaces,
                                          verbose=verbose)
  interactors <- interactors[ !duplicated(sapply(interactors, sourceId)) ] # remove the duplicates


  ## complex
  complexNodes <- getNodeSet(psiDoc, "//ns:interactionList/ns:interaction", namespaces)
  complexCount <- length(complexNodes)
  if (verbose) {statusDisplay(paste(complexCount, "complex(es) found\n"))}

  if (verbose) {statusDisplay("  Parsing complexes:\n")}
  complexList <- lapply(seq_along(complexNodes), function(cplxIx) {
    if (verbose) {statusIndicator(cplxIx, complexCount)}
    parseXmlComplexNode(complexNodes[[cplxIx]],
                        namespaces=namespaces,
                        psimi25source=psimi25source)
  })
  if (verbose) {statusDisplay("\n")}
  annotatedComplexList <- annotateComplexesWithInteractors(complexes=complexList,
                                                           interactors=interactors)
  
  if(verbose)
    statusDisplay("\n")
  free(psiDoc)
  
  new("psimi25ComplexEntry",
      interactors=interactors,
      complexes=annotatedComplexList,
      releaseDate=releaseDate)
}

## File parser: parsing file into graph
psimi25XML2Graph <- function(psimi25files,psimi25source,
                             type="interaction",
                             directed=TRUE,...) {
  stopifnot(type %in% c("interaction","complex"))

  if(type == "interaction"){
    result <- lapply(psimi25files, parsePsimi25Interaction, psimi25source,...)
    bpGraph <- interactionEntry2graph(result, directed=directed)
  }
  
  if (type == "complex"){
    result <- lapply(psimi25files, parsePsimi25Complex, psimi25source,...)
    bpGraph <- complexEntry2graph(result)
  }
  
  return(bpGraph)
}

interactionEntry2graph <- function(interactionEntry, directed=TRUE) {
  if(is(interactionEntry, "psimi25InteractionEntry")) {
    interactionEntry <- list(interactionEntry)
  }
  
  baitList <- lapply(interactionEntry, function(x){
    baits <- sapply(interactions(x), bait)
  })
    
  preyList <- lapply(interactionEntry, function(x){
    prey <- sapply(interactions(x), prey)
  })

  index <- sapply(baitList, class) == sapply(preyList, class)
  
  for(i in 1:length(index)){
    if(!index[i]){
      newBait <- vector(length=length(unlist(preyList[[i]])))
      k <- 1
      for(j in 1:length(preyList[[i]])){
        newBait[k:(k+length(preyList[[i]][[j]])-1)] <- rep(baitList[[i]][j], length(preyList[[i]][[j]]))
        k <- k + length(preyList[[i]][[j]])
      }
      baitList[[i]] <- newBait
    } 
  }
  
  
  b <- unlist(baitList); b[is.na(b)] <- "NA";
  p <- unlist(preyList); p[is.na(p)] <- "NA";
  
  bpList <- split(p,b)
  bpMat <- list2Matrix(bpList)
  
  bpG <- genBPGraph(bpMat, directed = directed)
  
  bpInteractors <- list()
  for(i in seq(interactionEntry)) {
    ints <- interactors(interactionEntry[[i]])
    newInts <- which(!names(ints) %in% bpInteractors)
    bpInteractors <- append(bpInteractors, ints[newInts])
  }
  bpGraph <- as(bpG, "psimi25Graph")
  bpGraph@interactors <- bpInteractors
  
  return(bpGraph)
}

complexEntry2graph <- function(complexEntry) {
  if(!is(complexEntry, "list")) {
    complexEntry <- list(complexEntry)
  }
  listOfListOfComps <- lapply(complexEntry, function(x){
    lapply(x@complexes, function(y){
      p <- as.vector(rep(y@members[,2], y@members[,3]))
      
      attr(p, "complexName") <- y@fullName 
      p
    })
  })
  he <- do.call(c, listOfListOfComps)
  nodes <- unique(unlist(he))
  hEdges <- lapply(he, function(x) {e <- Hyperedge(x, attr(x,"complexName"))})
  
  bpInteractors <- list()
  for(i in seq(complexEntry)) {
    ints <- interactors(complexEntry[[i]])
    newInts <- which(!names(ints) %in% bpInteractors)
    bpInteractors <- append(bpInteractors, ints[newInts])
  }
  
  bpGraph <- new("psimi25Hypergraph",
                 interactors=bpInteractors,
                 nodes=nodes,
                 hyperedges = hEdges)
  return(bpGraph)
}

buildPCHypergraph <- function(xmlFiles, psimi25source, split.by=c("none","organismName","taxId"), ...) {
  split.by <- match.arg(split.by)
  ie <- lapply(xmlFiles, parsePsimi25Complex, psimi25source, ...)
  hg <- complexEntry2graph(ie)

  if(split.by=="none")
    return(hg)
  
  hyperedges <- hyperedges(hg)
  interactors <- interactors(hg)
  hnodes <- nodes(hg)

  inOrgs <- sapply(interactors, organismName)
  inTax <- sapply(interactors, taxId)
  if(split.by == "organismName") {
    sf <- factor(inOrgs)
  } else if (split.by == "taxId") {
    sf <- factor(inTax)
  }

  hyperSf <- sapply(hyperedges, function(x) unique(sf[nodes(x)]))
  sfLevels <- levels(sf)

  hypers <- list()
  for(i in seq(along=sfLevels)) {
    le <- sfLevels[i]
    heOfLevel <- sapply(hyperSf, function(x) any(x %in% le))
    itOfLevel <- sf == le
    nodesOfLevel <- unique(unlist(sapply(hyperedges[heOfLevel],nodes)))
    
    hypers[[i]] <- new("psimi25Hypergraph",
                 interactors = interactors[itOfLevel],
                 nodes = nodesOfLevel,
                 hyperedges = hyperedges[heOfLevel])
                 
  }
  names(hypers) <- sfLevels

  return(hypers)
}

separateXMLDataByExpt <- function(xmlFiles, psimi25source, type = "direct", directed=TRUE, abstract=FALSE,...){
  

  if(!(type %in% c("direct","indirect", "eg"))){
    stop("The argument type needs to be either direct or indirect")
  }

  if(type == "direct"){
    interactionDetectionMethodWanted = c("two hybrid","two hybrid array", 
                               "2h fragment pooling", "2 hybrid",
                               "two hybrid pooling", "Two-hybrid")}
  if(type == "indirect"){
    interactionDetectionMethodWanted = c("coip","anti tag coip","pull down","tap",
        "anti bait coip","affinity chrom","chromatography","ion exchange chrom"
        ,"affinity techniques")}
  if(type == "eg"){
    interactionDetectionMethodWanted = c("spr")}
  
  #create a list of psimi25interaction objects corresponding to the list of xml files
  ie <- lapply(xmlFiles, parsePsimi25Interaction, psimi25source,...)

  
  #create an aggregate interactor list from all psimi25interactio objects
  interactorList <- lapply(ie, function(x) x@interactors)
  combinedList <- do.call(c, interactorList)
  protId <- unique(names(combinedList))
  uniqueCombList <- combinedList[protId]

  psiM <- lapply(ie, function(x){getDesired(x,interactionDetectionMethodWanted)})

  #now we have a kx3 matrix that has (bait,prey,pmid) info for all xml files
  psi <- do.call(rbind, psiM)

  if(is.null(psi)) {stop("There are no interactions of the type that you requested
                         within these XML files.")}

  #split the bait/prey info by pmid
  psiBPPairs <- split(data.frame(psi[,c("Bait","Prey")],
    stringsAsFactors=FALSE),psi[,"PMID"])

  #for each pmid, split the prey by the baits (we only want non-empty stuff)
  psiBPList <- lapply(psiBPPairs, function(x){split(x[,2], x[,1])})
  psiBPList <- lapply(psiBPList, function(x) {y = x[listLen(x)>0]; return(y)})

  #from each bpList for each pmid, we create adjMat and graphNEL
  psiBPMatrix <- lapply(psiBPList, list2Matrix)
  psiBPMatrix <- lapply(psiBPMatrix, t)
  psiBPGraphs <- lapply(psiBPMatrix, function(x) {genBPGraph(x, directed=directed)})

  
  
  #get abstract information for each dataset if desired
  if(abstract){
    pmids <- names(psiBPGraphs)
    abst <- getAbstractByPMID(pmids)

  #lastly, we create a psimi25Graph from each of the graphNELs
    psiGraphs <- mapply(function(x,w) {y <- as(x, "psimi25Graph");
                                       y@interactors <- uniqueCombList[nodes(x)];
                                       y@abstract = w; return(y)},
                        psiBPGraphs, abst)
  }

  else{

    psiGraphs <- lapply(psiBPGraphs, function(x,w) {
      y <- as(x, "psimi25Graph"); y@interactors <- uniqueCombList[nodes(x)];
      return(y)})

  }
  #return the graph objects
  return(psiGraphs)

  ###this works for intact, need to test on all other databases!!!
}

getDesired <- function(interactionEntry, intDetMethod, intType = NULL){

  x <- interactions(interactionEntry)
  
  dataL <- lapply(x, function(y){
    if((is.null(intDetMethod) || any(y@interactionDetectionMethod %in% intDetMethod)) &&
       (is.null(intType) || any(y@interactionType %in% intType))
    ) {
      z <- c(bait(y), prey(y), pubmedID(y), interactionType(y), interactionDetectionMethod(y));
      
      if(length(z)==5){
        if(!(is.na(z[1])) && !(is.na(z[2]))) {
          return(z)
        } else {
          return(NULL)
        }
      }
      
      if(length(z)>4){
        ## remove NA nodes
        bs <- bait(y); bs <- bs[!is.na(bs)]
        ps <- prey(y); ps <- ps[!is.na(ps)]
        ## if baits are all NAs, the item will be discarded, since later BP pairs will be indexed by bait
        if(all(is.na(bs))) {
          return(NULL)
        }
        
        nrow <- length(bs)*length(ps)
        bpm <- matrix(nrow=nrow, ncol=5)
        
        bpm[,1] <- rep(bs, each=length(ps));
        bpm[,2] <- rep(ps, length(bs))
        bpm[,3] <- rep(pubmedID(y), nrow);
        bpm[,4] <- rep(interactionType(y), nrow);
        bpm[,5] <- rep(interactionDetectionMethod(y), nrow);
        return(bpm)}
    }
  }
                  )
  
  dataL <- dataL[!(sapply(dataL, is.null))]
  if(length(dataL)>0){
    dataM <- do.call(rbind, dataL)
    colnames(dataM) <- c("Bait","Prey","PMID","Interaction Type","Interaction Detection Method")
    return(dataM)
  } else {
    return(NULL)
  }
  
}
  
