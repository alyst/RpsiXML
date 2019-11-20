##----------------------------------------------------------------------------##
##
## AllMethods: S4-class methods for RpsiXML package
##
## Author: Jitao David Zhang <jitao_david.zhang@roche.com>
##----------------------------------------------------------------------------##

##------------------------------------------------------------##
## auxilliary functions
##------------------------------------------------------------##
printManyItems <- function(title, x) {
  cat(title, ":", paste(x, collapse=", "), "\n")
}

quiteAsInteger <- function(x) {
  suppressWarnings(as.integer(x))
}

setMethod("name", "ANY", function(object) {
  object@name
})
setMethod("id", "ANY", function(object) {
  object@id
})
setReplaceMethod("id", c("ANY", "numeric"), function(object, value) {
  value <- as.integer(value)
  object@id <- value
  return(value)
})

setReplaceMethod("name", c("ANY","character"), function(object,value) {
  object@name <- value
  return(object)
})

## solve the generic conflicts between utils::head and hypergraph::head
setMethod("head", "ANY", function(.Object) {
  utils::head(.Object)
})
setMethod("tail", "ANY", function(.Object) {
  utils::tail(.Object)
})

##------------------------------------------------------------##
## sourceDbAndId methods
##------------------------------------------------------------##
setMethod("sourceDb", "sourceDbAndId", function(x) {
  x@sourceDb
})
setMethod("sourceId", "sourceDbAndId", function(x) {
  x@sourceId
})
setReplaceMethod("sourceDb", c("sourceDbAndId", "character"),
                 function(x, value) {
                   x@sourceDb <- value
                   return(x)
                 })
setReplaceMethod("sourceId", c("sourceDbAndId", "character"),
                 function(x, value) {
                   x@sourceId <- value
                   return(x)
                 })

##------------------------------------------------------------##
## organismIdAndName methods
##------------------------------------------------------------##
setMethod("taxId", "organismTaxIdAndName", function(x) {
  x@taxId
})
setMethod("organismName", "organismTaxIdAndName", function(x) {
  x@organismName
})
setReplaceMethod("taxId", c("organismTaxIdAndName", "character"), function(x, value) {
  x@taxId <- value
  return(x)
})
setReplaceMethod("organismName", c("organismTaxIdAndName", "character"), function(x, value) {
  x@organismName <- value
  return(x)
})

##------------------------------------------------------------##
## interactorListBase methods
##------------------------------------------------------------##
setMethod("interactors", "interactorListBase", function(x) {
  x@interactors
})
setMethod("numInteractors", "interactorListBase", function(x) {
  xit <- interactors(x)
  xitLength <- length(xit)
  return(xitLength)
})
setReplaceMethod("interactors", c("interactorListBase", "list"), function(x,value) {
  x@interactors <- value
  return(x)
})
setReplaceMethod("interactors", c("interactorListBase", "NULL"), function(x,value) {
  x@interactors <- list()
  return(x)
})

##------------------------------------------------------------##
## attributesList  methods
##------------------------------------------------------------##
setMethod("attributesList", "attributesListBase", function(x) {
  x@attributesList
})
setReplaceMethod("attributesList", c("attributesListBase", "list"), function(x, value) {
  x@attributesList <- value
  return(x)
})

##------------------------------------------------------------##
## psimi25Entry methods
##------------------------------------------------------------##
setMethod("releaseDate", "psimi25Entry", function(x) {
  x@releaseDate
})
setReplaceMethod("releaseDate", c("psimi25Entry", "character"), function(x,value) {
  x@releaseDate <- value
  return(x)
})

##------------------------------------------------------------##
## psimi25GraphBase methods
##------------------------------------------------------------##
setMethod("abstract", "psimi25GraphBase", function(object) {
  object@abstract
})

##------------------------------------------------------------##
## show methods
##------------------------------------------------------------##

setMethod("show", "psimi25Attribute", function(object) {
  cat(sprintf("%s: %s\n", object@name, object@attribute))
})

getHorizontalSepLine <- function(pattern="-", factor=0.8) {
  sepWidth <- round(options("width")[[1]] * factor)
  sep <- sprintf("%s\n",paste(rep(pattern, sepWidth), collapse=""))
  return(sep)
}
showHorizontalSepLine <- function(pattern="-", factor=0.8) {
  sepLine <- getHorizontalSepLine(pattern=pattern, factor=factor)
  cat(sepLine)
}

showListWithoutNames <- function(list, title="List", indent=2) {
  cat(title, "\n")
  for(i in seq(along=list)) {
    cat(rep(" ", indent))
    show(list[[i]])
  }
}

setMethod("show", "psimi25Interactor", function(object) {
  cat("interactor (", object@uniprotId, "):\n",
      getHorizontalSepLine(),
      "[ source database ]: ", object@sourceDb, "\n",
      "[ source ID ]: ", object@sourceId, "\n",
      "[ short label ]: ", object@shortLabel, "\n",
      "[ UniProt ID ]: ", object@uniprotId, "\n",
      "[ organism ]: ", paste(object@organismName, collapse=", "), "\n",
      "[ NCBI Taxonomy ID ]: ", object@taxId, "\n",
      "[ xrefs ]: ",nrow(xrefs(object)),"xrefs found, Use 'xrefs(x)' to see more annotation\n")
})

setMethod("show", "psimi25Interaction", function(object) {
  cat("interaction (", object@sourceId, "):\n",
      getHorizontalSepLine(),
      "[ source database ]:", object@sourceDb, "\n",
      "[ source experiment ID ]:", object@sourceId, "\n",
      "[ interaction type ]:", object@interactionType, "\n",
      "[ interaction detection method ]:", object@interactionDetectionMethod, "\n",
      "[ experiment ]: pubmed", object@expPubMed, "\n",
      "[ participants ]: there are", nrow(object@participants),
      "participants in total, here are the first few ones:\n")
  print(utils::head(participants(object)))
  cat("[ confidence value ]:", confidenceValue(object), "\n",
      "[ negative ]:", isNegative(object), "\n",
      "[ modeled ]:", isModeled(object), "\n",
      "[ intraMolecular ]:", isIntraMolecular(object), "\n"
      )
}
          )

setMethod("show", "psimi25InteractionEntry", function(object) {
  cat("==================================\n",
      "interaction entry (", object@releaseDate, "):\n",
      "==================================\n",
      "[ organism ]: ", paste(object@organismName, collapse=", "), "\n",
      "[ taxonomy ID ]: ", object@taxId, "\n",
      "[ interactors ]: there are ", numInteractors(object), 
      " interactors in total, here are the first few ones:\n")
  print(utils::head(interactorInfo(object)))
  cat("...\n", "[ interactions ]: there are ", 
      length(object@interactions), " interactions in total, here are the first few ones:\n")
  print(utils::head(object@interactions))
  cat("...\n")
  cat("Please use 'interactors' to view all interactors, and 'interactions' to view all interactions.\n")
})

setMethod("show","psimi25Complex", function(object) {
  cat(paste("complex (", object@sourceId, ")\n",
            getHorizontalSepLine(),
            "[ local ID ]: ", object@localId, "\n",
            "[ source database ]: ", object@sourceDb, "\n",
            "[ source ID ]: ", object@sourceId, "\n",
            "[ full name ]: ", object@fullName, "\n",
            "[ organism ]: ", paste(object@organismName, collapse=", "), "\n",
            "[ taxonomy ID ]: ", object@taxId, "\n",sep=""))
  showListWithoutNames(attributesList(object), title="[ attributesList ]:")
  cat("[ participants ]: \n")
  print(object@participants)
})

setMethod("show", "psimi25ComplexEntry", function(object) {
  cat("==================================\n",
      "complex entries (", object@releaseDate, "):\n",
      "==================================\n",
      "[ interactors ]: there are", length(object@interactors), 
      "interactors in total, here are the first few ones:\n")
  print(utils::head(interactorInfo(object)))
  cat("...\n", "[ complexes ]: there are", length(object@complexes), 
      "complexes in total, here are the first few ones:\n")
  print(utils::head(object@complexes))
  cat("...\n")
  cat("Please use 'interactors' to view all interactors, and use 'complexes' to view all complexes.\n")
})

          
setMethod("show", "psimi25Graph", function(object){
  print(as(object, "graph"))
  print(class(object))
})


setMethod("sourceDb", signature(x="psimi25Source"),
          function(x) x@sourceDb)

setMethod("uniqueIdentifierSymbol", signature("psimi25Source"),
          function(x) x@uniprotSymbol)

setMethod("uniprot", signature("psimi25Source"),
          function(x) x@uniprotSymbol)
setMethod("uniprot", signature("psimi25Interactor"),
          function(x) x@uniprotId)

setMethod("translateID", signature(r="psimi25Graph"),
          function(r, to){
            its <- interactors(r)
            
            trIds <- translateID(its, to=to)
            oldIds <- nodes(r)
            ##tc asks - is it safer to put nodes(r) = trIds[nodes(r)]
            ##since the names of trIds is the uniprot Ids?

            ## TODO: fix the "[[" 1: it should add a duplicate node automatically
            trNodes <- sapply(trIds[match(oldIds, names(trIds))], "[[", 1)

            ## for NAs: use the UniProt instead
            isFailedTranslation <- is.na(trNodes)
            trNodes[isFailedTranslation] <- oldIds[isFailedTranslation]
            
            newNodeNames <- trNodes
            ## handling duplicate nodes
            if(anyDuplicated(newNodeNames)) {
              manyToOne <- names(which(table(newNodeNames) > 1))
              toCollapse <- lapply(manyToOne, function(x) names(which(newNodeNames==x)))
              for(i in 1:length(toCollapse)) {
                r <- combineNodes(toCollapse[[i]], r, toCollapse[[i]][1])
              }
              newNodeNames <- newNodeNames[!duplicated(newNodeNames)]
            } 

            nodes(r) <- unname(newNodeNames)


            return(r)
          })



setMethod("translateID", signature(r="list"),
          function(r, to) {
            if (!all(sapply(r, is, "psimi25Interactor"))) {
              stop("'x' must be a list of psimi25Interactor objects!")
            }
            ## all available dbs
            ALLDB <- toupper(unique(unlist(sapply(r, xrefDbs))))
            TO <- toupper(to)
            if (!any(ALLDB %in% TO)) {
              stop(sprintf("%s not available! available xrefs are: %s\n",
                           TO,
                           paste(ALLDB, collapse=",")
                           ))
            }
            return(sapply(r, translateID, to=to))
          })

setMethod("translateID",
          signature(r="psimi25Interactor"),
          function(r, to) {
            if(missing(to)) {
              matchingRefs <- subset(xrefs(r), kind=="primaryRef")
            } else {
              matchingRefs <- subset(xrefs(r), grep(toupper(to), toupper(db)))
            }
            return(null2na(matchingRefs$id))
          })



setMethod("confidenceValue", "psimi25Interaction", function(x) {
   return(x@confidenceValue)
})

setMethod("isNegative", "psimi25Interaction", function(x) {
  return(x@isNegative)
})

setMethod("isModeled", "psimi25Interaction", function(x) {
  return(x@isModeled)
})

setMethod("isIntraMolecular", "psimi25Interaction", function(x) {
  return(x@isIntraMolecular)
})

### interactors
setMethod("xrefs", "psimi25Interactor", function(object) {
  return(get("xrefs", object@xrefs))
})
setMethod("xrefDbs", signature(x="psimi25Interactor"),
          function(x) {
            unique(xrefs(x)$db)
          })
setMethod("xrefDbs", signature(x="list"),
          function(x, intersect=FALSE) {
            isClass <- all(sapply(x, inherits, "psimi25Interactor"))
            if(!isClass) {
              stop("'x' must be a list of 'psimi25Interactor' objects!")
            }
            xrefs <- lapply(x, xrefDbs)
            nxrefs <- table(xrefs)
            if(!intersect) {
              return(names(nxrefs))
            } else {
              return(names(nxrefs)[nxrefs==length(x)])
            }
          })
setMethod("xrefDbs",
          signature(x="psimi25InteractionEntry"),
          function(x,...) {
            ints <- interactors(x)
            return(xrefDbs(ints, ...))
          })


setMethod("interactorInfo", signature(x="psimi25InteractionEntry"),
          function(x) {
            interactors <- interactors(x)
            interactorInfo(interactors)
            })
setMethod("interactorInfo", signature(x="psimi25Graph"),
          function(x) {
            interactors <- interactors(x)
            interactorInfo(interactors)
            })


setMethod("interactorInfo", signature(x="psimi25ComplexEntry"),
          function(x) {
            interactors <- interactors(x)
            interactorInfo(interactors)
            })

setMethod("interactorInfo", signature(x="list"),
          function(x) {
            stopifnot(all(sapply(x, is, "psimi25Interactor")))
            features <- c("localId", "sourceDb", "sourceId", "shortLabel",
                          "uniprotId",  "organismName", "taxId")
            res_list <- lapply(features, function(feat){
              sapply(x, function(xx) null2na(slot(xx, feat)))
            })
            names(res_list) <- features
            as.data.frame(res_list, stringsAsFactors = FALSE,
                          row.names = names(x))
          })

setMethod("interactions", signature(x="psimi25InteractionEntry"),
          function(x) x@interactions)
setReplaceMethod("interactions", signature(x="psimi25InteractionEntry", value="list"),
          function(x, value) {
            vlInteraction <- sapply(value, inherits, "psimi25Interaction")
            if(!all(vlInteraction))
              stop("'value' must be a list of psimi25Interaction")
            x@interactions <- value
            return(x)
            })

setMethod("complexes", signature(x="psimi25ComplexEntry"),
          function(x) x@complexes)


setMethod("participants", "psimi25Complex", function(x) {
  return(x@participants)
})
setReplaceMethod("participants", c("psimi25Complex","data.frame"), function(x,value) {
  x@participants <- value
  return(x)
})

## TODO: refactor to psimi25Names
setMethod("localId", "psimi25Complex", function(x) {
  return(x@localId)
})
setMethod("shortLabel", "psimi25Complex", function(x) {
  return(x@shortLabel)
})

setMethod("complexName", "psimi25Complex", function(x) {
  return(x@fullName)
})


## interactor accessors
setMethod("accession", "psimi25Interactor", function(x) {
  sourcedb <- sourceDb(x)
  accession <- translateID(x, sourceDb(x))
  return(accession)
})

## interaction accessors
setMethod("bait", "psimi25Interaction", function(x) {
  return(subset(x@participants, role=="bait")$uniprotId)
})
setMethod("prey", "psimi25Interaction", function(x)  {
  return(subset(x@participants, role=="prey")$uniprotId)
})
setMethod("participants", "psimi25Interaction", function(x) {
  return(x@participants)
})
setMethod("pubmedID", "psimi25Interaction", function(x) {
  return(x@expPubMed)
})



setMethod("interactionType", signature(object="psimi25Interaction"),
          function(object) object@interactionType)

setMethod("interactionDetectionMethod", signature(object="psimi25Interaction"),
          function(object) object@interactionDetectionMethod)

setMethod("numInteractions", signature(x="psimi25InteractionEntry"), function(x) {
  xit <- interactions(x)
  return(length(xit))
})

setMethod("pubmedID", signature(x="psimi25InteractionEntry"), function(x) {
  xit <- interactions(x)  
  pubmedid <- sapply(xit, pubmedID)
  uniquePubmedID <- unique(pubmedid)
  return(uniquePubmedID)
})

##--------------------##
## hypergraph methods
##--------------------##
setMethod("initialize",
          signature=signature(
            .Object="psimi25Hypergraph"),
          function(.Object, interactors, ...) {
            .Object@interactors <- interactors
            callNextMethod(.Object, ...)
          })

setMethod("hyperedgeNodes", "Hypergraph", function(x) {
  hyperEdges <- hyperedges(x)
  res <- lapply(hyperEdges, nodes)
  names(res) <- hyperedgeLabels(x)
  return(res)
})

setMethod("edgeLabel", "psimi25Hypergraph",function(x) {
  hyperedgeLabels(x)
})

setMethod("complexes", "psimi25Hypergraph", function(x) {
  hyperedgeNodes(x)
})


setMethod("numEdges", "psimi25Hypergraph", function(object) {
  return(length(hyperedges(object)))
})

setMethod("show", "psimi25Hypergraph", function(object){
  cat(paste("========================\n",
      "psimi25Hypergraph\n",
      "========================\n",
      "[ complexes ]: there are ", numEdges(object), " complexes (hyperedges) with ", numNodes(object) ," proteins (nodes) in the hypergraph, here are the first new ones:\n",sep=""))
  print(utils::head(complexes(object),3))
  cat("...\n")
  cat("Please use 'edgeLabel' to list only complexes names, or use 'complexes' to view all the complexes in list form\n")
})

setMethod("interactorInfo", signature(x="psimi25Hypergraph"),
          function(x) {
            interactors <- interactors(x)
            interactorInfo(interactors)
          })
setMethod("translateID", signature(r="psimi25Hypergraph"),
          function(r, to, uniprotId){
            its <- interactors(r)
            if(missing(uniprotId)) {
              return(translateID(its, to))
            }
            if(length(uniprotId) == 1 && is.na(uniprotId)) {
              nits <- is.na(names(its))
              return(translateID(its[nits], to))
            } else if (!uniprotId %in% names(its)) {
              warning("'uniprotId':", uniprotId, " not found in the hypergraph!")
              return(NULL)
            } else {
              translateID(its[[uniprotId]], to)
            }
          })

setMethod("revInciMat", signature(x="matrix"), function(x) {
  hyperNodes <- rownames(x)
  hyperEdgeNames <- colnames(x)
  hEdges <- list()
  for(i in seq(along=hyperEdgeNames)) {
    hEdges[[i]] <- hyperNodes[x[,i] == 1]
  }
  names(hEdges) <- hyperEdgeNames
  hyperEdges <- lapply(hEdges, "Hyperedge")

  hg <- new("Hypergraph", nodes=hyperNodes, hyperedges=hyperEdges)
  return(hg)
})


##------------------------------------------------------------##
## Methods for new XML interfaces
##------------------------------------------------------------##

####------------------------------------------------------------##
#### psimi25NamesType (and psimi25NamesAlis)
####------------------------------------------------------------##
##setMethod("psimi25NamesAlias",
##          c("character", "ANY","ANY"),
##          function(iValue, type, typeAc) {
##            obj <- new("psimi25NamesAlias",
##                       .Data=iValue,
##                       type=type,
##                       typeAc=typeAc)
##            return(obj)
##          })
##
##setMethod("psimi25NamesType",
##          c("ANY","ANY","ANY"),
##          function(shortLabel, fullName, alias) {
##            if(missing(shortLabel))
##              shortLabel <- as.character(NA)
##            if(missing(fullName))
##              fullName <- as.character(NA)
##            if(missing(alias))
##              alias <- list()
##            
##            obj <- new("psimi25NamesType",
##                       shortLabel=shortLabel,
##                       fullName=fullName,
##                       alias=alias)
##            return(obj)
##          })
##
##setMethod("show", "psimi25NamesType", function(object) {
##  printManyItems("Short Label", object@shortLabel)
##  printManyItems("Full name", object@fullName)
##  printManyItems("Alias", object@alias)
##
##})
##
##
####------------------------------------------------------------##
#### psimi25Attribute
####------------------------------------------------------------##
##setMethod("psimi25Attribute",
##          c("character","ANY","ANY"),
##          function(iValue, name, nameAc) {
##            if(missing(name))
##              name <- as.character(NA)
##            if(missing(nameAc))
##              nameAc <- as.character(NA)
##            obj <- new("psimi25Attribute",
##                       .Data=iValue,
##                       name=name,
##                       nameAc=nameAc)
##            return(obj)
##          })
##setMethod("show", "psimi25Attribute", function(object) {
##  cat(object@name, "[", object@nameAc, "]", ":", object@.Data, "\n")
##})
##setMethod("name", "psimi25Attribute", function(object) {
##  return(object@name)
##})
##setReplaceMethod("name", c("psimi25Attribute", "character"), function(object,value) {
##  object@name <- value
##  return(object)
##})
##setMethod("nameAc", "psimi25Attribute", function(object) {
##  return(object@nameAc)
##})
##setReplaceMethod("nameAc", "psimi25Attribute", function(object,value) {
##  object@nameAc <- value
##  return(object)
##})
##setMethod("iValue", "psimi25Attribute", function(object) {
##  return(object@.Data)
##})
##setReplaceMethod("iValue", "psimi25Attribute", function(object,value) {
##  object@.Data <- value
##  return(object)
##})
##
####setMethod("psimi25AttributeListType",
####          "list",
####          function(list) {
####            obj <- new("psimi25AttributeListType",
####                       .Data=list)
####            return(obj)
####          })
##
####------------------------------------------------------------##
#### psimi25DbReferenceType
####------------------------------------------------------------##
##setMethod("psimi25DbReferenceType",
##          c("ANY", "ANY", "ANY","ANY","ANY","ANY","ANY","ANY"),
##          function(list, db, dbAc, id, secondary, version, refType, refTypeAc) {
##            if(missing(list)) list <- list()
##            if(missing(db)) db <- as.character(NA)
##            if(missing(dbAc)) dbAc <- as.character(NA)
##            if(missing(id)) id <- as.character(NA)
##            if(missing(secondary)) secondary <- as.character(NA)
##            if(missing(version)) version <- as.character(NA)
##            if(missing(refType)) refType <- as.character(NA)
##            if(missing(refTypeAc)) refTypeAc <- as.character(NA)
##
##            obj <- new("psimi25DbReferenceType",
##                       .Data=list,
##                       db=db, dbAc=dbAc, id=id, secondary=secondary,
##                       version=version, refType=refType, refTypeAc=refTypeAc)
##          })
##
##setMethod("attrInfo", "psimi25DbReferenceType", function(object) {
##   att <- c("db","dbAc","id","secondary",
##                                 "version","refType","refTypeAc")
##   info <- data.frame(Value=c(object@db, object@dbAc, object@id, object@secondary,
##                        object@version, object@refType, object@refTypeAc),
##                      row.names=att)
##   return(info)
##})
##
##setMethod("show", "psimi25DbReferenceType",
##          function(object) {
##            cat("Attributes:\n")
##            att <- attrInfo(object)
##            show(format(att, justify="left"))
##            show(object@.Data)
##          })
##
####------------------------------------------------------------##
#### psimi25XrefType
####------------------------------------------------------------##
##setMethod("psimi25XrefType",
##          c("psimi25DbReferenceType","list"),
##          function(primaryRef,secondaryRef) {
##            obj <- new("psimi25XrefType",
##                       primaryRef=primaryRef,
##                       secondaryRef=secondaryRef)
##            return(obj)
##          })
##
####------------------------------------------------------------##
#### psimi25BibrefType
####------------------------------------------------------------##
##setMethod("psimi25BibrefType",
##          c("psimi25XrefType","list"),
##          function(xref, attributeList) {
##            obj <- new("psimi25BibrefType",
##                       xref=xref,
##                       attributeList=attributeList)
##            return(obj)
##          })
##
####------------------------------------------------------------##
#### psimi25Availability
####------------------------------------------------------------##
##setMethod("psimi25AvailabilityType",
##          c("character","ANY"),
##          function(iValue, id) {
##            suppressWarnings(id <- as.integer(id)) ## reason: some repositories do not adhere to the specification of id
##            obj <- new("psimi25AvailabilityType",
##                       .Data=iValue, id=id)
##            return(obj)
##          })
##setMethod("show",
##          "psimi25AvailabilityType",
##          function(object) {
##            printManyItems(object@id, object@.Data)
##          })
##setMethod("iValue", "psimi25AvailabilityType",
##          function(object) {
##            return(object@.Data)
##          })
##setReplaceMethod("iValue", "psimi25AvailabilityType",
##                 function(object,value) {
##                   object@.Data <- value
##                   return(object)
##                 })
##
####------------------------------------------------------------##
#### psimi25CvType and psimi25OpenCvType
####------------------------------------------------------------##
##setMethod("name", "psimi25CommonNameRef", function(object) {
##  object@name
##})
##
##setReplaceMethod("name", c("psimi25CommonNameRef","psimi25NamesType"), function(object,value) {
##  object@name <- value
##  return(object)
##})
##
##setMethod("xref", "psimi25CommonNameRef", function(object) {
##  object@xref
##})
##
##setReplaceMethod("xref", c("psimi25CommonNameRef","psimi25XrefType"), function(object,value) {
##  object@xref <- value
##  return(object)
##})
##
##setMethod("attributeList", "psimi25CommonNameRefAttr", function(object) {
##  object@attributeList
##})
##
##
##setReplaceMethod("attributeList", c("psimi25CommonNameRefAttr","list"), function(object,value) {
##  object@attributeList <- value
##  return(object)
##})
##
##
##setMethod("psimi25CvType",
##          c("psimi25NamesType","psimi25XrefType"),
##          function(name, xref) {
##            obj <- new("psimi25CvType",
##                       name=name,
##                       xref=xref)
##            return(obj)
##          })
##setMethod("psimi25OpenCvType",
##          c("psimi25NamesType", "psimi25XrefType","list"),
##          function(name, xref, attributeList) {
##            obj <- new("psimi25OpenCvType",
##                       name=name,
##                       xref=xref,
##                       attributeList=attributeList)
##            return(obj)
##          })
##
####------------------------------------------------------------##
#### psimi25BiosourceType
####------------------------------------------------------------##
##setMethod("psimi25BioSourceType",
##          c("psimi25NamesType","psimi25OpenCvType","psimi25OpenCvType", "psimi25OpenCvType", "numeric"),
##          function(name, cellType, compartment, tissue, ncbiTaxId) {
##            ncbiTaxId <- as.integer(ncbiTaxId) ## in case of numeric
##            obj <- new("psimi25BioSourceType",
##                       name=name,
##                       cellType=cellType,
##                       compartment=compartment,
##                       tissue=tissue,
##                       ncbiTaxId=ncbiTaxId)
##            return(obj)
##          })
##setMethod("cellType", "psimi25BioSourceType", function(object) object@cellType)
##setMethod("compartment", "psimi25BioSourceType", function(object) object@compartment)
##setMethod("tissue", "psimi25BioSourceType", function(object) object@tissue)
##setMethod("ncbiTaxId", "psimi25BioSourceType", function(object) object@ncbiTaxId)
##
##setReplaceMethod("cellType", c("psimi25BioSourceType","psimi25OpenCvType"),
##                 function(object,value) { object@cellType <- value; return(object) })
##setReplaceMethod("compartment", c("psimi25BioSourceType","psimi25OpenCvType"),
##                 function(object,value) { object@compartment <- value; return(object) })
##setReplaceMethod("tissue", c("psimi25BioSourceType","psimi25OpenCvType"),
##                 function(object,value) { object@tissue <- value; return(object) })
##setReplaceMethod("ncbiTaxId", c("psimi25BioSourceType","numeric"),
##                 function(object,value) {
##                   value <- as.integer(value)
##                   object@ncbiTaxId <- value;
##                   return(object)
##                 })
##
####------------------------------------------------------------##
#### psimi25ConfidenceType
####------------------------------------------------------------##
##setMethod("psimi25ConfidenceType",
##          c("psimi25OpenCvType","character"),
##          function(unit, value) {
##            obj <- new("psimi25ConfidenceType",
##                       unit=unit,
##                       value=value)
##            return(obj)
##          })
##setMethod("unit", "psimi25ConfidenceType", function(object) object@unit)
##setMethod("value", "psimi25ConfidenceType", function(object) object@value)
##
##setReplaceMethod("unit", c("psimi25ConfidenceType","psimi25OpenCvType"),
##                 function(object,value) { object@unit <- value; return(object) })
##setReplaceMethod("value", c("psimi25ConfidenceType","character"),
##                 function(object,value) { object@value <- value; return(object) })
##
####------------------------------------------------------------##
#### psimi25InteractorElementType
####------------------------------------------------------------##
##setMethod("psimi25InteractorElementType",
##          c("psimi25NamesType","psimi25XrefType","list",
##            "psimi25CvType","psimi25BioSourceType","character", "numeric"),
##          function(name, xref, attributeList,
##                   interactorType, organism, sequence, id) {
##            obj <- new("psimi25InteractorElementType",
##                       name=name,
##                       xref=xref,
##                       attributeList=attributeList,
##                       interactorType=interactorType,
##                       organism=organism,
##                       sequence=sequence,
##                       id=id)
##            return(obj)
##          })
##
##setMethod("interactorType", "psimi25InteractorElementType", function(object) object@interactorType)
##setMethod("organism", "psimi25InteractorElementType", function(object) object@organism)
##setMethod("Sequence", "psimi25InteractorElementType", function(object) object@sequence)
##
##setReplaceMethod("interactorType", c("psimi25InteractorElementType","psimi25CvType"),
##                 function(object,value) { object@interactorType <- value; return(object) })
##setReplaceMethod("organism", c("psimi25InteractorElementType","psimi25BioSourceType"),
##                 function(object,value) { object@organism <- value; return(object) })
##setReplaceMethod("Sequence", c("psimi25InteractorElementType","character"),
##                 function(object,value) { object@sequence <- value; return(object) })
##
####------------------------------------------------------------##
#### psimi25ParticipantType
####------------------------------------------------------------##
##setMethod("psimi25HostOrganism",
##          c("psimi25BioSourceType","list"),
##          function(bioSourceType, experimentRefList) {
##            obj <- new("psimi25HostOrganism",
##                       bioSourceType,
##                       experimentRefList=experimentRefList)
##            return(obj)
##          })
##setMethod("interactorRef", "psimi25ParticipantType", function(object) object@interactorRef)
##setMethod("interactor", "psimi25ParticipantType", function(object) object@interactor)
##setMethod("interactionRef", "psimi25ParticipantType", function(object) object@interactionRef)
##setMethod("participantIdentificationMethodList", "psimi25ParticipantType", function(object) object@participantIdentificationMethodList)
##setMethod("biologicalRole", "psimi25ParticipantType", function(object) object@biologicalRole)
##setMethod("experimentalRoleList", "psimi25ParticipantType", function(object) object@experimentalRoleList)
##setMethod("experimentalPreparationList", "psimi25ParticipantType", function(object) object@experimentalPreparationList)
##setMethod("experimentalInteractorList", "psimi25ParticipantType", function(object) object@experimentalInteractorList)
##setMethod("featureList", "psimi25ParticipantType", function(object) object@featureList)
##setMethod("hostOrganismList", "psimi25ParticipantType", function(object) object@hostOrganismList)
##setMethod("confidenceList", "psimi25ParticipantType", function(object) object@confidenceList)
##setMethod("parameterList", "psimi25ParticipantType", function(object) object@parameterList)
##
##setReplaceMethod("interactorRef", c("psimi25ParticipantType","numeric"),
##                 function(object,value) {
##                   value <- as.integer(value)
##                   object@interactorRef <- value; return(object)
##                 })
##setReplaceMethod("interactor", c("psimi25ParticipantType","psimi25InteractorElementType"),
##                 function(object,value) {
##                   object@interactor <- value; return(object)
##                 })
##setReplaceMethod("interactionRef", c("psimi25ParticipantType","integer"),
##                 function(object,value) {
##                   value <- as.integer(value)
##                   object@interactionRef <- value; return(object)
##                 })
##setReplaceMethod("participantIdentificationMethodList", c("psimi25ParticipantType","list"),
##                 function(object,value) {
##                   object@participantIdentificationMethodList <- value; return(object)
##                 })
##setReplaceMethod("biologicalRole", c("psimi25ParticipantType","psimi25CvType"),
##                 function(object,value) {
##                   object@biologicalRole <- value; return(object)
##                 })
##setReplaceMethod("experimentalRoleList", c("psimi25ParticipantType","list"),
##                 function(object,value) {
##                   object@experimentalRoleList <- value; return(object)
##                 })
##setReplaceMethod("experimentalPreparationList", c("psimi25ParticipantType","list"),
##                 function(object,value) {
##                   object@experimentalPreparationList <- value; return(object)
##                 })
##setReplaceMethod("experimentalInteractorList", c("psimi25ParticipantType","psimi25ExperimentInteractor"),
##                 function(object,value) {
##                   object@experimentalInteractorList <- value; return(object)
##                 })
##setReplaceMethod("featureList", c("psimi25ParticipantType","list"),
##                 function(object,value) {
##                   object@featureList <- value; return(object)
##                 })
##setReplaceMethod("hostOrganismList", c("psimi25ParticipantType","list"),
##                 function(object,value) {
##                   object@hostOrganismList <- value; return(object)
##                 })
##setReplaceMethod("confidenceList", c("psimi25ParticipantType","list"),
##                 function(object,value) {
##                   object@confidenceList <- value; return(object)
##                 })
##setReplaceMethod("parameterList", c("psimi25ParticipantType","list"),
##                 function(object,value) { object@parameterList <- value; return(object) })
##
####------------------------------------------------------------##
#### psimi25InteractionElementType
####------------------------------------------------------------##
##setMethod("psimi25ExperimentRefListType", "numeric", function(object) {
##  object <- as.integer(object)
##  obj <- new("psimi25ExperimentRefListType", .Data=object)
##  return(obj)
##})
##setMethod("psimi25ExperimentInteractor",
##          c("numeric","psimi25InteractorElementType","list"),
##          function(interactorRef, interactor, experimentRefList) {
##            obj <- new("psim25ExperimentInteractor",
##                       interactorRef=interactorRef,
##                       interactor=interactor,
##                       experimentRefList=experimentRefList)
##            return(obj)
##          })
##
##setMethod("participantList", "psimi25InteractionElementType", function(object) object@participantList)
##setMethod("inferredInteractionList", "psimi25InteractionElementType", function(object) object@inferredInteractionList)
##setMethod("interactionType", "psimi25InteractionElementType", function(object) object@interactionType)
##setMethod("modelled", "psimi25InteractionElementType", function(object) object@modelled)
##setMethod("experimentList", "psimi25InteractionElementType", function(object) object@experimentList)
##setMethod("intraMolecular", "psimi25InteractionElementType", function(object) object@intraMolecular)
##setMethod("negative", "psimi25InteractionElementType", function(object) object@negative)
##setMethod("confidenceList", "psimi25InteractionElementType", function(object) object@confidenceList)
##setMethod("parameterList", "psimi25InteractionElementType", function(object) object@parameterList)
##setMethod("imexId", "psimi25InteractionElementType", function(object) object@imexId)
##
##setReplaceMethod("participantList", c("psimi25InteractionElementType","list"),
##                 function(object,value) { object@participantList <- value; return(object) })
##setReplaceMethod("inferredInteractionList", c("psimi25InteractionElementType","psimi25InferredInteraction"),
##                 function(object,value) { object@inferredInteractionList <- value; return(object) })
##setReplaceMethod("interactionType", c("psimi25InteractionElementType","psimi25CvType"),
##                 function(object,value) { object@interactionType <- value; return(object) })
##setReplaceMethod("modelled", c("psimi25InteractionElementType","logical"),
##                 function(object,value) { object@modelled <- value; return(object) })
##setReplaceMethod("intraMolecular", c("psimi25InteractionElementType","logical"),
##                 function(object,value) { object@intraMolecular <- value; return(object) })
##setReplaceMethod("negative", c("psimi25InteractionElementType","logical"),
##                 function(object,value) { object@negative <- value; return(object) })
##setReplaceMethod("confidenceList", c("psimi25InteractionElementType","list"),
##                 function(object,value) { object@confidenceList <- value; return(object) })
##setReplaceMethod("parameterList", c("psimi25InteractionElementType","list"),
##                 function(object,value) { object@parameterList <- value; return(object) })
##setReplaceMethod("imexId", c("psimi25InteractionElementType","character"),
##                 function(object,value) { object@imexId <- value; return(object) })
##
##setMethod("psimi25ExperimentList", c("numeric","psimi25ExperimentType"),
##          function(experimentRef, experimentDescription) {
##            experimentRef <- suppressWarnings(as.integer(experimentRef))
##            obj <- new("psimi25ExperimentList",
##                       experimentRef=experimentRef,
##                       experimentDescription=experimentDescription)
##            return(obj)
##          })
####------------------------------------------------------------##
#### psimi25ExperimentType
####------------------------------------------------------------##
##setMethod("psimi25CvExperimentRefs",
##          c("psimi25CvType","list"),
##          function(cv, experimentRefList) {
##            obj <- new("psimi25CvExperimentRefs",
##                       cv,
##                       experimentRefList=experimentRefList)
##            return(obj)
##          })
##setMethod("psimi25CvExperimentRefsList",
##          "list",
##          function(object) {
##            obj <-  new("psimi25CvExperimentRefsList", object)
##            return(obj)
##          })
##
##setMethod("psimi25ExperimentType",
##          signature=signature("ANY","ANY","ANY",
##            "ANY","ANY","ANY",
##            "ANY","ANY","ANY", "numeric"),
##          function(name, bibref, xref, hostOrganismList,
##                   interactionDetectionMethod, participantIdentificationMethod,
##                   featureDetectionMethod, confidenceList,attributeList,id) {
##            if(missing(id))
##              id <- as.integer(NA)
##            id <- as.integer(id)
##            if(missing(name))
##              name <- new("psimi25NamesType")
##            if(missing(bibref))
##              bibref <- new("psimi25BibrefType")
##            if(missing(xref))
##              xref <- new("psimi25XrefType")
##            if(missing(bibref))
##              bibref <- new("psimi25XrefType")
##            if(missing(hostOrganismList))
##              hostOrganismList <- new("list")
##            if(missing(interactionDetectionMethod))
##              interactionDetectionMethod <- new("psimi25CvType")
##            if(missing(participantIdentificationMethod))
##              participantIdentificationMethod <- new("psimi25CvType")
##            if(missing(featureDetectionMethod))
##              featureDetectionMethod <- new("psimi25CvType")
##            if(missing(confidenceList))
##              confidenceList <- new("list")
##            if(missing(attributeList))
##              attributeList <- new("list")
##            
##            obj <- new("psimi25ExperimentType",
##                       name=name,
##                       bibref=bibref,
##                       xref=xref,
##                       hostOrganismList=hostOrganismList,
##                       interactionDetectionMethod=interactionDetectionMethod,
##                       participantIdentificationMethod=participantIdentificationMethod,
##                       featureDetectionMethod=featureDetectionMethod,
##                       confidenceList=confidenceList,
##                       attributeList=attributeList,
##                       id=id)
##            return(obj)
##          })
##
####------------------------------------------------------------##
#### psimi25Source
####------------------------------------------------------------##
##setMethod("psimi25Source",
##          c("ANY","ANY","ANY","ANY","character","character"),
##          function(name, xref, bibref, attributeList, release, releaseDate ){
##            if(missing(name)) name <- new("psimi25NamesType")
##            if(missing(xref)) xref <- new("psimi25XrefType")
##            if(missing(bibref)) bibref <- new("psimi25BibrefType")
##            if(missing(attributeList)) attributeList <- new("list")
##            obj <- new("psimi25Source",
##                       name=name, xref=xref,
##                       bibref=bibref, attributeList=attributeList,
##                       release=release, releaseDate=releaseDate)
##            return(obj)
##                       
##          })
####------------------------------------------------------------##
#### entry
####------------------------------------------------------------##
##setMethod("psimi25Entry",
##          c("psimi25Source","list","list",
##            "list","list",
##            "list"),
##          function(source, availabilityList, experimentList,
##                   interactorList, interactionList, attributeList) {
##            obj <- new("psimi25Entry",
##                       source=source,
##                       availabilityList=availabilityList,
##                       experimentList=experimentList,
##                       interactorList=interactorList,
##                       interactionList=interactionList,
##                       attributeList= attributeList)
##            return(obj)
##          })
##
##setMethod("level", "psimi25Entry", function(object) object@level)
##setMethod("version", "psimi25Entry", function(object) object@version)
##setMethod("minorVersion", "psimi25Entry", function(object) object@minorVersion)
##
##setReplaceMethod("level", c("psimi25Entry", "integer"), function(object, value) {
##  object@level <- value
##  return(object)
##})
##setReplaceMethod("version", c("psimi25Entry", "integer"), function(object, value) {
##  object@version <- value
##  return(object)
##})
##setReplaceMethod("minorVersion", c("psimi25Entry", "integer"), function(object, value) {
##  object@minorVersion <- value
##  return(object)
##})
##

