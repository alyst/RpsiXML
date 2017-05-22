##----------------------------------------------------------------------------##
## Generic functions for the RpsiXML package
##
## Author: Jitao David Zhang <jitao_david.zhang@roche.com>
##----------------------------------------------------------------------------##

##----------------------------------------##
## RpsiXML object generics
##----------------------------------------##

setGeneric("sourceDb",function(x) standardGeneric("sourceDb"))
setGeneric("sourceDb<-",function(x, value) standardGeneric("sourceDb<-"))
setGeneric("sourceId",function(x) standardGeneric("sourceId"))
setGeneric("sourceId<-",function(x, value) standardGeneric("sourceId<-"))

setGeneric("taxId", function(x) standardGeneric("taxId"))
setGeneric("taxId<-", function(x, value) standardGeneric("taxId<-"))
setGeneric("organismName", function(x) standardGeneric("organismName"))
setGeneric("organismName<-", function(x, value) standardGeneric("organismName<-"))

setGeneric("interactions", function(x) standardGeneric("interactions"))
setGeneric("interactions<-", function(x, value) standardGeneric("interactions<-"))
setGeneric("interactors", function(x) standardGeneric("interactors"))
setGeneric("interactors<-", function(x, value) standardGeneric("interactors<-"))
setGeneric("numInteractors", function(x) standardGeneric("numInteractors"))
setGeneric("numInteractions", function(x) standardGeneric("numInteractions"))

setGeneric("releaseDate", function(x) standardGeneric("releaseDate"))
setGeneric("releaseDate<-", function(x, value) standardGeneric("releaseDate<-"))

setGeneric("uniqueIdentifierSymbol", function(x) standardGeneric("uniqueIdentifierSymbol"))
           
setGeneric("xrefDbs", function(x,...) standardGeneric("xrefDbs"))
setGeneric("translateID", function(r, ...) standardGeneric("translateID"))


setGeneric("complexes", function(x) standardGeneric("complexes"))
setGeneric("participants", function(x) standardGeneric("participants"))
setGeneric("participants<-", function(x,value) standardGeneric("participants<-"))
setGeneric("shortLabel", function(x) standardGeneric("shortLabel"))
setGeneric("shortLabel<-", function(x,value) standardGeneric("shortLabel<-"))
setGeneric("fullName", function(x) standardGeneric("fullName"))

setGeneric("attributesList", function(x) standardGeneric("attributesList"))
setGeneric("attributesList<-", function(x,value) standardGeneric("attributesList<-"))

setGeneric("complexName", function(x,...) standardGeneric("complexName"))
setGeneric("interactorRef", function(x) standardGeneric("interactorRef"))
setGeneric("interactorRef<-", function(x,value) standardGeneric("interactorRef<-"))


setGeneric("interactorInfo", function(x) standardGeneric("interactorInfo"))
setGeneric("uniprot", function(x) standardGeneric("uniprot"))

## count methods
## edge/node methods
setGeneric("edgeLabel", function(x,...) standardGeneric("edgeLabel"))
setGeneric("hyperedgeNodes", function(x,...) standardGeneric("hyperedgeNodes"))

## accessors
setGeneric("accession", function(x) standardGeneric("accession"))
setGeneric("bait", function(x,...) standardGeneric("bait"))
setGeneric("prey", function(x,...) standardGeneric("prey"))
setGeneric("pubmedID", function(x,...) standardGeneric("pubmedID"))
setGeneric("confidenceValue", function(x,...) standardGeneric("confidenceValue"))
setGeneric("isNegative", function(x, ...) standardGeneric("isNegative"))
setGeneric("isModeled", function(x, ...) standardGeneric("isModeled"))
setGeneric("isIntraMolecular", function(x, ...) standardGeneric("isIntraMolecular"))

setGeneric("revInciMat", function(x,...) standardGeneric("revInciMat"))

##----------------------------------------##
## XML interfacing methods
##----------------------------------------##
##setGeneric("psimi25NamesAlias", function(iValue, type, typeAc) standardGeneric("psimi25NamesAlias"))
##setGeneric("psimi25NamesType", function(shortLabel, fullName, alias) standardGeneric("psimi25NamesType"))
##setGeneric("psimi25Attribute", function(iValue, name, nameAc) standardGeneric("psimi25Attribute"))
##setGeneric("psimi25DbReferenceType",
##           function(list, db, dbAc, id, secondary, version, refType, refTypeAc)
##           standardGeneric("psimi25DbReferenceType"))
##setGeneric("psimi25XrefType", function(primaryRef, secondaryRef) standardGeneric("psimi25XrefType"))
##setGeneric("psimi25BibrefType", function(xref, attributeList) standardGeneric("psimi25BibrefType"))
##setGeneric("psimi25AvailabilityType", function(iValue, id) standardGeneric("psimi25AvailabilityType"))
##setGeneric("psimi25CvType", function(name, xref) standardGeneric("psimi25CvType"))
##setGeneric("psimi25Source", function(name, xref, bibref, attributeList, release, releaseDate) standardGeneric("psimi25Source"))
##setGeneric("psimi25OpenCvType", function(name, xref, attributeList) standardGeneric("psimi25OpenCvType"))
##setGeneric("psimi25BioSourceType", function(name, cellType, compartment, tissue, ncbiTaxId) standardGeneric("psimi25BioSourceType"))
##setGeneric("psimi25ConfidenceType", function(unit, value) standardGeneric("psimi25ConfidenceType"))
##setGeneric("psimi25InteractorElementType", function(name, xref, attributeList, interactorType, organism, sequence, id)
##           standardGeneric("psimi25InteractorElementType"))
##setGeneric("psimi25ExperimentType", function(name, bibref, xref, hostOrganismList,
##                                             interactionDetectionMethod, participantIdentificationMethod,
##                                             featureDetectionMethod, confidenceList, attributeList, id) standardGeneric("psimi25ExperimentType"))
##setGeneric("psimi25ExperimentList", function(experimentRef, experimentDescription)
##           standardGeneric("psimi25ExperimentList"))
##setGeneric("psimi25ExperimentRefListType", function(object) standardGeneric("psimi25ExperimentRefListType"))
##setGeneric("psimi25InteractionElementType",
##           function(name, xref,
##                    availabilityRef, availability,
##                    experimentList, participantList,
##                    interredInteractionList, interactionType,
##                    modelled, intraMolecular,
##                    negative, confidenceList,
##                    parameterList, attributeList, imexId, id)
##           standardGeneric("psimi25InteractionElementType"))
##
##setGeneric("psimi25CvExperimentRefs", function(cv, experimentRefList) standardGeneric("psimi25CvExperimentRefs"))
##setGeneric("psimi25CvExperimentRefsList", function(object) standardGeneric("psimi25CvExperimentRefsList"))
##setGeneric("psimi25HostOrganism", function(bioSourceType, experimentRefList) standardGeneric("psimi25HostOrganism"))
##setGeneric("psimi25ExperimentInteractor", function(interactorRef, interactor, experimentRefList) standardGeneric("psimi25ExperimentInteractor"))
##setGeneric("psimi25Entry", function(source, availabilityList, experimentList,
##                                    interactorList, interactionList, attributeList)
##           standardGeneric("psimi25Entry"))

##----------------------------------------##
## accessors (internal)
##----------------------------------------##
setGeneric("name", function(object) standardGeneric("name"))
setGeneric("name<-", function(object, value) standardGeneric("name<-"))
setGeneric("nameAc", function(object) standardGeneric("nameAc"))
setGeneric("nameAc<-", function(object, value) standardGeneric("nameAc<-"))
setGeneric("iValue", function(object) standardGeneric("iValue")) ## inherited value (.Data)
setGeneric("iValue<-", function(object, value) standardGeneric("iValue<-")) ## inherited value (.Data)
setGeneric("id", function(object) standardGeneric("id"))
setGeneric("id<-", function(object, value) standardGeneric("id<-"))
setGeneric("xrefs", function(object) standardGeneric("xrefs"))
setGeneric("xrefs<-", function(object,value) standardGeneric("xrefs<-"))
setGeneric("attributeList", function(object) standardGeneric("attributeList"))
setGeneric("attributeList<-", function(object,value) standardGeneric("attributeList<-"))
setGeneric("cellType", function(object) standardGeneric("cellType"))
setGeneric("cellType<-", function(object, value) standardGeneric("cellType<-"))
setGeneric("compartment", function(object) standardGeneric("compartment"))
setGeneric("compartment<-", function(object, value) standardGeneric("compartment<-"))
setGeneric("tissue", function(object) standardGeneric("tissue"))
setGeneric("tissue<-", function(object, value) standardGeneric("tissue<-"))
setGeneric("ncbiTaxId", function(object) standardGeneric("ncbiTaxId"))
setGeneric("ncbiTaxId<-", function(object, value) standardGeneric("ncbiTaxId<-"))
setGeneric("unit", function(object) standardGeneric("unit"))
setGeneric("unit<-", function(object, value) standardGeneric("unit<-"))
setGeneric("value", function(object) standardGeneric("value"))
setGeneric("value<-", function(object, value) standardGeneric("value<-"))
setGeneric("interactorType", function(object) standardGeneric("interactorType"))
setGeneric("interactorType<-", function(object, value) standardGeneric("interactorType<-"))
setGeneric("organism", function(object) standardGeneric("organism"))
setGeneric("organism<-", function(object, value) standardGeneric("organism<-"))
setGeneric("Sequence", function(object) standardGeneric("Sequence"))
setGeneric("Sequence<-", function(object, value) standardGeneric("Sequence<-"))

setGeneric("experimentList", function(object) standardGeneric("experimentList"))
setGeneric("experimentList<-", function(object, value) standardGeneric("experimentList<-"))

setGeneric("participantList", function(object) standardGeneric("participantList"))
setGeneric("participantList<-", function(object, value) standardGeneric("participantList<-"))
setGeneric("inferredInteractionList", function(object) standardGeneric("inferredInteractionList"))
setGeneric("inferredInteractionList<-", function(object,value) standardGeneric("inferredInteractionList<-"))
setGeneric("interactionType", function(object) standardGeneric("interactionType"))
setGeneric("interactionType<-", function(object, value) standardGeneric("interactionType<-"))
setGeneric("interactionDetectionMethod", function(object) standardGeneric("interactionDetectionMethod"))
setGeneric("interactionDetectionMethod<-", function(object, value) standardGeneric("interactionDetectionMethod<-"))
setGeneric("modelled", function(object) standardGeneric("modelled"))
setGeneric("modelled<-", function(object,value) standardGeneric("modelled<-"))
setGeneric("intraMolecular", function(object) standardGeneric("intraMolecular"))
setGeneric("intraMolecular<-", function(object,value) standardGeneric("intraMolecular<-"))
setGeneric("negative", function(object) standardGeneric("negative"))
setGeneric("negative<-", function(object,value) standardGeneric("negative<-"))
setGeneric("confidenceList", function(object) standardGeneric("confidenceList"))
setGeneric("confidenceList<-", function(object,value) standardGeneric("confidenceList<-"))
setGeneric("parameterList", function(object) standardGeneric("parameterList"))
setGeneric("parameterList<-", function(object,value) standardGeneric("parameterList<-"))
setGeneric("imexId", function(object) standardGeneric("imexId"))
setGeneric("imexId<-", function(object,value) standardGeneric("imexId<-"))

## pariticpant type

setGeneric("interactor", function(object) standardGeneric("interactor"))
setGeneric("interactor<-", function(object,value) standardGeneric("interactor<-"))
setGeneric("interactionRef", function(object) standardGeneric("interactionRef"))
setGeneric("interactionRef<-", function(object,value) standardGeneric("interactionRef<-"))
setGeneric("participantIdentificationMethodList", function(object) standardGeneric("participantIdentificationMethodList"))
setGeneric("participantIdentificationMethodList<-", function(object,value) standardGeneric("participantIdentificationMethodList<-"))
setGeneric("biologicalRole", function(object) standardGeneric("biologicalRole"))
setGeneric("biologicalRole<-", function(object,value) standardGeneric("biologicalRole<-"))
setGeneric("experimentalRoleList", function(object) standardGeneric("experimentalRoleList"))
setGeneric("experimentalRoleList<-", function(object,value) standardGeneric("experimentalRoleList<-"))
setGeneric("experimentalPreparationList", function(object) standardGeneric("experimentalPreparationList"))
setGeneric("experimentalPreparationList<-", function(object,value) standardGeneric("experimentalPreparationList<-"))
setGeneric("experimentalInteractorList", function(object) standardGeneric("experimentalInteractorList"))
setGeneric("experimentalInteractorList<-", function(object,value) standardGeneric("experimentalInteractorList<-"))
setGeneric("featureList", function(object) standardGeneric("featureList"))
setGeneric("featureList<-", function(object,value) standardGeneric("featureList<-"))
setGeneric("hostOrganismList", function(object) standardGeneric("hostOrganismList"))
setGeneric("hostOrganismList<-", function(object,value) standardGeneric("hostOrganismList<-"))
setGeneric("confidenceList", function(object) standardGeneric("confidenceList"))
setGeneric("confidenceList<-", function(object,value) standardGeneric("confidenceList<-"))
setGeneric("parameterList", function(object) standardGeneric("parameterList"))
setGeneric("parameterList<-", function(object,value) standardGeneric("parameterList<-"))

setGeneric("level", function(object) standardGeneric("level"))
setGeneric("level<-", function(object,value) standardGeneric("level<-"))
setGeneric("version", function(object) standardGeneric("version"))
setGeneric("version<-", function(object,value) standardGeneric("version<-"))
setGeneric("minorVersion", function(object) standardGeneric("minorVersion"))
setGeneric("minorVersion<-", function(object,value) standardGeneric("minorVersion<-"))

setGeneric("attrInfo", function(object,...) standardGeneric("attrInfo"))


##------------------------------------------------------------##
## deprecated
##------------------------------------------------------------##
##setGeneric("typedList", function(..., type) standardGeneric("typedList"))
##setGeneric("ttapply", function(X, INDEX, FUN, ..., simplify) standardGeneric("ttapply"))
##setGeneric("tlapply", function(X, FUN, ... ) standardGeneric("tlapply"))
##

##setGeneric("psimi25AttributeListType", function(list) standardGeneric("psimi25AttributeListType"))
