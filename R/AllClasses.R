##----------------------------------------------------------------------------##
## Class definitions of RpsiXML package
##
## Author: Jitao David Zhang
## Description: Data structures of RpsiXML package
##              See documentation and help pages for details
##              The PSI-MI 2.5 schema: http://psidev.sourceforge.net/mi/rel25/doc/
##----------------------------------------------------------------------------##


##------------------------------------------------------------##
## RpsiXML Objects: essential data structures for interactions and
## graphs
##------------------------------------------------------------##
##----------------------------------------##
## Private and Virtual Classes
##----------------------------------------##

## created helper classes
setClass("sourceDbAndId",
         representation(sourceDb="character",
                        sourceId="character"),
         contains="VIRTUAL")
setClass("organismTaxIdAndName",
         representation(taxId="character",
                        organismName="character"),
         contains="VIRTUAL")
setClass("interactorListBase",
         representation(interactors="list"),
         contains="VIRTUAL")
setClass("attributesListBase",
         representation(attributesList="list"),
         contains="VIRTUAL")
setClass("psimi25Entry",
         representation(releaseDate="character"),
         contains=c("organismTaxIdAndName","interactorListBase", "VIRTUAL"))

setClass("psimi25GraphBase",
         representation(abstract = "pubMedAbst"),
         contains=c("interactorListBase", "VIRTUAL"))

## basic classes defined by PSI-MI 2.5
setClass("psimi25Names",
         representation(shortLabel="character",
                        fullName="character",
                        alias="character"))
setClass("psimi25Xref",
         contains="environment")
setClass("psimi25Cv",
         representation(psimi25Names="psimi25Names",
                        xref="psimi25Xref"))
setClass("psimi25Attribute",
         representation(name="character",
                        nameAc="character",
                        attribute="character"),
         prototype=list(name="name",nameAc="nameAc", attribute="attribute"))
setClass("psimi25Confidence",
         representation(unit="psimi25Cv",
                        value="character"))
setClass("psimi25HostOrganism",
         representation(ncbiTaxId="character",
                        psimi25Names="psimi25Names"
                        ## cellType,
                        ## compartment,
                        ## tissue
                        ))
setClass("psimi25Feature",
         representation(psimi25Names="psimi25Names",
                        xref="psimi25Xref",
                        featureType="psimi25Cv",
                        ## featureDetectionMethod,
                        ## experimentRefList,
                        featureRangeList="list"
                        ),
         contains="attributesListBase")


##----------------------------------------##
## Public Classes
##----------------------------------------##

##------------------------------##
## Interaction entry
##------------------------------##

## psimi25 experiment
setClass("psimi25Experiment",
         representation(interactionType = "character",
                        expPubMed = "character"),
         contains=c("sourceDbAndId")
         )

setClass("psimi25Source",
         representation(label="character",
                        sourceDb = "character",
                        uniprotSymbol = "character"),
         prototype = list(
           label = "Label for the database",
           sourceDb = "source database identifier in the psi-mi 25 file",
           uniprotSymbol = "NA"))

## PSI-MI 25 Interaction Entry, roughly corresponding to to the 'entry' definition in MIF25.xsd
setClass("psimi25InteractionEntry",
         representation(interactions = "list"),
         contains=c("psimi25Entry")
         )

setClass("psimi25Participant",
         representation(interactorRef="character",
                        psimi25Names="psimi25Names",
                        xref="psimi25Xref",
                        biologicalRole="psimi25Cv",
                        participantIdentificationList="list",
                        experimentalRoleList="list",
                        experimentalPrepartionList="list",
                        featureList="list",
                        hostOrganismList="list",
                        confidenceList="list"),
         contains="attributesListBase")

## psimi25 interaction, roughly matching the 'interactionElementType'
setClass("psimi25Interaction",
         representation(interactionType = "character",
                        expPubMed = "character",
                        sourceId = "character",
                        confidenceValue = "character",
                        participant = "character",
                        bait = "character",
                        baitUniProt = "character",
                        prey = "character",
                        preyUniProt = "character",
                        inhibitor = "character",
                        neutralComponent = "character",
                        isNegative = "logical",
                        isModeled = "logical",
                        isIntraMolecular = "logical"),
         contains=c("sourceDbAndId")
         )

## psimi25 interactor, roughly matching the 'interactorElementType'
setClass("psimi25Interactor",
         representation(shortLabel = "character",
                        uniprotId = "character",
                        xref = "environment"),
         contains=c("sourceDbAndId", "organismTaxIdAndName")
         )

##------------------------------##
## Complex entry
##------------------------------##
setClass("psimi25ComplexEntry",
         representation(complexes = "list"),
         contains=c("psimi25Entry")
         )

setClass("psimi25Complex",
         representation(shortLabel = "character",
                        fullName = "character",
                        members = "data.frame",
                        interactorRef="character"
                        ),
         contains=c("sourceDbAndId", "organismTaxIdAndName", "attributesListBase")
         )

##------------------------------##
## Graph: graph, hypergraph, etc
##------------------------------##
setClass("psimi25Graph",
         contains = c("psimi25GraphBase", "graphNEL"))

setClass("psimi25Hypergraph",
         representation(interactors = "list"),
         contains = c("psimi25GraphBase", "Hypergraph"))

#### ATTENTION: new XML interfaces not in use
####------------------------------------------------------------##
#### PSI-MI 25 Elementary parent classes
#### 24 complex types (R name, wcw. with convenient wrapper) in dependence-order
#### 1.  namesType (psimi25NamesType)
#### 2.  attributeListType (psimi25AttributeListType)
#### 3.  availabilityType (psimi25AvailabilityType)
#### 4.  dbReferenceType (psimi25DbRereferenceType)
#### 5.  xrefType (psimi25XrefType)
#### 6.  bibrefType (psimi25BibrefType)
#### 7.  cvType (psimi25CvType, wcw. psimi25CvTypeList)
#### 8.  bioSourceType (psimi25BioSourcType, wcw. psimi25BioSourceTypeList)
#### 9.  openCvType (psimi25OpenCvType, wcw. psimi25OpenCvTypelist)
#### 10.  experimentRefListType
#### 11. confidenceType (psimi25ConfidenceType)
#### 12. confidenceListType (psimi25ConfidenceListType)
#### 13. interactorElementType (psimi25InteractorElementType)
#### 14. experimentType (psimi25ExperimentType)
#### 15. fullNameType (psimi25FullNameType)
#### 16. labelType (psimi25LabelType)
#### 17. parameterType (psimi25ParameterType)
#### 18. intervalType (psimi25IntervalType)
#### 19. positionType (psimi25PositionType)
#### 20. baseLocationType (psimi25BaseLocationType) --> needs testing
#### 21. featureElementType (psimi25FeatureElementType)
#### 22. participantType (psimi25PariticipantType)
#### 23. interactionElementType (psimi25InteractionElementType)
####------------------------------------------------------------##
##
####--------------------##
#### namesType
####--------------------##
##setClass("psimi25NamesAlias",
##         representation(typeAc="character",
##                        type="character"),
##         contains="character",
##         prototype=prototype(as.character(NA),
##                    type=as.character(NA),
##                    typeAc=as.character(NA)))
##
##
#### 'alias' must be list of 'psimi25NamesAlias' object
##setClass("psimi25NamesType",
##         representation(shortLabel="character",
##                        fullName="character",
##                        alias="list"),
##         prototype=prototype(shortLabel=as.character(NA),
##           fullName=as.character(NA)
##         ))
##GC_mockNamesType <- new("psimi25NamesType")
##
####--------------------##
#### attributeListType
####--------------------##
##setClass("psimi25Attribute",
##         representation(name="character", ## required
##                        nameAc="character"),
##         prototype=prototype(name=as.character(NA), nameAc=as.character(NA)),
##         contains="character")
##
##GC_mockAttributeListType <- list()
##
####--------------------##
#### availabilityType
####--------------------##
##setClass("psimi25AvailabilityType",
##         representation(id="integer"), 
##         contains="character")
##
##
##GC_mockAvailabilityType <- new("psimi25AvailabilityType")
##GC_mockAvailabilityTypeList <- list()
##
####--------------------##
#### dbReferenceType
####--------------------##
##setClass("psimi25DbReferenceType",
##         representation(db="character", ## required
##                        dbAc="character",
##                        id="character", ## required
##                        secondary="character",
##                        version="character",
##                        refType="character",
##                        refTypeAc="character"),
##         prototype=prototype(list(),
##           db=as.character(NA),dbAc=as.character(NA),id=as.character(NA),
##           secondary=as.character(NA), version=as.character(NA), refType=as.character(NA),
##           refTypeAc=as.character(NA)),
##         contains="list")
##
##GC_mockDbReferenceType <- new("psimi25DbReferenceType")
##
##
##
####--------------------##
#### xrefType
####--------------------##
##setClass("psimi25XrefType",
##         representation(primaryRef="psimi25DbReferenceType",
##                        secondaryRef="list"),
##         )
##
####--------------------##
#### bibrefType
####--------------------##
##GC_mockXrefType <- new("psimi25XrefType")
##setClass("psimi25BibrefType",
##         representation(xref="psimi25XrefType",
##                        attributeList="list")) ## ONLY one of the xref/attributelist is allowed
##GC_mockBibrefType <- new("psimi25BibrefType")
##
####--------------------##
#### cv and opencv
####--------------------##
###### common combinations
##setClass("psimi25CommonNameRef",
##         representation(name="psimi25NamesType",
##                        xref="psimi25XrefType"),
##         prototype=prototype(name=new("psimi25NamesType"),
##           xref=new("psimi25XrefType"))
##         )
##
##setClass("psimi25CommonNameRefAttr",
##         representation(attributeList="list"),
##         contains="psimi25CommonNameRef"
##         )
##setClass("psimi25CvType",
##         contains="psimi25CommonNameRef")
##
##setClass("psimi25OpenCvType",
##         contains="psimi25CommonNameRefAttr")
##
##GC_mockCvType <- new("psimi25CvType")
##GC_mockOpenCvType <- new("psimi25OpenCvType")
##GC_mockCvTypeList <- list()
##
####--------------------##
#### bioSourceType
####--------------------##
##setClass("psimi25BioSourceType",
##         representation(name="psimi25NamesType",
##                        cellType="psimi25OpenCvType",
##                        compartment="psimi25OpenCvType",
##                        tissue="psimi25OpenCvType",
##                        ncbiTaxId="integer")
##         )
##GC_mockBioSourceType <- new("psimi25BioSourceType")
##
####--------------------##
#### experimentRefListType
####--------------------##
##setClass("psimi25ExperimentRef",
##         contains="integer")
##
##GC_mockExperimentRefListType <- list()
##
####----------------------------------------##
#### confidenceType/confidenceListType
####----------------------------------------##
##setClass("psimi25ConfidenceType",
##         representation(unit="psimi25OpenCvType",
##                        value="character"))
##
##GC_mockConfidenceType <- new("psimi25ConfidenceType")
##
##setClass("psimi25Confidence",
##         representation(experimentRefList="list"),
##         contains="psimi25ConfidenceType")
##GC_mockConfidence <- new("psimi25Confidence")
##
##
##GC_mockConfidenceListType <- list()
##
####--------------------##
#### interactorElementType
####--------------------##
#### Attention: since 'names' is keyword in 'representation' function, the slot is renamed as 'name'
##setClass("psimi25InteractorElementType",
##         representation(interactorType="psimi25CvType",
##                        organism="psimi25BioSourceType",
##                        sequence="character", ## will be replaced by BioString! TODO
##                        id="integer"),
##         contains="psimi25CommonNameRefAttr")
##GC_mockInteractorElementType <- new("psimi25InteractorElementType")
##
##setClass("psimi25InteractorElementTypeList",
##         prototype=list(),
##         contains="list")
##GC_mockInteractorElementTypeList <- new("psimi25InteractorElementTypeList")
##
####--------------------##
#### experimentType
####--------------------##
#### Attention: since 'names' is keyword in 'representation' function, the slot is renamed as 'name'
##setClass("psimi25ExperimentType",
##         representation(bibref="psimi25BibrefType",
##                        hostOrganismList="list",
##                        interactionDetectionMethod="psimi25CvType",
##                        participantIdentificationMethod="psimi25CvType",
##                        featureDetectionMethod="psimi25CvType",
##                        confidenceList="list",
##                        id="integer"),
##         contains="psimi25CommonNameRefAttr"
##         )
##GC_mockExperimentType <- new("psimi25ExperimentType")
##
##
####--------------------##
#### other simple types
####--------------------##
##setClass("psimi25FullNameType", contains="character")
##setClass("psimi25LabelType", contains="character")
##setClass("psimi25PositionType", representation(position="integer"))
##setClass("psimi25ParameterType",
##         representation(term="character", ## required
##                        termAc="character",
##                        unit="character",
##                        unitAc="character",
##                        base="integer",
##                        exponent="integer",
##                        factor="numeric"),
##         contains="numeric")
##GC_mockParameterType <- new("psimi25ParameterType")
##
##
##GC_mockParameterTypeList <- list()
##
##setClass("psimi25IntervalType",
##         representation(begin="integer", ## required
##                        end="integer"), ## required
##         )
####--------------------##
#### baseLocationType
####--------------------##
##setClass("psimi25BaseLocationTypeStartAtom",
##         representation(startStatus="psimi25CvType",
##                        begin="integer",
##                        beginInterval="psimi25IntervalType")) ## only one of begin/beginInterval is allowed
##
##setClass("psimi25BaseLocationTypeEndAtom",
##         representation(endStatus="psimi25CvType",
##                        end="integer",
##                        endInterval="psimi25IntervalType")) ## only one of begin/beginInterval is allowed
##
##setClass("psimi25BaseLocationType",
##         representation(start="list",
##                        end="list",
##                        isLink="logical"),
##         prototype=prototype(isLink=FALSE))
##
####--------------------##                        
#### featureElementType
####--------------------##
#### Attention: since 'names' is keyword in 'representation' function, the slot is renamed as 'name'
##setClass("psimi25FeatureRangeList",
##          contains="list",
##         prototype=list())
##
##setClass("psimi25FeatureElementType",
##         representation(name="psimi25NamesType",
##                        xref="psimi25XrefType",
##                        featureType="psimi25CvType",
##                        featureDetectionMethod="psimi25CvType",
##                        experimentRefList="list",
##                        featureRangeList="psimi25FeatureRangeList",
##                        attributesList="list"))
##
##GC_mockFeatureElementType <- new("psimi25FeatureElementType")
##GC_mockFeatureElementTypeList <- list()
##
####--------------------##
#### participantType
####--------------------##
#### Attention: since 'names' is keyword in 'representation' function, the slot is renamed as 'name'
##setClass("psimi25CvExperimentRefs",
##         representation(experimentRefList="list"),
##         contains="psimi25CvType")
##GC_mockCvExperimentRefs <- new("psimi25CvExperimentRefs")
##
##GC_mockCvExperimentRefsList <- list()
##
##setClass("psimi25ExperimentInteractor",
##         representation(interactorRef="integer",
##                        interactor="psimi25InteractorElementTypeList",## only one of interactorRef or interactor
##                        experimentRefList="list"))
##GC_mockExperimentInteractor <- new("psimi25ExperimentInteractor")
##
##GC_mockExperimentInteractorList <- list()
##
##setClass("psimi25HostOrganism",
##         representation(experimentRefList="list"),
##         contains="psimi25BioSourceType")
##GC_mockHostOrganism <- new("psimi25HostOrganism")
##
##
##GC_mockHostOrganismList <- list()
##
##setClass("psimi25ParticipantType",
##         representation(interactorRef="integer",
##                        interactor="psimi25InteractorElementType",
##                        interactionRef="integer", ## ONLY one of the interactorRef/interactor/interactionRef,
##                        participantIdentificationMethodList="list",
##                        biologicalRole="psimi25CvType",
##                        experimentalRoleList="list",
##                        experimentalPreparationList="list",
##                        experimentalInteractorList="list",
##                        featureList="list",
##                        hostOrganismList="list",
##                        confidenceList="list",
##                        parameterList="list",
##                        id="integer"),
##         contains="psimi25CommonNameRefAttr")
####--------------------##
#### interactionElementType
####--------------------##
#### Attention: since 'names' is keyword in 'representation' function, the slot is renamed as 'name'
##setClass("psimi25ExperimentList",
##         representation(experimentRef="integer",
##                        experimentDescription="psimi25ExperimentType" ## ONLY one of the experimentRef/experimentDescription is needed
##                        ))
##GC_mockExperimentList <- new("psimi25ExperimentList")
#### TODO: not to be confused with the element of entry
##
##setClass("psimi25InferredInteractionParticipantAtom",
##         representation(participantRef="integer",
##                        participantFeatureRef="integer") ## ONLY one of the participantRef/participantFeatureRef is needed
##         )
##
##setClass("psimi25InferredInteraction",
##         representation(participant="list",
##                        experimentRefList="list")
##         )
##GC_mockInferredInteraction <- new("psimi25InferredInteraction")
##
##
##GC_mockInferredInteractionList <- list()
##
##setClass("psimi25InteractionElementType",
##         representation(availabilityRef="integer",
##                        availability="psimi25AvailabilityType", ## ONLY one of the availabilityRef/availability is needed
##                        experimentList="psimi25ExperimentList",
##                        participantList="list",
##                        inferredInteractionList="list",
##                        interactionType="list",
##                        modelled="logical",
##                        intraMolecular="logical",
##                        negative="logical",
##                        confidenceList="list",
##                        parameterList="list",
##                        imexId="character",
##                        id="integer"),
##         prototype=prototype(negative=FALSE, intraMolecular=FALSE),
##         contains="psimi25CommonNameRefAttr"
##         )
##GC_mockInteractionElementType <- new("psimi25InteractionElementType")
##
####----------------------------------------------------------------------------##
#### PSI-MI XML Elements
####----------------------------------------------------------------------------##
##setClass("psimi25Source", ## ATTENTION: old API is affected!
##         representation(bibref="psimi25BibrefType",
##                        release="character",
##                        releaseDate="character"),
##         contains="psimi25CommonNameRefAttr"
##         )
##
##setClass("psimi25Entry",
##         representation(source="psimi25Source",
##                        availabilityList="list",
##                        experimentList="list",
##                        interactorList="list",
##                        interactionList="list",
##                        attributeList="list")
##         )
##
##setClass("psimi25EntrySet",
##         representation(level="integer",
##                        version="integer",
##                        minorVersion="integer"),
##         contains="list")
##
##
####------------------------------------------------------------##
#### obsolete classes
####------------------------------------------------------------##
##
####validityTypedList <- function(object) {
####  if(length(object)==0) {
####    ##    warning("Empty object! It is not possible to determine type sanity\n")
####    return(TRUE)
####  }
####  
####  classes <- lapply(object, class)
####  uClass <- unique(classes)
####
####  if(length(uClass)==1) {
####    val <- TRUE
####    expType <- object@type
####    if(length(expType)>1 || !is.na(expType)) {
####      if(!all(uClass[[1]] == expType)) {
####        val <- paste(expType, " Type(s) expected, but got ", uClass[[1]], " type (s)\n", sep="")
####      }
####    }
####  } else {
####    val <- "Hetero. types detected\n"
####  }
####  return(val)
####}
####
####setClass("typedList",
####         representation(type="character"),
####         prototype=prototype(type=as.character(NA)),
####         contains="list",
####         validity=validityTypedList)
####
##
####setClass("psimi25NamesAliasList",
####          contains ="list",
####          prototype=list())
####
##
####setClass("psimi25AttributeListType",
####         prototype=list(),
####         contains="list")
####GC_mockAttributeListType <- new("psimi25AttributeListType")
####
##
####setClass("psimi25AvailabilityTypeList",
####         prototype=list(),
####         contains="list")
####
##
####GC_mockAvailabilityTypeList <- new("psimi25AvailabilityTypeList")
##
####setClass("psimi25DbReferenceTypeList",
####         prototype=list(),
####         contains="list")
####
##
####setClass("psimi25BibrefType",
####         representation(xref="psimi25XrefType",
####                        attributeList="list"),
####         validity=function(object) {
####           attLength <- length(object@attributeList)
####           newXref <- !identical(GC_mockXrefType, length(object@xref@primaryRef) == 1)
####           if (attLength && newXref)
####             return("Choice between xref and attributeList in psimi25XrefType!\n")
####           return(TRUE)
####         }## only ONE of xref/attributeList
####         )
##
##
####setClass("psimi25CvTypeList",
####         prototype=list(),
####         contains="list")
####
##
##
####setClass("psimi25OpenCvTypeList",
####         prototype=list(),
####        contains="list")
##
##
####GC_mockOpenCvType <- new("psimi25OpenCvType")
####GC_mockCvTypeList <- new("psimi25CvTypeList")
##
####setClass("psimi25BioSourceTypeList",
####         prototype=list(),
####         contains="list")
####
##
####setClass("psimi25ExperimentRefListType",
####         prototype=list(),
####         contains="list")
####GC_mockExperimentRefListType <- new("psimi25ExperimentRefListType")
##
####setClass("psimi25ConfidenceListType",
####         prototype=list(),
####         contains="list")
####GC_mockConfidenceListType <- new("psimi25ConfidenceListType")
##
####setClass("psimi25ExperimentTypeList",
####         prototype=list(),
####         contains="list")
####
##
##
####setClass("psimi25ParameterTypeList",
####         contains="list",
####         prototype=list())
####GC_mockParameterTypeList <- new("psimi25ParameterTypeList")
####
##
####setClass("psimi25BaseLocationTypeStart",
####         contains="list",
####         prototype=list())
##
####setClass("psimi25BaseLocationTypeEnd",
####         contains="list",
####         prototype=list())
####
##
####setClass("psimi25FeatureElementTypeList",
####         prototype=list(),
####         contains="list")
####GC_mockFeatureElementTypeList <- new("psimi25FeatureElementTypeList")
##
##
####setClass("psimi25CvExperimentRefsList",
####         prototype=list(),
####         contains="list")
####GC_mockCvExperimentRefsList <- new("psimi25CvExperimentRefsList")
##
####setClass("psimi25ExperimentInteractorList",
####         contains="list",
####         prototype=list())
####GC_mockExperimentInteractorList <- new("psimi25ExperimentInteractorList")
##
####setClass("psimi25HostOrganismList",
####         contains="list",
####         prototype=list())
####GC_mockHostOrganismList <- new("psimi25HostOrganismList")
##
####setClass("psimi25ParticipantTypeList",
####          prototype=list(),
####         contains="list")
##
####setClass("psimi25InferredInteractionParticipant",
####         contains="list",
####         prototype=list(),
####         )
##
####setClass("psimi25InferredInteractionList",
####         contains="list",
####         prototype=list())
####GC_mockInferredInteractionList <- new("psimi25InferredInteractionList")
##
####setClass("psimi25InteractionElementTypeList",
####         prototype=list(),
####         contains="list")
