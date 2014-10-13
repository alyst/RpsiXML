library(RpsiXML)
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

## HPRD: schema check passed
hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)
hprdInteractors <- interactors(hprdSet)
hprdInteractorInfo <- interactorInfo(hprdSet)
pubmedID(hprdSet)

## BiOGRID: schema check failed, trivial
gridxml <- file.path(xmlDir, "biogrid_200804_test.xml")
gridSet <- parsePsimi25Interaction(gridxml, BIOGRID.PSIMI25)

## large files
gridLargexml <- file.path(xmlDir, "biogrid_2008_test_large.xml")
gridHumanxml <- file.path(xmlDir, "biogrid_2008_test_human.xml")

if(file.exists(gridLargexml)) {
  Rprof("gridLarge.prof")
  gridLargeSet <- parsePsimi25Interaction(gridLargexml, BIOGRID.PSIMI25, verbose=TRUE)
  Rprof(NULL)
}
if(file.exists(gridHumanxml)){
  Rprof("gridHuman.prof")
  gridHumanSet <- parsePsimi25Interaction(gridHumanxml, BIOGRID.PSIMI25)
  Rprof(NULL)
}

## MINT: schema check passed
mintxml <- file.path(xmlDir, "mint_200711_test.xml")
mintSet <- parsePsimi25Interaction(mintxml, MINT.PSIMI25)

## IntAct: schema check passed
intactxml <- file.path(xmlDir, "intact_2008_test.xml")
intactSet <- parsePsimi25Interaction(intactxml, INTACT.PSIMI25)
intactGraph <- psimi25XML2Graph(intactxml, INTACT.PSIMI25)
intactGraphNew <- translateID(intactGraph,"sourceId")## translate the nodes of the graph to another identifier
intactTestInteraction <- interactions(intactSet)[[1]]
bait(intactTestInteraction)
prey(intactTestInteraction)
participant(intactTestInteraction)
inhibitor(intactTestInteraction)
neutralComponent(intactTestInteraction)
pubmedID(intactTestInteraction)
intactTestInteractor <- interactors(intactSet)[[1]]
availableXrefs(intactTestInteractor)

## which cross references does the set provide?
intactXrefs <- availableXrefs(intactSet)
## which references exist for every interactor?
intactXrefIns <- availableXrefs(intactSet, intersect=TRUE)

intactSetInteractors <- interactors(intactSet)
intactXrefExample <- xref(intactSetInteractors[[1]])
translateID(intactSetInteractors,"intact")
translateID(intactSetInteractors[[1]],"intact")

intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexSet <- parsePsimi25Complex(intactComplexxml, INTACT.PSIMI25)
intactComplexGraph <- psimi25XML2Graph(intactComplexxml, INTACT.PSIMI25, type="complex")
translateID(intactComplexGraph, "intact", "P49432")
translateID(intactComplexGraph, "intact", NA)

complexSample <- complexes(intactComplexSet)[[2]]
complexName(complexSample)
attributesList(complexSample)
members(complexSample)

## DIP: schema check failed, namespace not unique, modified
dipxml <- file.path(xmlDir, "dip_2008_test.xml")
dipSet <- parsePsimi25Interaction(dipxml, DIP.PSIMI25)
dipInteractions <- interactions(dipSet)
dipInteractionTypes <- sapply(dipInteractions, interactionType)
stopifnot(!any(is.na(dipInteractionTypes)))

## MatrixDB
matrixdbxml <- file.path(xmlDir, "matrixdb_20080609.xml")
matrixdbSet <- parsePsimi25Interaction(matrixdbxml, MATRIXDB.PSIMI25)

## MIPS: schema check failed
##mipsxml <- file.path(xmlDir, "mips_2007_test.xml")
##mipsSet <- parsePsimi25Interaction(mipsxml, MIPS.PSIMI25) ## needs implementation

## CORUM:
corumxml <- file.path(xmlDir,"corumPsimi")
if(file.exists(corumxml)) {
corumSet <- parsePsimi25Interaction(corumxml,CORUM.PSIMI25, verbose=FALSE)
}

##############################
## Validating PSI-MI 2.5 Files
##############################
#okFile <- system.file("extdata/psi25files/intact_2008_test.xml",
#                      package="RpsiXML")
#errorFile <- system.file("extdata/psi25files/mips_2007_test.xml",
#                         package="RpsiXML")
#mif25Schema <- system.file("extdata/schemas/MIF25.xsd",package="RpsiXML")
#stopifnot(file.exists(okFile) &&
#          file.exists(errorFile) &&
#          file.exists(mif25Schema))
#
#validatePSIMI25(okFile)
#validatePSIMI25(errorFile)
#validatePSIMI25(errorFile, ignore.stderr=FALSE)

