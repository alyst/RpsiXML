## psimi25 data source

HPRD.PSIMI25 <- new("psimi25Source",label="HPRD",sourceDb="HPRD",uniprotSymbol = "uniprot")

BIOGRID.PSIMI25 <- new("psimi25Source", label="BioGRID", sourceDb="", uniprotSymbol =  "UNIPROTKB")

MINT.PSIMI25 <-  new("psimi25Source", label="MINT", sourceDb="mint", uniprotSymbol =  "uniprotkb")

INTACT.PSIMI25 <- new("psimi25Source", label="IntAct", sourceDb="intact", uniprotSymbol =  "uniprotkb")

DIP.PSIMI25 <- new("psimi25Source", label="DIP", sourceDb="dip",uniprotSymbol = "uniprot knowledge base")

MIPS.PSIMI25 <- new("psimi25Source", label="MIPS", sourceDb="MI")

MATRIXDB.PSIMI25 <- new("psimi25Source", label="MatrixDB", sourceDb="matrixdb",uniprotSymbol="uniprotkb")

CORUM.PSIMI25 <- new("psimi25Source", label="CORUM", sourceDb="corum", uniprotSymbol="uniprotkb")

MPACT.PSIMI25 <- new("psimi25Source", label="mips",sourceDb="MIPS", uniprotSymbol="uniprotkb")
