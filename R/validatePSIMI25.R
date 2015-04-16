validatePSIMI25 <- function(file,
                            schema=system.file("extdata/schemas/MIF25.xsd",package="RpsiXML"),
                            ignore.stderr=TRUE) {

  ## backup file name
  ofile <- file
  
  ## if file a URL?
  isURL <- length(grep("^(ftp|http|file)://", file)) > 0

  if(isURL) {
    t <- tempfile()
    download.file(file, destfile=t)
    file <- t
  }

  ## so far the routine depends on the tool (xmllint) coming together with libxml2
  ## later it may be replaced by C routines
  if(.Platform$OS.type == "windows") {
    if(interactive()) {
      cat("Oops! So far only UNIX is supported!")
    }
    return(invisible(0))
  }
  
  val <- system(paste("xmllint --noout --schema", schema, file),ignore.stderr=T)

  if ( val == 0 ) {
    cat(paste(ofile, "valid!\n"))
  } else {
    cat(paste(ofile, "invalid! Error code:",val,"\n"))
    if(!ignore.stderr) {
      system(paste("xmllint --noout --schema", schema, file))
    } else {
      cat(paste("Set ignore.stderr=FALSE to see the errors\n"))
    }
  }
  return(invisible(val))
}

