graphConverter <- function(graph,
                          srcSpecies,
                          destSpecies,
                          srcIDType="UNIPROT",
                          destIDType="UNIPROT",
                          keepMultGeneMatches=FALSE,
                          keepMultProtMatches=FALSE){
  ##argument checking.
  if(!is(graph, "graphNEL")){stop("parameter graph needs to be a graphNEL.")}
  
  ##Format as a matrix for easy accessibility
  mat = as(graph, "matrix")
  if(dim(mat)[1]==dim(mat)[2]){
    ##ReMap the col IDs        
    mcnames = idConverter(colnames(mat),
      srcSpecies,
      destSpecies,
      srcIDType,
      destIDType,
      keepMultGeneMatches,
      keepMultProtMatches,
      keepMultDestIDMatches = FALSE)#you always need to limit what is returned if things are not an EG
    ##If the names have dupicates, then we want to remove those too.
    ##This means that entries with mutiple results will only ever return one hit per.
    ##This filtering should only really matter when we are mapping intra-species.
    dupNameIndex = duplicated(names(mcnames))
    mcnames = mcnames[!dupNameIndex]    
    
    if( length(mcnames)==0 ){warning("There are no matches in the destination species."); return(NA)}
    ##Do some indexing to adjust for the fact that not everything will map over
    cnames = colnames(mat)
    names(cnames) = colnames(mat)
    cindex = cnames[names(mcnames)]
    
    ##select out the data from the matrix that we want.
    if(length(cindex)>1){
      mat = mat[cindex,cindex]
    }else{warning("There are not enough salvageable nodes remaining to make a network (map-able nodes must be > 1)."); return(NA)}
    ##then finally reassign to the "dest species" names. #fails for 4,5
    colnames(mat) = rownames(mat) = mcnames[cindex]
    
  }else{warning("graph matrix is not symmetrical."); return(NA)}
  
  ##Reformat as a graphNEL
  as(mat,"graphNEL")
}



#inciMat() ##gives the matrix of the hypergraph
#mention that if they pass an incidence matrix

hyperGraphConverter <- function(graph,
                                srcSpecies,
                                destSpecies,
                                srcIDType="UNIPROT",
                                destIDType="UNIPROT",
                                mapCols=FALSE,
                                keepMultGeneMatches=FALSE,
                                keepMultProtMatches=FALSE){
  ##argument checking:
  if( !is(graph,"Hypergraph") && !is(graph,"matrix") ){
    stop("parameter graph needs to be a Hypergraph or an incidence matrix.")}
  
  ##graph could be a matrix or a Hypergraph.  If its a Hypergraph, then we want to convert it for now
  if(is(graph, "Hypergraph")){
    mat = inciMat(graph)
  }else{
    mat = graph
  }

  ##Its ALWAYS going to be the rows that are the nodes when you call inciMat()    
  ##ReMap the row IDs        
  mrnames = idConverter(rownames(mat),
    srcSpecies,
    destSpecies,
    srcIDType,
    destIDType,
    keepMultGeneMatches,
                          keepMultProtMatches,
    keepMultDestIDMatches = FALSE)#you always need to limit what is returned if things are not an EG
  if( length(mrnames)==0 ){stop("There are no matches in the destination species (matching on the rownames).")}
  
  if(mapCols){
    mcnames = idConverter(colnames(mat),
      srcSpecies,
      destSpecies,
                              srcIDType,
      destIDType,
      keepMultGeneMatches,
      keepMultProtMatches,
      keepMultDestIDMatches = FALSE)#you always need to limit what is returned if things are not an EG
    if( length(mcnames)==0 ){stop("There are no matches in the destination species (matching on the rownames).")}
    }
  ##Do some indexing to adjust for the fact that not everything will map over
  rnames = rownames(mat)
  names(rnames) = rownames(mat)
  rindex = rnames[names(mrnames)]
  
  if(mapCols){
    cnames = colnames(mat)
    names(cnames) = colnames(mat)
    cindex = cnames[names(mcnames)]
  }
  ##select out the data from the matrix that we want.
  mat = mat[rindex,]
  if(mapCols){
    mat = mat[,cindex]
    }    
  ##then finally reassign to the "dest species" names.
  rownames(mat) = mrnames[rindex]
  if(mapCols){
    colnames(mat) = mcnames[cindex]
  }
  ##return the correct format
  if(is(graph, "Hypergraph")){
    ##Reformat as a Hypergraph and return
    return(revInciMat(mat)) 
  }else{
    ##return incidence matrix
    return(mat)
  }
}

