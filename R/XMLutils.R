##----------------------------------------##
## Attributes Parser
##----------------------------------------##
parseXmlAttribute <- function(node, namespaces) {
  attributeName <- xmlGetAttr(node, name="name", default=as.character(NA))
  attributeNameAc <- xmlGetAttr(node, name="nameAc", default=as.character(NA))
  attributeValue <- xmlValue(node)
  
  attribute <- new("psimi25Attribute",
                   name=attributeName,
                   nameAc=attributeNameAc,
                   attribute=attributeValue)
}

parseXmlAttributeNodeSet <- function(nodes, namespaces) {
  if(length(nodes)==0) {
    return(list())
  } else {
    attributes <- sapply(nodes, parseXmlAttribute, namespaces=namespaces)
    return(attributes)
  }
}

parseXmlAttributesListByPath <- function(doc,
                                         path,
                                         namespaces) {
  attributeNodeSets <- getNodeSet(doc=doc,
                                  path=path,
                                  namespaces=namespaces)
  attributes <- parseXmlAttributeNodeSet(attributeNodeSets)
  return(attributes)
}


