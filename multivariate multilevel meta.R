rm.allrowNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE])
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}

#===============================================================================================================================
       
rm.allcolNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[, colSums(is.na(i) | i == "") != nrow(i), drop = FALSE])
    
  } else { X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE] }
}

#===============================================================================================================================
           
rm.colrowNA <- function(X){

r <- rm.allrowNA(X)
rm.allcolNA(r)  

}                                      

#===============================================================================================================================

data.tree_ <- function(data, toplab = NULL, cex = 1, auto = FALSE, ...){
  
  if(auto){    
    cats <- sapply(data, Negate(is.numeric))  
    data <- data[cats]
  }
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  plotrix::sizetree(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, ...)
  
} 

#===============================================================================================================================

data_tree <- function(data, vars=c("study","group","outcome","time"),outcome=1, study="Sheen et al"){
  
  data <- rm.colrowNA(trim_(data))
  
  G <- if(!is.null(outcome))substitute(outcome) else NULL
  H <- if(!is.null(study))substitute(study) else NULL
  
  d <- subset(data, study %in% study[outcome > G])
  
  dd <- subset(d,study==H)
  
  data.tree_(dd[vars])
  
  unique(d$study)
}
