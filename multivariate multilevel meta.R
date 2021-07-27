data.tree_ <- function(data, toplab = NULL, cex = 1, auto = FALSE, ...){
  
  if(auto){    
    cats <- sapply(data, Negate(is.numeric))  
    data <- data[cats]
  }
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  plotrix::sizetree(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, ...)
  
} 



data_tree <- function(data, vars=c("study","group","outcome","time"),outcome=1, study="Sheen et al"){
  
  data <- rm.colrowNA(trim_(data))
  
  G <- if(!is.null(outcome))substitute(outcome) else NULL
  H <- if(!is.null(study))substitute(study) else NULL
  
  d <- subset(data, study %in% study[outcome > G])
  
  dd <- subset(d,study==H)
  
  data.tree_(dd[vars])
  
  unique(d$study)
}
