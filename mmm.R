
#===============================================================================================================================

trim_ <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

#===============================================================================================================================


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
  
  rm.allcolNA(rm.allrowNA(X))  
  
}                     
#===============================================================================================================================

data.tree_ <- function(data, toplab = NULL, cex = 1, ...){
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  plotrix::sizetree(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, ...)
  
} 

#===============================================================================================================================

study_tree_limited <- function(data, study_col, grp2_col, grp3_col, time_col, study_name = NULL, reset = TRUE,
                       structure = c("simple","typical","complex"), output_studies = FALSE,
                       toplab = NULL, cex = 1)
{
  
  data <- rm.colrowNA(trim_(data))
  
  ss <- substitute(study_col)
  sss <- deparse(ss)
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  } 
  
  data <- data %>%
    dplyr::select({{study_col}},{{grp2_col}},{{grp3_col}},{{time_col}})
  
  if(is.null(study_name)){
    
    struc <- match.arg(structure)
    
    hlist <- data %>%
      dplyr::group_by({{study_col}}) %>%
      dplyr::mutate(grp = stringr::str_c(n_distinct({{grp2_col}}) == 1, 
                                         n_distinct({{grp3_col}}) == 1 )) %>%
      dplyr::ungroup(.) %>%
      dplyr::group_split(grp, .keep = FALSE) 
    
    h <- setNames(rev(hlist), c("g_cons_o_cons","g_cons_o_vary","g_vary_o_cons","g_vary_o_vary"))
    
    res <- Filter(NROW, h)
    
    main_title <- sapply(res, function(i) length(unique(i[[sss]])))
    
    typic <- function(vec) vec[ceiling(length(vec)/2)]
    
    nms <- lapply(res, function(i){
      nr <- sapply(split(i, i[[sss]]), nrow);
      study_type <- if(struc == "typical") {typic(as.numeric(names(table(nr))))
      } else if(struc == "simple") {min(as.numeric(names(table(nr))))
      } else {max(as.numeric(names(table(nr))))};
      names(nr)[nr == study_type][1]
    })
    
    list2plot <- lapply(seq_along(res),function(i) subset(res[[i]], eval(ss) == nms[i]))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL))}
    
    main <- paste(main_title,"studies")
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex)))
    
    if(output_studies) res
    
  } else {
    
    study_name <- trimws(study_name)
    study_names <- unique(data[[sss]])
    
    idx <- study_name %in% study_names 
    
    if(!all(idx)) stop(dQuote(toString(study_name[!idx]))," not found in the studies.", call. = FALSE)
    
    list2plot <- lapply(study_name, function(i) subset(data, eval(ss) == i))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL))}
    
    invisible(lapply(list2plot, data.tree_, toplab, cex))
  }
}

#===============================================================================================================================

study_tree <- function(data, highest_level, ..., highest_level_name = NULL, reset = TRUE,
                 structure = c("simple","typical","complex"), output_highest_level = FALSE,
                 toplab = NULL, cex = 1) 
  {
  
  data <- rm.colrowNA(trim_(data))
  
  dot_cols <- rlang::ensyms(...)
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)
  
  ss <- substitute(highest_level)
  sss <- deparse(ss)
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }
  
  data <- data %>%
    dplyr::select({{highest_level}}, !!! dot_cols)
    
  if(is.null(highest_level_name)){
    
    struc <- match.arg(structure) 
  
 hlist <- data %>%
    dplyr::group_by({{highest_level}}) %>%
    dplyr::mutate(grp = across(all_of(str_cols), ~ n_distinct(.) == 1) %>%
                    purrr::reduce(stringr::str_c, collapse="")) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_split(grp, .keep = FALSE)
 
 res <- Filter(NROW, rev(hlist))
 
 main_title <- sapply(res, function(i) length(unique(i[[sss]])))
 
 typic <- function(vec) vec[ceiling(length(vec)/2)]
 
 nms <- lapply(res, function(i){
   nr <- sapply(split(i, i[[sss]]), nrow);
   study_type <- if(struc == "typical") {typic(as.numeric(names(table(nr))))
   } else if(struc == "simple") {min(as.numeric(names(table(nr))))
   } else {max(as.numeric(names(table(nr))))};
   names(nr)[nr == study_type][1]
 })
 
 list2plot <- lapply(seq_along(res),function(i) subset(res[[i]], eval(ss) == nms[i]))
 
 LL <- length(list2plot)
 
 if(LL > 1L) { par(mfrow = n2mfrow(LL)) }
 
 main <- paste(main_title,"studies")
 
 invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex)))
 
 if(output_highest_level) res
 
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) stop(dQuote(toString(highest_level_name[!idx]))," not found in the ", paste0("'",sss,"'", " column."), call. = FALSE)
    
    list2plot <- lapply(highest_level_name, function(i) subset(data, eval(ss) == i))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL))}
    
    invisible(lapply(list2plot, data.tree_, toplab, cex))
  }
}                        
                        
#===============================================================================================================================
                        
needzzsf <- c("plotrix","tidyverse")    

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)


suppressWarnings(
  suppressMessages({ 
    
    for(i in needzzsf){
      library(i, character.only = TRUE)
    }
  }))

options(dplyr.summarise.inform = FALSE)
