
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

#===== Hack to fix the bug in plotrix package regarding its 'showcount' argument (reported the bug to the package author as well)
#===== Credit to: MrFlick (https://stackoverflow.com/a/68570496/7223434) 

sizetree_ <- plotrix::sizetree
environment(sizetree_) <- globalenv()
# This "path" navigates the AST for the function to find the offending line of code
path <- c(8, 3, 5, 4, 2, 3, 2, 3, 2, 3, 8, 3, 5)
orig <- body(sizetree_)[[path]]
orig
## Problem line, no showcount= parameter
# sizetree(nextx, right, top, right + 1, lastcenter = top - xfreq[bar]/2, 
#     showval = showval, stacklabels = stacklabels, firstcall = FALSE, 
#     col = newcol, border = border, base.cex = base.cex)
## fix it up
scall <- orig
scall$showcount <- quote(showcount)
body(sizetree_)[[path]] <- scall         
           
#===============================================================================================================================
           
data.tree_ <- function(data, toplab = NULL, cex = 1, rowcount = FALSE, ...){
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  sizetree_(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, showcount = rowcount, ...)
} 

#===============================================================================================================================

pluralify_ <- function (x, keep.original = FALSE, 
                         irregular = lexicon::pos_df_irregular_nouns) {
  
  stopifnot(is.data.frame(irregular))
  
  hits <- match(tolower(x), tolower(irregular[[1]]))
  
  ends <- "(sh?|x|z|ch)$"
  plural_ <- ifelse(grepl(ends, x), "es", "s")
  out <- gsub("ys$", "ies", paste0(x, plural_))
  out[which(!is.na(hits))] <- irregular[[2]][hits[which(!is.na(hits))]]
  
  c(if (keep.original) {
    x
  }, out)
  
}           
           
#===============================================================================================================================

meta_tree <- function(data, highest_level, ..., highest_level_name = NULL, reset = TRUE,
                      structure = c("simple","typical","complex"), output_highest_level = FALSE,
                      toplab = NULL, cex = 1, main = NULL, rowcount = FALSE) 
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
    
    main_no. <- sapply(res, function(i) length(unique(i[[sss]])))
    
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
    
    main <- if(is.null(main)) ifelse(main_no. > 1, pluralify_(sss), sss) else main
    
    main <- paste(main_no., main)
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount)))
    
    if(output_highest_level) res
    
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) stop(dQuote(toString(highest_level_name[!idx]))," not found in the ", paste0("'",sss,"'", " column."), call. = FALSE)
    
    list2plot <- lapply(highest_level_name, function(i) subset(data, eval(ss) == i))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL)) }
    
    invisible(lapply(list2plot, data.tree_, toplab, cex, rowcount))
  }
}                                                
                        
#===============================================================================================================================

cat_overlap <- function(data, study_id, cat_mod){
  
  study_id <- rlang::ensym(study_id)
  cat_mod <- rlang::ensym(cat_mod)
  
  studies_cats <- 
  data %>%
  dplyr::group_by(!!study_id, !!cat_mod) %>%
  dplyr::summarise(effects = n(), .groups = "drop_last")

  cat_names <- paste0(rlang::as_string(cat_mod), c(".x", ".y"))
  
  studies_cats <- 
    studies_cats %>%
    dplyr::inner_join(studies_cats, by = rlang::as_string(study_id)) %>%
    dplyr::group_by(!!!rlang::syms(cat_names)) %>%
    dplyr::summarise(
      studies = n(),
      effects = sum(effects.x),
      .groups = "drop_last") %>% 
    dplyr::mutate(n = paste0(studies, " (", effects, ")") )
  
  studies_cats %>%
    dplyr::select(-studies, -effects) %>%
    tidyr::pivot_wider(names_from = cat_names[2], values_from = n) %>%
    dplyr::rename(`Moderator Category` = cat_names[1])
}                        
                        
#================================================================================================================================
 
pt.curve <- function(X, adjust = 1, compact = NULL, pch = 16, col = 2, cex = .7, seed = 0, reset = TRUE, add = FALSE, na.rm = TRUE, ...) {
  
  if(na.rm) X <- na.omit(X)  
  n.target <- length(X)
  
  d <- density(X, adjust = adjust, n = n.target)
  
  n <- if(!is.null(compact)) { 
    
    auc <- sum(d$y*median(diff(d$x)))/(diff(range(d$x))*max(d$y))
    
    compact*ceiling(n.target/auc)
    
  } else { n.target }
  
  set.seed(seed)
  pts <- data.frame(x = runif(n, min(d$x), max(d$x)), y = runif(n, 0, max(d$y)))
  
  pts <- pts[pts$y < approx(d$x, d$y, xout = pts$x)$y, ]
  
  if(nrow(pts) == 0) stop("Increase the size of sample 'X' OR use 'compact = NULL'.", call. = FALSE)
  
  pts <- pts[sample(seq_len(nrow(pts)), n, replace = TRUE), ]
  
  if(!add){
  
  if(reset) graphics.off()    
  plot(pts, pch = pch, col = col, cex = cex, ...)
    
  } else {
    
  points(pts, pch = pch, col = col, cex = cex, ...)
    
  }
}             

                        
#=================================================================================================================================
                        
needzzsf <- c("plotrix","lexicon","tidyverse")    

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)


suppressWarnings(
  suppressMessages({ 
    
    for(i in needzzsf){
      library(i, character.only = TRUE)
    }
  }))

options(dplyr.summarise.inform = FALSE)
