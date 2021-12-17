
#===============================================================================================================================

trim_ <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

#===============================================================================================================================              
              
change_case <- function(X, except = NULL, fun = tolower){
  y <- names(Filter(function(i) is.character(i) | is.factor(i), X[setdiff(names(X), except)]))
  X[y] <- lapply(X[y], fun)
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
           
na_cols <- function(data) names(which(colSums(is.na(data)) > 0))    
           
#===============================================================================================================================

full_clean <- function(data) rm.colrowNA(trim_(data))           
   
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
                                 
#================================================================================================================================   
           
cat_overlap <- function(data, study_id, ...){
  
  study_id <- rlang::ensym(study_id)
  cat_mod <- rlang::ensyms(...)
  cat_nms <- purrr::map_chr(cat_mod, rlang::as_string)
  
  idx <- cat_nms %in% names(data)
  if(!all(idx)) stop(toString(dQuote(cat_nms[!idx]))," not found in the 'data'.", call. = FALSE)
  
  setNames(purrr::map(cat_mod,  ~ {
    
    studies_cats <- 
      data %>%
      dplyr::group_by(!!study_id, !!.x) %>%
      dplyr::summarise(effects = n(), .groups = 'drop')
    nm1 <- rlang::as_string(.x)
    cat_names <- paste0(nm1, c(".x", ".y"))
    
    studies_cats <- 
      studies_cats %>%
      dplyr::inner_join(studies_cats, by = rlang::as_string(study_id)) %>%
      dplyr::group_by(!!!rlang::syms(cat_names)) %>%
      dplyr::summarise(
        studies = n(),
        effects = sum(effects.x), .groups = 'drop') %>% 
      dplyr::mutate(n = paste0(studies, " (", effects, ")") )
    
    out1 <- studies_cats %>%
      dplyr::select(-studies, -effects) %>%        
      tidyr::pivot_wider(names_from = cat_names[2], 
                         values_from = n, names_sort = TRUE) %>%
      dplyr::rename_with(~nm1,  cat_names[1]) %>%
      dplyr::arrange(dplyr::across(tidyselect::all_of(nm1))) %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(nm1), as.character))  
 
   out2 <- out1[-1]
   out2[upper.tri(out2)] <- "-"
   dplyr::bind_cols(out1[1],out2)
    
  }), cat_nms)
}           
      
           
#================================================================================================================================           
           
mod_for_un <- function(data,left,right) crossprod(table(data[c(right,left)])>0) 
                      
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

#===============================================================================================================================
                        
#===== Hack to fix the bug in plotrix package regarding its 'showcount' argument (reported the bug to the package author as well)
#===== Credit to: MrFlick (https://stackoverflow.com/a/68570496/7223434) 

sizetree <- plotrix::sizetree
environment(sizetree) <- globalenv()
# This "path" navigates the AST for the function to find the offending line of code
path <- c(8, 3, 5, 4, 2, 3, 2, 3, 2, 3, 8, 3, 5)
orig <- body(sizetree)[[path]]
orig
## Problem line, no showcount= parameter
# sizetree(nextx, right, top, right + 1, lastcenter = top - xfreq[bar]/2, 
#     showval = showval, stacklabels = stacklabels, firstcall = FALSE, 
#     col = newcol, border = border, base.cex = base.cex)
## fix it up
scall <- orig
scall$showcount <- quote(showcount)
body(sizetree)[[path]] <- scall         
 
#===============================================================================================================================
           
data.tree_ <- function(data, toplab = NULL, cex = 1, rowcount = FALSE, ...){
  
  toplab <- if(is.null(toplab)) names(data) else toplab
  
  sizetree(data, toplab = toplab, stacklabels = FALSE, border = 0, base.cex = cex, showcount = rowcount, ...)
}           
           
#===============================================================================================================================

meta_tree2 <- function(data, highest_level, ..., highest_level_name = NULL, reset = TRUE,
                      structure = c("simple","typical","complex"), output_highest_level = FALSE,
                      toplab = NULL, cex = 1, main = NULL, rowcount = FALSE, main_extra_name = FALSE) 
{
  
  data <- rm.colrowNA(trim_(data))
  
  dot_cols <- rlang::ensyms(...)
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)

  idx <- str_cols %in% names(data)
  if(!all(idx)) stop(toString(dQuote(str_cols[!idx]))," not found in the 'data'.", call. = FALSE)
  
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
                      purrr::reduce(stringr::str_c, collapse = "")) %>%
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
    
    if(main_extra_name) main <- paste0(main, " [",nms,"]")
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount)))
    
    if(output_highest_level) res
    
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) stop(toString(dQuote(highest_level_name[!idx]))," not found in the ", paste0("'",sss,"'", " column."), call. = FALSE)
    
    list2plot <- lapply(highest_level_name, function(i) subset(data, eval(ss) == i))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL)) }
    
    invisible(lapply(list2plot, data.tree_, toplab, cex, rowcount))
  }
}                      

#=======================================================================================================================================================
                        
meta_tree <- function(data, highest_level, ..., highest_level_name = NULL, reset = TRUE,
                      structure = c("simple","typical","complex"), output_highest_level = FALSE,
                      toplab = NULL, cex = 1, main = NULL, rowcount = TRUE, main_extra_name = FALSE) 
{
  
  data <- rm.colrowNA(trim_(data)) %>%
    mutate(row_id = row_number())
  
  dot_cols <- rlang::ensyms(...)
  str_cols <- purrr::map_chr(dot_cols, rlang::as_string)
  
  idx <- str_cols %in% names(data)
  if(!all(idx)) stop(toString(dQuote(str_cols[!idx]))," not found in the 'data'.", call. = FALSE)
  
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
      dplyr::mutate(grp = dplyr::across(tidyselect::all_of(str_cols), ~ {
        tmp <- dplyr::n_distinct(.)
        #dplyr::case_when(tmp  == 1 ~ 1, tmp == n() ~ 2, TRUE ~ 3)
         dplyr::case_when(tmp  == 1 ~ 1, tmp == n() ~ 2, tmp > 1 & tmp < n() ~ 3,  TRUE ~ 4)
      }) %>%
        purrr::reduce(stringr::str_c, collapse = "")) %>%
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
    
    if(main_extra_name) main <- paste0(main, " [",nms,"]")
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount)))
    
    if(output_highest_level) res
    
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) stop(toString(dQuote(highest_level_name[!idx]))," not found in the ", paste0("'",sss,"'", " column."), call. = FALSE)
    
    list2plot <- lapply(highest_level_name, function(i) subset(data, eval(ss) == i))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL)) }
    
    invisible(lapply(list2plot, data.tree_, toplab, cex, rowcount))
  }
}                   
                        
#=====================================================================================================================================================
                        
latent_metareg <- function(fit, formula, group.id = c("first","second"), 
                           tol = 1e-12, std = FALSE) 
{

  has_G <- fit$withG
  has_H <- fit$withH
  group.id <- match.arg(group.id)
  
  if(!has_G) stop("Model must contain at least one correlated random effect term.", call. = FALSE)
  if(!all(fit$struct %in% c("UN", "GEN"))) stop("Model must use 'UN' or 'GEN' structures.", call. = FALSE)
  if(group.id != "first" & !has_H) stop("No second grouping variable detected.", call. = FALSE)
  
  C <- if(group.id == "first") fit$G else fit$H
  R <- stats::cov2cor(C)
  vars <- colnames(C)
  formula <- as.character(formula)
  out_loc <- which(vars == formula[2])
  pred.names <- trimws(strsplit(formula[3], "+", fixed = TRUE)[[1]])
  pred_loc <- which(vars %in% pred.names)
  if (tol < 1e-12) warning("Do not reduce tol, solution may not be numerically stable.")
 
  if (Matrix::det(R) < tol) {
      warning("Singular covariance matrix detected.", call. = FALSE) 
  }

  res <- array(dim = length(pred_loc))
  
  if (std) {
    res <- base::solve(R[pred_loc, pred_loc], R[out_loc, pred_loc], tol = tol)
  } else {
    res <- base::solve(C[pred_loc, pred_loc], C[out_loc, pred_loc], tol = tol)
  }

  names(res) <- if(length(pred.names) > 1) colnames(R[pred_loc, pred_loc]) else pred.names
  return(res)
}                        
                        
#================================================================================================================================================
                        
meta_struct <- function(..., raw = FALSE, row_id = FALSE, missing = FALSE, seed = NULL, comparison = TRUE){
  
  dots <- rlang::list2(...) 
  inp_ls <- purrr::map(dots, ~ purrr::map(.x, seq_len)) %>% transpose %>% 
    purrr::map(purrr::compact)
  
  high_name <- names(dots)[1]
  
  data_ls <- purrr::map(inp_ls, 
                        ~ tidyr::expand_grid(!!! .x,
                                             comparison = c("control","treatment")))
  
  out <- if(!raw) { purrr::map2(data_ls, inp_ls, ~ .x %>% 
                                  dplyr::group_by(!!! rlang::syms(names(.y))) %>%
                                  dplyr::summarise(comparison = stringr::str_c(sort(comparison, decreasing = TRUE), 
                                                                               collapse = ' vs '), .groups = 'drop'))
  } else { data_ls }
  
  res <- dplyr::bind_rows(out)
  
  res <- if(comparison) res else dplyr::select(res, -comparison)
  
  res <- if(row_id) dplyr::bind_cols(res, row_id = seq_len(nrow(res))) else res
  
  if(missing){
    
    set.seed(seed)
    
    res %>% dplyr::group_by(!!!rlang::syms(high_name)) %>% 
      dplyr::sample_n(sample(dplyr::n(), 1)) %>% 
      dplyr::ungroup()
    
  } else {
    
   as.data.frame(res)
  }
} 
      
#===============================================================================================================================================
                     
dcohen <- function(x, pop.d = 0, n1, n2 = NA, log = FALSE){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- pop.d*sqrt(N)
  
  suppressWarnings(dt(x*sqrt(N), df, ncp, log = log)*sqrt(N))
}

#=======================================================================================================================================


qcohen <- function(p, pop.d = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  q <- Vectorize(function(p, pop.d, n1, n2, lower.tail, log.p){
    
    N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    ncp <- pop.d*sqrt(N)
    
    suppressWarnings(qt(p, df, ncp, lower.tail = lower.tail, log.p = log.p)/sqrt(N))
  })
  q(p = p, pop.d = pop.d, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}


#=======================================================================================================================================


pcohen <- function(q, pop.d = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  p <- Vectorize(function(q, pop.d, n1, n2, lower.tail, log.p){
    
    N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    ncp <- pop.d*sqrt(N)
    
    suppressWarnings(pt(q*sqrt(N), df, ncp, lower.tail = lower.tail, log.p = log.p))
  })
  p(q = q, pop.d = pop.d, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}


#=======================================================================================================================================

rcohen <- function(n, pop.d = 0, n1, n2 = NA){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- pop.d*sqrt(N)
  
  suppressWarnings(rt(n, df, ncp)/sqrt(N))
}
   
#===============================================================================================================================================                        
                        
pca_rma <- function(fit, digits = 4, cum = TRUE){
  
  if(!inherits(fit, "rma.mv")) stop("Model must be 'rma.mv'.", call. = FALSE) 
  if(is.null(fit$random)) stop("Model must contain at least one random effect term.", call. = FALSE)
  
  has_G <- fit$withG
  has_H <- fit$withH
  has_S <- fit$withS
  
  if(has_S){
   S <- fit$sigma2
   S <- if(length(S) < 2) diag(S, 1, 1) else diag(S)
   colnames(S) <- rownames(S) <- fit$s.names
   sds <- setNames(svd(chol(S))$d, colnames(S))
   pov_S <- round(sds^2 / sum(sds^2), digits = digits)
   pov_S <- if(cum) cumsum(pov_S) else pov_S
   S.id.name <- "Ave_ESs (Intercepts)"
  }
  
  if(has_G){
   G <- fit$G
   sds <- setNames(svd(chol(G))$d, colnames(G))
   pov_G <- round(sds^2 / sum(sds^2), digits = digits)
   pov_G <- if(cum) cumsum(pov_G) else pov_G
   G.id.name <- tail(fit$g.names, 1)
  }
  
  if(has_H){
    H <- fit$H
    sds <- setNames(svd(chol(H))$d, colnames(H))
    pov_H <- round(sds^2 / sum(sds^2), digits = digits)
    pov_H <- if(cum) cumsum(pov_H) else pov_H
    H.id.name <- tail(fit$h.names, 1)
  }
  
    if(has_S & !has_G)
    setNames(list(pov_S), paste("proprtion of variance:",dQuote(S.id.name)))
  else
    if(has_G & !has_H & !has_S)
    setNames(list(pov_G), paste("proprtion of variance:",dQuote(G.id.name))) 
  else
    if(has_G & has_H & !has_S)
    setNames(list(pov_G, pov_H), paste("proprtion of variance:",dQuote(c(G.id.name, H.id.name))))
  else
    if(has_S & has_G)
    setNames(list(pov_S, pov_G), paste("proprtion of variance:",dQuote(c(S.id.name, G.id.name))))
  else 
    setNames(list(pov_S, pov_G, pov_H), paste("proprtion of variance:",dQuote(c(S.id.name, G.id.name, H.id.name))))
}
   
#================================================================================================================================================
    
prob_treat <- function(a = NULL, legend = "topright"){
  
  from <- -4
  to <- 4
  b <- 0
  c <- 0 
  
  no_a <- is.null(a)
  
  a <- if(no_a) c(-1e9,-.12, -.35, -.8,0) else a
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  
  I <- eq(a, b, c)
  
  a = I[[1]] ; b = I[[2]] ; c = I[[3]]
  
  loop <- seq_along(a)
  
  h <- list()
  graphics.off()  
  for(i in loop){
    p <- function(x) c[i] + ((1 - c[i])/(1 + exp(-a[i]*(x - b[i]))))  
    h[[i]] <- curve(p, from, to, add = i!= 1, n = 1e3, las = 1, ylim = 0:1,
                    font.lab = 2, xlab = "Pre-treatment (Pre-test)", xaxt = "n",
                    type = "n", ylab = "Prob. of Treatment", mgp = c(2, .5, 0), tck = -.015)
  }

  if(!no_a) {u <- par('usr')
     abline(h = c(0, 1, .5), col = 8, lty = c(3, 3, 2))}
  
  axis(1, at = 0, labels = "Mean (cut-off)", mgp = c(2, .2, 0), tck = -.015)
  axis(2, at = .5, col = 2, col.axis = 2, mgp = c(2, .5, 0), tck = -.015, las = 1)
  
  for(i in loop){
    lines(h[[i]], col = i, lwd = 2)
   # segments(b[i], u[3], b[i], .5, col = i, lty = 3)  
    points(b[i], .5, pch = 21, col = i, font = 2, cex = 1.5, bg = "cyan", lwd = 2)
  }
  
 if(no_a) legend(legend, c("RD Design","NEGD: Near RD", "NEGD: Bet. RD & RE", "NEGD: Near RE","RE"), 
         lty = 1, pt.bg = loop, col = c(1,4,3,2,5), cex = .7, pt.cex = .6, 
         bg = 0, box.col = 0, x.intersp = .5)
  box()
}        

#================================================================================================================================================
                              
num2fac <- function(data, ..., caps = FALSE, reverse = FALSE){
  
  var_ls <- purrr::map(rlang::ensyms(..., .named = TRUE), as.character)
  l_varls <- length(var_ls) 
  
  if(l_varls != length(caps)) {
    caps <- rep(caps, l_varls)
  }
  if(l_varls != length(reverse)) {
    reverse <- rep(reverse, l_varls)
  }
  
  let <- ifelse(caps, list(base::LETTERS), list(base::letters))
  names(let) <- names(var_ls)
  FUN <- ifelse(reverse, list(utils::tail), list(utils::head))
  names(FUN) <- names(var_ls)
  
  dplyr::mutate(data,
         purrr::map_dfc(var_ls,
                        ~ factor(FUN[[.x]](let[[.x]],max(data[[.x]]))[data[[.x]]])))
}    
                              
#===============================================================================================================================================
                              
do_context <- function(data, context_vars, group_id){
  
  data <- full_clean(data)
  all_names <- names(data)
  
  id <- grep("id|group|grp|study|study_id|studyid", all_names, value = TRUE, ignore.case = TRUE)
  
  if(!all(group_id %in% all_names) || missing(group_id)) { 
    
    stop(paste(toString(dQuote(group_id)), "not found for 'group_id' in the 'data'.", if(length(id)>0) 
      paste("\nPossibly you meant to use", toString(dQuote(id)), "as 'group_id', no?")), call. = FALSE) 
    
    }
  
  ok <- context_vars %in% all_names
  
  if(!all(ok)) message(paste("\n", toString(dQuote(context_vars[!ok])), "not found in the 'data' thus ignored."))
  
  context_vars <- context_vars[ok] 
  
  dum_vars <- all_names[sapply(data, function(i)is.character(i)|is.factor(i))]
  
  dum_names <- context_vars[context_vars %in% dum_vars]
  
  is_dum <- length(dum_names) > 0
  
  num_names <- context_vars[!(context_vars %in% dum_vars)]
  
  is_num <- length(num_names) > 0
  
  if(is_num){
    data <- data %>%
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(all_of(num_names), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  if(is_dum){
    data <- data %>%
      dummy_cols(select_columns = dum_names) %>% 
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(starts_with(paste0(dum_names, "_")), list(wthn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  return(data)
}
     
                               
#================================================================================================================================================                               
                               
mask <- function(data, what, full = FALSE){
  
  data[] <- lapply(data, function(x) type.convert(as.character(x), as.is = TRUE))
  
  f1 <- function(x) as.numeric(factor(x, levels = unique(x)))
  f2 <- function(x) {
    temp <- substr(x, 1, 1)
    paste0(temp, ave(x, temp, FUN = function(y) match(y, unique(y))))
  }
  cols <- names(data)[sapply(data, is.numeric)]
  num.cols <- cols[cols %in% what]
  cols <- names(data)[sapply(data, is.character)]
  char.cols <- cols[cols %in% what]
  
  if(!full){
    
    if(length(num.cols))  data[num.cols] <- lapply(data[num.cols], f1)
    if(length(char.cols)) data[char.cols] <- lapply(data[char.cols], f2)
    
  }else{
    
    data[what] <- lapply(data[what], f1)
  }
  return(data)
}  

                     
#================================================================================================================================================
                     
is_nestor <- function(data, ...){
  
  dot_cols <- rlang::ensyms(...)
  IDs <- purrr::map_chr(dot_cols, rlang::as_string)
  
  if(length(unique(IDs)) < 2) stop("At least 2 ID variables required.", call. = FALSE)
  
 data <- na.omit(data)
 
is_nest <- function(data, nestor, nested){

f <- function(data, nesto, neste){
groupings <- tapply(data[[nesto]], data[[neste]], function(x) length(unique(x)))
all(groupings == 1L)
}
setNames(sapply(nested, function(i)f(data=data, nestor, i)), nested)
}

setNames(lapply(seq_along(IDs), function(i) is_nest(data, IDs[i], IDs[-i])), IDs)
}                     
    
#================================================================================================================================================
                
mod_levels <- function(data, ...){
  
  dot_cols <- rlang::ensyms(...)
  what <- purrr::map_chr(dot_cols, rlang::as_string)
  
f <- function(dat, wht){
sort(unique(na.omit(unlist(dat[wht]))))
}
setNames(lapply(what, function(i) f(data, i)), what)
}                

#================================================================================================================================================
                
winsor_tukey <- function(data, yi, coef = 3){

yi <- deparse(substitute(yi))

first_third_QR <- fivenum(data[[yi]], na.rm = TRUE)[c(2,4)]

inter_QR <- diff(first_third_QR)

Tukey_min <- first_third_QR[1] - (coef * inter_QR )
Tukey_max <- first_third_QR[2] + (coef * inter_QR )

ind1 <- data[[yi]] < Tukey_min
ind2 <- data[[yi]] > Tukey_max

data[[yi]][ which(ind1 & !is.na(data[[yi]])) ] <- Tukey_min
data[[yi]][ which(ind2 & !is.na(data[[yi]])) ] <- Tukey_max

return(data)
}                                
                
#================================================================================================================================================   
  
interactive_outlier <- function(fit, cook = NULL, st_del_res_z = NULL, 
                                cex_add_point = .5, cex_multi_point = 6,
                                whisker_coef = 2.5, cex_text_outlier = .6,
                                cex_main = .9, parallel = "snow", ncpus = 4, 
                                reestimate = FALSE, save = FALSE, 
                                file_name_cook = "cooks1",
                                file_name_res_z = "rstudent1",
                                view = 1, pos = 2)
{
  
  if(!inherits(fit,c("rma.mv"))) stop("Model is not 'rma.mv()'.", call. = FALSE)
  datziola <- clubSandwich:::getData(fit) %>%
  mutate(obsss = row_number())
  
  
  # Check Hat values (weight-leveraging effects)
  hat <- hatvalues.rma.mv(fit)
  
  if(is.null(cook)){
    
    cook <- cooks.distance.rma.mv(fit, progbar=TRUE,
                                  parallel = parallel, 
                                  ncpus = ncpus, 
                                  reestimate = reestimate)
    
    if(save){
      
      filenm <- paste0(file_name_cook,".rds")
      saveRDS(cook, filenm)
      
      message("\nNote: Check folder '", basename(getwd()),"' for the ", dQuote(filenm)," file.\n") 
    }
  }
  
  if(is.null(st_del_res_z)){
    
    st_del_res_z <- rstudent.rma.mv(fit, progbar=TRUE,
                                    parallel = parallel, 
                                    ncpus = ncpus, 
                                    reestimate = reestimate)$z
    
    if(save){
      
      filenm <- paste0(file_name_res_z,".rds")
      saveRDS(st_del_res_z, filenm)
      
      message("\nNote: Check folder '", basename(getwd()),"' for the ", dQuote(filenm)," file.\n") 
    }
  }
  
  # Make visual size of effects proportional to their cook's distance (estimate influence)
  cex <- cex_add_point+cex_multi_point*sqrt(if(view == 1)cook else hat)
  
  # Plot Leverage against Studentized residuals proportioned on cook's distances
  plot(if(view == 1) hat else cook, st_del_res_z, cex=cex, las=1, mgp=c(1.5,.3,0),
       xlab=if(view == 1) "Leverage (Hat Value)" else "Effect Influence (Cook's Dis.)", 
       ylab="Outlier (Studentized Del. Value)",pch=19,cex.axis = .9,tcl = -.3,
       col = adjustcolor(1, .5))
  title(if(view == 1) "Size of points denote \nestimate-influencing effects\n (Cook's distances)" 
        else "Size of points denote \nleverage effects\n (Hat value)", 
        cex.main = cex_main, line = .3)
  
  outlier_limits <- qnorm(c(.025,.5,.975))
  
  abline(h=outlier_limits, lty=c(3,1,3), lwd=c(1,2,1))
  
  max_hat <- max(mean(range(hat)), boxplot.stats(hat, coef = whisker_coef)$stats[5])
  
  max_cook <- max(mean(range(cook)), boxplot.stats(cook, coef = whisker_coef)$stats[5])
  
  if(view == 1) abline(v = max_hat, col=2) else abline(v = max_cook, col=2)
  
  # To be outlier, an estimate must simultaneously (a) be outlying (per studentized value)
  # (b) have high leverage (hat value), and (c) high model influence (cook's distance)
  
  i <- abs(st_del_res_z) > outlier_limits[3]  
  j <- hat > max_hat
  k <- cook > max_cook
  L <- which(i & j & k)
  
  if(length(L)==0) { 
    
    message("Note: No interactive outlier detected.") 
    
    return(NA)
  }
  
  u <- par()$usr
  
  if(any(st_del_res_z[L]>0)) rect(if(view == 1) max_hat else max_cook, outlier_limits[3], u[2], u[4], col = adjustcolor(2, .2), border = NA)
  if(any(st_del_res_z[L]<0)) rect(if(view == 1) max_hat else max_cook, outlier_limits[1], u[2], u[3], col = adjustcolor(2, .2), border = NA)
  
  # Show which effects meet all the three conditions
  text(if(view == 1) hat[L] else cook[L], st_del_res_z[L], labels = names(L), pos = pos, col = "magenta", cex = cex_text_outlier, xpd = NA)
  points(if(view == 1) hat[L] else cook[L], st_del_res_z[L], cex=cex[L])
  
  LL <- as.numeric(names(L))
  
  removed <- filter(datziola, obsss %in% LL)
  new_data <- filter(datziola, !obsss %in% LL)
  
  return(list(removed = removed, 
              new_data = new_data))
}   

#=================================================================================================================================================
                
maxs <- function(xx, ns){

f <- function(x, n){
  
  len <- length(x)
  x[x == sort(x, partial = len-n+1)[len-n+1]]
}

sapply(seq_len(ns), function(i) f(xx, i))
}


#================================================================================================================================================================

mins <- function(xx, ns){
  
 f <- function(x, n){
  
  sort(x)[n]
}                          
 sapply(seq_len(ns), function(i) f(xx, i))
}                

#================================================================================================================================================================        
pair_list_ <- function(data, groups) {
  
  unique_counts <- table(data[[groups]])
  single_count <- names(unique_counts[unique_counts == 1])
  
  if(length(single_count)) {
    message(sprintf('Note: %s %s with one row, was excluded.', groups, toString(single_count)))
  }
  multiple_count <- names(unique_counts[unique_counts > 1])
  
  lst <- combn(multiple_count, 2, function(x) {
    data[data[[groups]] %in% x, ]
  }, simplify = FALSE)  

return(lst)
}

#=================================================================================================================================================

DEF <- function(observations, groups, nested = TRUE){

  if(nested){
  icc <- ICC::ICCbare(as.factor(groups), observations)
  
  nObs <- tapply(observations, groups, length)
  
  (mean(nObs) - 1) * icc + 1
  } else 1
}

#=================================================================================================================================================

RR_ <- function(observations, groups, 
                ref_group = NULL, conf_level = 0.95, 
                bias_adjust = TRUE, raw = TRUE, 
                DEF = 1) 
{
  
  groups <- factor(groups)
  treat_level <- levels(groups)[(ref_group != levels(groups))]
  level_labs <- c(ref_group, treat_level)
  
  nObs <- table(groups)[level_labs]
  means <- tapply(observations, groups, base::mean)[level_labs]
  variances <- tapply(observations, groups, stats::var)[level_labs]
  
  if (!all(means > 0)) 
    stop("The mean of one or both groups is at the floor of 0.", call. = FALSE)
  if (bias_adjust == TRUE) {
    BC <- log(means) - variances/(2 * nObs * means^2)
    lRR <- as.numeric(BC[2] - BC[1])
  }
  else {
    lRR <- log(means[2]) - log(means[1])
  }
  
  V_lRR <- sum(variances/(nObs * means^2))
  
  percent_dif <- 100*(exp(lRR)-1)
  
    V_lRR <- DEF * V_lRR
    
  CI <- lRR + c(-1, 1) * stats::qnorm(1 - (1 - conf_level)/2)*sqrt(V_lRR)
  
  row_names <- paste(treat_level,ref_group, sep="/")
  
  if (raw) {
    return(data.frame(RR = exp(lRR), percent_dif = percent_dif, 
                      CI_lower = exp(CI[1]), CI_upper = exp(CI[2]), 
                      row.names = row_names))
  }
  else {
    return(data.frame(LRR = lRR, percent_dif = percent_dif, 
                      SE = sqrt(V_lRR), CI_lower = CI[1], CI_upper = CI[2], 
                      row.names = row_names))
  }
}

#=================================================================================================================================================

response_ratio <- function(data, observations, groups, 
                           nested = FALSE, raw = TRUE,
                           bias_adjust = TRUE,
                           conf_level = 0.95){

  data <- rm.colrowNA(trim_(data))
  data_na_rm <- na.omit(data)
  
  if(nrow(data_na_rm) != nrow(data)) message("Note: Rows with missing values were excluded.")
  
  data <- data_na_rm  
  
observations <- deparse(substitute(observations))
groups <- deparse(substitute(groups))

def <- DEF(data[[observations]], data[[groups]], nested = nested)

lstt <- pair_list_(data, groups)

map_dfr(lstt, function(i) 
  map_dfr(1:2, function(j) 
    RR_(i[[observations]], i[[groups]], 
        ref_group = as.vector(unique(i[[groups]]))[j],
        DEF = def, raw = raw, bias_adjust = bias_adjust,
        conf_level = conf_level)))

}

#=================================================================================================================================================
          
fixed_form_rma <- function(fit){ 
   
  a <- fit$formula.yi
  b <- fit$formula.mods
  y <- fit$call$yi
  
   if(!is.null(a) & !is.null(b)) a 
  else if(!is.null(a) & is.null(b)) a
  else if(is.null(a) & !is.null(b)) as.formula(paste(as.character(y), paste(as.character(b), collapse = "")))
}

#=================================================================================================================================================  

plot.efflist <- function (x, selection, rows, cols, graphics = TRUE, 
                          lattice, rug = FALSE, multiline = TRUE, ...) 
{
   lattice <- if (missing(lattice)) 
      list()
   else lattice
   if (!missing(selection)) {
      if (is.character(selection)) 
         selection <- gsub(" ", "", selection)
      pp <- plot(x[[selection]], lattice = lattice, rug = rug, multiline=multiline, ...)
      pp$x.scales$tck=c(1,0)
      pp$y.scales$tck=c(1,0)
      return(pp)
   }
   effects <- gsub(":", "*", names(x))
   neffects <- length(x)
   mfrow <- mfrow(neffects)
   if (missing(rows) || missing(cols)) {
      rows <- mfrow[1]
      cols <- mfrow[2]
   }
   for (i in 1:rows) {
      for (j in 1:cols) {
         if ((i - 1) * cols + j > neffects) 
            break
         more <- !((i - 1) * cols + j == neffects)
         lattice[["array"]] <- list(row = i, col = j, 
                                    nrow = rows, ncol = cols, more = more)
         pp <- plot(x[[(i - 1) * cols + j]], lattice = lattice, rug = rug, multiline = multiline,
                    ...)
         # hack to turn off opposite side tick marks
         pp$x.scales$tck=c(1,0)
         pp$y.scales$tck=c(1,0)
         print(pp)
      }
   }
}
environment(plot.efflist) <- asNamespace("effects")
  
#======================================================================================================================================================  
  
plot_rma_effect <- function(fit, full=FALSE, multiline=TRUE, dots=FALSE,
                            confint = list(style="auto"), x.var, 
                            key.args= list(space="top",cex=.7,cex.title=.8), main=NA,
                            index=NULL, xlab, ylab, z.var, colors, cex, lty, 
                            lwd, ylim, xlim, factor.names, band.transparency, 
                            band.colors, grid=TRUE, axes, lattice, rotx, roty,
                            symbols=list(pch = 19), ticks.x, lines=TRUE, robust=FALSE,
                            cluster, plot=TRUE, ...) 
{
   
   if(!inherits(fit,c("rma.mv","rma","rma.uni"))) stop("Model is not 'rma()' or 'rma.mv()'.", call. = FALSE)
   
   data_. <- eval(fit$call$data)
   form_. <- fixed_form_rma(fit)
   
   fit2 <- nlme::gls(form_., data = data_., method = fit$method, na.action = "na.omit")
   
   fit2$call$model <- form_.
   fit2$call$data <- data_.
   
   fit2$coefficients <- unlist(data.frame(t(fit$b))) 
   fit2$varBeta <- fit$vb
   fit2$df.residual <- .Machine$double.xmax
   #fit2$dims$N <- .Machine$double.xmax
  
   if(robust) { 
      
   is_rmv <- inherits(fit,"rma.mv") 
   
   if(is_rmv){   
      mm <- try(clubSandwich:::findCluster.rma.mv(fit), silent = TRUE)
      
      is_bad_for_sandwich <- inherits(mm, "try-error")
      
   } else { is_bad_for_sandwich <- FALSE }
      
      if(is_bad_for_sandwich) message("Robust estimation was infeasible.\n")
      
      if(!is_rmv) cluster <- rownames(data_.)
      
      if(!is_bad_for_sandwich) fit2$varBeta <- vcovCR(fit, cluster=cluster, type="CR2")
   }
   
   x <- if(!full) allEffects(fit2, ...) else predictorEffects(fit2, ...)
   
   if(!is.null(index)) { 
      
     x <- if(!full) x[[index[1]]] else x[index]
   
   }

if(plot){   
xcv <- plot(x, multiline=multiline, main=main, rug=FALSE,
        confint=confint, x.var=x.var, z.var=z.var, key.args=key.args, 
        xlab=xlab, ylab=ylab, colors=colors, cex=cex, lty=lty, 
        lwd=lwd, ylim=ylim, xlim=xlim, factor.names=factor.names, band.transparency=band.transparency, 
        band.colors=band.colors, grid=grid, lattice=lattice, axes=axes, rotx=rotx, roty=roty, symbols=symbols,
        ticks.x=ticks.x, lines=lines)
   
   xcv$x.scales$tck=c(1,0)
   xcv$y.scales$tck=c(1,0)

   print(xcv)
   
   return(x)
   
} else {
   
   return(x)
   }
} 

#=================================================================================================================================================
  
cor_comb <- function(data, var1, var2) { 
   
   uni <- unique(union(data[[var1]],data[[var2]]))
   
   dat <- data.frame(t(combn(uni,2)))
   colnames(dat) <- c(var1,var2)
   dat
}  
 
  
#=================================================================================================================================================
  
add_blank_row <- function(data, n_blank_row = 1, by = "study", 
                          write = FALSE, file_name = "blank")
  {
  
 data <- full_clean(data)
 dat <- group_split(data, !!rlang::sym(by)) %>% 
  map_dfr(~ .x[1:(nrow(.x) + n_blank_row),])
 
 if(write){
   file_name <- paste0(file_name, ".csv")
   write_csv(dat, file_name, na = "")
   }
 return(dat)
}  
  
#=================================================================================================================================================  
  
# Functions for group assignment differences in studies
  
d_cluster <- function(yi_adj,n,N,icc) { yi_adj*sqrt( 1-((2*(n-1)*icc)/(N-2)) ) }

#=================================================================================================================================================
  
d_vi_cluster <- function(n, N, nT, nC, icc, w, yi){
  
  eta <-  1 + ((n - 1)*icc)
  
  #  z <- ( ((N-2)-2*(n-1)*icc)^2 + n*(N-2*n)*icc^2 + 2*(N-2*n)*icc*(1-icc) )  / 
   # (2*((N-2*n)*icc*(1-icc))^2)
    
    z <- (((N-2)*(1-icc)^2) + (n*(N-2*n)*icc^2) + (2*(N-2*n)*icc*(1-icc)))   / 
         (2*((N-2) - 2*(n-1)*icc)^2)
    
    
    (w*sqrt( (((nT+nC)/(nT*nC))*eta) + yi^2*z))^2
}

#=================================================================================================================================================
  
LRR_vi_cluster <- function(n, icc, vi){
  
  DEF <- (n - 1) * icc + 1
  DEF*vi
}  
  
#=================================================================================================================================================
  
# Correction factor for hand caculation of Hedges' g
cfactor <- function(df) metafor:::.cmicalc(df)  
  
  
#=================================================================================================================================================  
# To get SDs for each group, the best choice is the first formula:
  
range_iqr_n2sd <- function(min, max, q1, q3, n, dat) { ( (max-min) / (4*qnorm( (n-.375)/(n+.25)))) + 
  ( (q3-q1) / (4*qnorm( (n*.75-.125)/(n+.25))) ) } 

range_n2sd <- function(min, max, n) ( (max-min) / (2*qnorm( (n-.375)/(n+.25))))

iqr_n2sd <- function(q1, q3, n) ( (q3-q1) / (2*qnorm( (n*.75-.125)/(n+.25))) )

iqr_sd_cochrane <- function(q1, q3) (q3-q1)/1.35

range_f2sd <- function(min, max, f) f * (max-min) 
  
#=================================================================================================================================================  
  
t2d <- function(t, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  t/sqrt(N)
}

#=================================================================================================================================================
  
vi_d <- function(yi, n1i, n2i) { ((n1i + n2i)/(n1i * n2i)) + ((yi^2)/(2 * (n1i+n2i))) }
  
#=================================================================================================================================================
   
# To get Means for each group, the best choice is the first formula:
  
med_range_iqr_n2mean <- function(min, max, q1, q3, n, median) {
  
((2.2/(2.2+n^.75))*((min+max)/2))+((.7-(.72/n^.55))*((q1+q3)/2))+((.3+(.72/n^.55)-(2.2/(2.2+n^.75)))*median)

}  
  
#=================================================================================================================================================
  
# Missing standard deviation and no other measure of variability: The Cochrane Handbook (6.5.2.3)
# Note that this SD is the average of the SDs of the two groups and so it this same SD should be 
# inputted into the meta-analysis for both groups.
  
mdifSE_n2sd <- function(mdifSE, n1, n2)  { mdifSE / ( sqrt ( (1/n1) + (1/n2) ) ) }  
  
#=================================================================================================================================================  
  
pairwise_rma <- function(fit, type = "Tukey"){
  
mat <- setNames(rep(1,length(coef(fit))), names(coef(fit)))

summary(glht(fit, linfct=contrMat(mat, type=type)), test=adjusted("none"))
}
  
#=================================================================================================================================================                                
  
needzzsf <- c('metafor', 'clubSandwich', 'nlme', 'effects', 'lexicon', 'plotrix', 'rlang', 'fastDummies', 'multcomp', 'tidyverse')      

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))
  
options(dplyr.summarise.inform = FALSE)                        
