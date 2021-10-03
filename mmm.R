
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

cat_overlap1 <- function(data, study_id, ...){
  
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
    cat_names <- paste0(rlang::as_string(.x), c(".x", ".y"))
    
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
      tidyr::pivot_wider(names_from = cat_names[2], values_from = n) %>%
      dplyr::rename(`Moderator Category` = cat_names[1]) %>% 
      dplyr::mutate(`Moderator Category` = as.character(`Moderator Category`)) %>% 
      #dplyr::slice(match(names(.)[-1], `Moderator Category`))
      dplyr::arrange(factor(`Moderator Category`, levels = names(.)[-1]))
    
    out2 <- out1[-1]
    out2[upper.tri(out2)] <- "-"
    dplyr::bind_cols(out1[1],out2)
    
  }), cat_nms)
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
    
    #====
    hlist <- data %>%
    dplyr::group_by({{highest_level}}) %>%
      dplyr::mutate(grp = dplyr::across(tidyselect::all_of(str_cols), ~ {
        tmp <- dplyr::n_distinct(.)
        #dplyr::case_when(tmp  == 1 ~ 1, tmp == n() ~ 2, TRUE ~ 3)
        dplyr::case_when(tmp  == 1 ~ 1, tmp == n() ~ 2, tmp >1 & tmp < n() ~ 3,  TRUE ~ 4)
      }) %>%
        purrr::reduce(stringr::str_c, collapse = "")) %>%
      dplyr::ungroup(.) %>%
      dplyr::group_split(grp, .keep = FALSE)
    #====
    
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

  has_G <- isTRUE(fit$withG)
  has_H <- isTRUE(fit$withH)
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
                        
study_struct <- function(..., raw = FALSE, row_id = FALSE, missing = FALSE, seed = NULL, comparison = TRUE){
  
  dots <- rlang::list2(...) 
  inp_ls <- purrr::map(dots, ~ purrr::map(.x, seq_len)) %>% transpose %>% 
    purrr::map(purrr::compact)
  
  data_ls <- purrr::map(inp_ls, 
                        ~ tidyr::expand_grid(!!! .x,
                                             comparison = c("control","treatment")))
  
  out <- if(!raw) { purrr::map2(data_ls, inp_ls, ~ .x %>% 
                                  dplyr::group_by(!!! rlang::syms(names(.y))) %>%
                                  dplyr::summarise(comparison = stringr::str_c(sort(comparison, decreasing = TRUE), 
                                                                          collapse = ' vs '), .groups = 'drop'))
  } else { data_ls }
  
  
  res <- dplyr::bind_rows(out, .id = "study")
  
  res <- if(comparison) res else dplyr::select(res, -comparison)
  
  res <- if(row_id) dplyr::bind_cols(res, row_id = seq_len(nrow(res))) else res
  
  if(missing){
    
    set.seed(seed)
    
    res %>% dplyr::group_by(study) %>% 
      dplyr::sample_n(sample(dplyr::n(), 1)) %>% 
      dplyr::ungroup()
    
  } else {
    
    res
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
                    font.lab = 2, xlab = "Pre-test", xaxt = "n",
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
  
  all_names <- names(data)
  
  id <- grep("id|group|grp", all_names, value = TRUE, ignore.case = TRUE)
  
  if(!all(group_id %in% all_names)) { 
    
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
                        
needzzsf <- c('metafor', 'clubSandwich', 'lexicon', 'plotrix', 'rlang', 'fastDummies', 'tidyverse')      

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))
  
options(dplyr.summarise.inform = FALSE)                        
