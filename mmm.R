
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
  
X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE]  
  
}

#===============================================================================================================================

rm.allcolNA <- function(X) { 
  
X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE]
  
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
                    
odiag <- function(x) x[col(x) != row(x)]

                    
#===============================================================================================================================
                    
shift_rows <- function(data, user_index, up = TRUE){

indx <- seq_len(nrow(data))
remain <- indx[!indx %in% user_index]

data[if(up) c(user_index, remain) else c(remain, user_index), ]
}                    
                    
#===============================================================================================================================
                    
get_error_rho <- function(fit){
  
is_V <- any(odiag(fit$V) != 0)  

if(is_V) {
  
  u <- unique(odiag(round(cov2cor(fit$V), 4)))
  u[u!=0]
  
} else 0
  
}                    

#===============================================================================================================================                    
                    
is_crossed <- function(obj){
  
  mod_struct <- clubSandwich:::parse_structure(obj)
  highest_cluster <- names(mod_struct$level_dat)[which.min(mod_struct$level_dat)]
  cluster <- mod_struct$cluster_dat[[highest_cluster]]
  out <- !clubSandwich:::test_nested(cluster, fac = mod_struct$cluster_dat)
  out[names(out) %in% obj$s.names]
}
     
#===============================================================================================================================
     
inter_index <- function(fit){
  nm <- rownames(fit$b)
  grep(":", nm)
}     

#===============================================================================================================================
     
bet_with <- function(data, contrast_col = Hypothesis) 
  {
  
out <- data %>%
    as_tibble %>%
    separate({{contrast_col}}, into = c('pre', 'post'), sep = "\\s+-\\s+", remove = FALSE) %>%
    mutate(across(pre:post, ~ map(str_extract_all(., "[A-Za-z0-9-]+\\s*\\d*"), trimws))) %>% 
    filter(lengths(map2(pre, post, intersect)) > 0) %>%
    select(-pre, -post)

as.data.frame(out)
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
  
#================================================================================================================================     
     
add_sig_funnel <- function(funnel, level=.95, col="magenta", 
                           pch=21, bg="cyan", digits=2, ...){
  
  right <- (1 - level)/2   
  
  x <- funnel$x
  y <- funnel$y
  crit <- qnorm(right, lower.tail = FALSE)
  
  ci1 <- crit*y
  x1 <- x[x > ci1]
  y1 <- y[x > ci1]
  
  suppressMessages(suppressWarnings(points(x1, y1, col=col, pch=pch, bg=bg, ...)))
  
  ci2 <- -crit*y
  x2 <- x[x < ci2]
  y2 <- y[x < ci2]
  
  total <- c(x1,x2)  
  
  suppressMessages(suppressWarnings(points(x2, y2, col=col, pch=pch, bg=bg, ...)))
  
  out <- data.frame(total = length(total), total_perc=length(total)/length(x)*100, left = length(x2), perc_left = length(x2) / length(x)*100,
                    right = length(x1), perc_right = length(x1) / length(x)*100)
  
  round(setNames(out, c("Total","Total(%)","Left","Left(%)","Right","Right(%)")), digits = digits)
} 

     
#================================================================================================================================
     
contour_funnel <- function(fit = NULL, x, y, level = c(95, 99),
                           shade = c("white","gray50"),
                           xlab = "Effect Size", sig = TRUE,
                           yaxis = "sei", ...){
yaxis <- "sei"
x <- if(!is.null(fit)) fit$yi else x
y <- if(!is.null(fit)) fit$vi else y

f1 <- metafor::funnel.default(x = x,
                             vi = y,
                          level = level, 
                          shade = shade,
                           xlab = xlab, 
                          yaxis = yaxis, ...)
  
if(sig) suppressMessages(suppressWarnings(add_sig_funnel(f1, ...)))  
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
                              
do_context <- function(data, context_vars, group_id = NULL){
  
  data <- full_clean(data)
  all_names <- names(data)
  
  id <- grep("id|group|grp|study", all_names, value = TRUE, ignore.case = TRUE)
  
  if(is.null(group_id) || !all(group_id %in% all_names)) { 
    
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
      mutate(across(all_of(num_names), list(whn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
      as.data.frame()
  }
  
  if(is_dum){
    data <- data %>%
      dummy_cols(select_columns = dum_names) %>% 
      group_by(across(all_of(group_id))) %>%                          
      mutate(across(starts_with(paste0(dum_names, "_")), list(whn = ~ . - mean(.), btw = ~ mean(.)))) %>% 
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
  
  if(!inherits(fit,c("rma.mv"))) stop("Model is not 'rma.mv'.", call. = FALSE)
  datziola <- clubSandwich:::getData(fit) %>%
    mutate(obsss = fit$slab)
  
  
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
  
  # Make visual size of effects proportional to their hat/cook's distance (estimate influence)
  cex <- cex_add_point+cex_multi_point*sqrt(if(view == 1)cook else hat)
  
  outlier_limits <- qnorm(c(.025,.5,.975))
  ylim <- range(c(outlier_limits, st_del_res_z))
  
  # Plot Leverage against Studentized residuals proportioned on cook's distances
  plot(if(view == 1) hat else cook, st_del_res_z, cex=cex, las=1, mgp=c(1.5,.3,0),
       xlab = if(view == 1) "Leverage (Hat Value)" else "Effect Influence (Cook's Dis.)", 
       ylab = "Outlier (Standardized Del. Value)",pch=19,cex.axis = .9,tcl = -.3,
        col = adjustcolor(1, .5),
       ylim = ylim)
  
  title(if(view == 1) "Size of points denote \nestimate-influencing effects\n (Cook's distances)" 
        else "Size of points denote \nleverage effects\n (Hat value)", 
        cex.main = cex_main, line = .3)
  
  abline(h=outlier_limits, lty=c(3,1,3), lwd=c(1,2,1))
  
  max_hat <- max(mean(range(hat)), boxplot.stats(hat, coef = whisker_coef)$stats[5])
  
  max_cook <- max(mean(range(cook)), boxplot.stats(cook, coef = whisker_coef)$stats[5])
  
  abline(v = if(view == 1) max_hat else max_cook, col=2)
  
  # To be outlier, an estimate must simultaneously (a) be outlying (per studentized value)
  # (b) have high leverage (per hat value), and (c) high model influence (per cook's distance)
  
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
  
  LL <- names(L)
  
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
  else as.formula(paste(as.character(y), "~ 1"))
}

#=================================================================================================================================================  
  
shorten_ <- function(vec, n = 3) { 
  
  gsub("\\s+", "", 
       sub(sprintf("^(([^,]+, ){%s}).*, ([^,]+)$", n), "\\1...,\\3", toString(vec)))  
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
                            cluster, plot=TRUE, int_only = FALSE, ...) 
{
  
  if(!inherits(fit,c("rma.mv","rma","rma.uni"))) stop("Model is not 'rma()' or 'rma.mv()'.", call. = FALSE)
  
  lm_fit <- lm(fixed_form_rma(fit), data = eval(fit$call$data), na.action = "na.omit")
  
  is_singular <- anyNA(coef(lm_fit))
  
  fit2 <- rma2gls(fit)
  
  if(is_singular) { 
    
    fit2$coefficients <- replace(lm_fit$coefficients, !is.na(lm_fit$coefficients), fit2$coefficients)
  }
  
  x <- if(!full) allEffects(fit2, ...) else predictorEffects(fit2, ...)
  
  ok <- sapply(x, function(i) !all(is.na(i$fit)))

  x <- x[ok]

  x <- if(!is.null(index)) x[index] else if(int_only) x[grep(":",names(x))] else x
  
  if(plot){   
    xcv <- plot(x, multiline=multiline, main=main, rug=FALSE,
                confint=confint, x.var=x.var, z.var=z.var, key.args=key.args, 
                xlab=xlab, ylab=ylab, colors=colors, cex=cex, lty=lty, 
                lwd=lwd, ylim=ylim, xlim=xlim, factor.names=factor.names, band.transparency=band.transparency, 
                band.colors=band.colors, grid=grid, lattice=lattice, axes=axes, rotx=rotx, roty=roty, symbols=symbols,
                ticks.x=ticks.x, lines=lines)
    
    xcv$x.scales$tck=c(1,0)
    xcv$y.scales$tck=c(1,0)
    
    invisible(return(x))
    
    xcv 
    
  } else {
    
    invisible(return(x))
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
  
random_left <- function(random_fml) {
    
  as.formula(as.character(parse(text = sub("\\|.*", "", random_fml))))
  
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
  
pairwise_rma <- function(fit, type = "Tukey", adjust = "none"){
  
  fit <- clean_reg_names(fit)
  
  mat <- setNames(rep(1,length(coef(fit))), names(coef(fit)))
  
  summary(glht(fit, linfct=contrMat(mat, type=type)), test=adjusted(adjust))
}
  

#================================================================================================================================================  
  
rma2gls <- function(fit){
  
  data_. <- eval(fit$call$data)
  form_. <- fixed_form_rma(fit)
  
  fit2 <- nlme::gls(form_., data = data_., method = fit$method, na.action = "na.omit",
            control = glsControl(singular.ok = TRUE))
  
  fit2$call$model <- form_.
  fit2$call$data <- data_.
  
  fit2$coefficients <- unlist(data.frame(t(fit$b))) 
  fit2$varBeta <- fit$vb
  
  return(fit2)
} 
  
#================================================================================================================================================
  
plot_model <- function(fit, coef = 1:5, xlab = "", ylab = "", cont_name = NULL,
                       labels = NULL, main = "", ylim = NULL, ...){
  
  fit <- clean_reg_names(fit)
  
  ff <- conf_int(fit, "CR2")
  ci <- c(ff$CI_L[coef], ff$CI_U[coef])
  
  cont <- function(fit, coef, cont_name, ...){
    
    data <- clubSandwich:::getData(fit)
    
    curve(ff$beta[1] + x*ff$beta[coef], from=min(data[[cont_name]],na.rm = TRUE), 
          to=max(data[[cont_name]],na.rm = TRUE), xlab = xlab, ylab = ylab, 
          mgp=c(1.5,.4,0), ...)
    title(main = main, cex.main=.9, line = .2)
  }
  
  if(!is.null(cont_name)) return(
    
    
    cont(fit, coef, cont_name, ...)
    
  )
  
  plot(coef-1,ff$beta[coef], type="b", ylim = 
         if(is.null(ylim)) range(ci) else ylim, xaxt="n",
       ylab = ylab, xlab = xlab, mgp = c(1.5,.4,0), ...)
  
  abline(h=0,col=8,lty=3)
  segments(coef-1, ci[coef], y1=ci[coef+length(coef)])
  
  axis(1, at=coef-1, labels = if(is.null(labels)) rownames(ff)[coef] else labels,
       mgp=c(1.5,.4,0), ...)
  
  points(coef-1,ff$beta[coef], pch=19)
  title(main,cex.main=.9,line = .2)
}  

#=================================================================================================================================================

lrr2percent <- function(lrr) (exp(lrr) - 1)*1e2       
      
#=================================================================================================================================================       
       
lrr <- function(mT, nT, sdT, mC, nC, sdC){
  
f <- function(mT, nT, sdT, mC, nC, sdC){
  
  means <- c(mT, mC)
  ns <- c(nT, nC)
  variances <- c(sdT, sdC)^2
  
  bc <- log(means) - variances/(2 * ns * means^2)
  LRR <- as.numeric(bc[1] - bc[2])
  
  LRR_vi <- sum(variances/(ns * means^2))
  
  percent_dif <- (exp(LRR)-1)*1e2
  
  c(LRR = LRR, percent_dif = percent_dif, 
     LRR_vi = LRR_vi)
}

pmap_dfr(list(mT=mT, nT=nT, sdT=sdT, mC=mC, nC=nC, sdC=sdC), f)

}       
     
#=================================================================================================================================================       
       
clean_reg_names <- function(fit) {
  
  fmla <- fixed_form_rma(fit)
  vec <- rownames(fit$b)
  
  vec <- clean_reg(fmla, vec)
  if(fit$int.only) vec[vec=="Intercept"||vec==""] <- "Overall Effect"
  rownames(fit$b) <- vec
  rownames(fit$beta) <- vec
  rownames(fit$vb) <- colnames(fit$vb) <- vec
  return(fit)
}     
  
#=================================================================================================================================================       
       
results_rma2 <- function(fit, digits = 3, robust = FALSE, blank_sign = "", 
                        cat_shown = 1, shift_up = NULL, shift_down = NULL, 
                        drop_rows = NULL, drop_cols = NULL, QE = FALSE, sig = FALSE,
                        clean_names = TRUE){
  
  if(!inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  fixed_eff <- is.null(fit$random)
  
  if(clean_names) fit <- clean_reg_names(fit)
  
  cr <- if(!fixed_eff) is_crossed(fit) else FALSE
  
  if(robust & any(cr) || robust & fixed_eff) { 
    
    robust <- FALSE
    message("Robust estimation not available for models with", if(any(cr))" crossed random-" else " only fixed-", "effects.")
  }
  
  res <- if(!robust) { 
    
    a <- coef(summary(fit))
    
    nm <- c("Estimate","SE","t-value","Df","p-value","Lower","Upper")
    
    nm <- if(fit$test == "t") { nm } else { nm[3] <- "z-value";
    nm[-4] }
    
    setNames(a, nm)
    
  } else {
    
    a <- as.data.frame(conf_int(fit, vcov = "CR2"))
    
    a$p_Satt <- coef_test(fit, vcov = "CR2")$p_Satt
    
    a <- a[c(1:3,6,4:5)]
    
    setNames(a, c("Estimate","SE","Df","p-value","Lower","Upper"))
  }
  
  u <- get_error_rho(fit)
  cte <- length(u) == 1
  
  d6 <- data.frame(r = if(cte) u else mean(u, na.rm = TRUE), 
                   row.names = paste0("Within Corr.(",if(cte) "constant" else "average",")"))
  
  if(QE){
    qe <- data.frame(Estimate = fit$QE, Df = nobs.rma(fit), 
                     pval = fit$QEp, row.names = "QE") %>%
      dplyr::rename("p-value"="pval") 
    
    res <- bind_rows(res, qe)
  }
  
  
  if(fixed_eff) { 
    
    out <- roundi(dplyr::bind_rows(res, d6), digits = digits)
    out[is.na(out)] <- blank_sign
    
    if(sig){ 
      
      p.values <- as.numeric(out$"p-value")
      
      Signif <- symnum(p.values, corr = FALSE, 
                       na = FALSE, cutpoints = 
                         c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " ")) 
      
      out <- add_column(out, Sig. = Signif, .after = "p-value")
    }
    
    if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
    if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
    if(!is.null(drop_rows)) out <- out[-drop_rows, ]
    
    out <- dplyr::select(out, -dplyr::all_of(drop_cols))
    
    return(out)
  }
  
  res <- rbind(res, "(RANDOM)" = NA)
  
  if(fit$withS){
    
    d1 <- data.frame(Sigma = sqrt(fit$sigma2), 
                     row.names = paste0(names(cr), ifelse(cr,"(Crossed)","(Intercept)"))) 
    
  } else { d1 <- NULL}
  
  if(fit$withG){
    
    h <- paste(fit$struct[1], "Corr.")
    is_un <- fit$struct[1] == "UN"
    is_gen <- fit$struct[1] == "GEN"
    is_diag <- fit$struct[1] == "DIAG"
    is_simple <- length(fit$tau2) == 1
    
    rnm <- paste("Level:", tail(fit$g.names,1))
    clnm <- clean_GH_names(fit)
    
    d2 <- data.frame(Tau = sqrt(fit$tau2), 
                     row.names = paste0(if(!is_simple) clnm else fit$g.names[1],
                                        paste0(if(is_diag)"(Uncor." 
                                               else "(Cor.",if(!is_simple & !is_gen) paste0(" ", fit$g.names[1]),")")))
    
    d2 <- rbind(NA, d2)
    rownames(d2)[1] <- rnm
    
    d3 <- data.frame(Rho = fit$rho, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~") 
                     else paste0(h,"(",shorten_(clnm, cat_shown),")")) 
    
  } else { d2 <- NULL; d3 <- NULL}
  
  if(fit$withH){
    
    h <- paste(fit$struct[2], "Corr.")
    is_un <- fit$struct[2] == "UN"
    is_gen <- fit$struct[2] == "GEN"
    is_diag <- fit$struct[2] == "DIAG"
    is_simple <- length(fit$gamma2) == 1
    
    rnm <- paste("Level:", paste0(tail(fit$h.names,1)," "))
    
    clnm <- clean_GH_names(fit, G=FALSE)
    
    d4 <- data.frame(Gamma = sqrt(fit$gamma2), 
                     row.names = paste0(if(!is_simple) clnm else fit$h.names[1],
                                        paste0(if(is_diag)"(Uncor." 
                                               else "(Cor.",if(!is_simple) paste0(" ",if(!is_gen)fit$h.names[1]),") "))) 
    
    d4 <- rbind(NA, d4)
    rownames(d4)[1] <- rnm
    
    d5 <- data.frame(Phi = fit$phi, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~ ")
                     else paste0(h,"(",shorten_(clnm, cat_shown),") "))
    
  } else { d4 <- NULL; d5 <- NULL}
  
  out <- roundi(dplyr::bind_rows(res, d1, d2, d3, d4, d5, d6), digits = digits)
  
  out[is.na(out)] <- blank_sign
  
  if(sig){ 
    
    p.values <- as.numeric(out$"p-value")
    
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " ")) 
    
    out <- add_column(out, Sig. = Signif, .after = "p-value")
  }
  
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  
  out <- dplyr::select(out, -dplyr::all_of(drop_cols))
  
  return(out)
}
                     
#=================================================================================================================================================
                     
results_rma <- function(fit, digits = 3, robust = FALSE, blank_sign = "", 
                        cat_shown = 1, shift_up = NULL, shift_down = NULL, 
                        drop_rows = NULL, drop_cols = NULL, QE = FALSE, sig = FALSE,
                        clean_names = TRUE){
  
  if(!inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  fixed_eff <- is.null(fit$random)
  
  if(clean_names) fit <- clean_reg_names(fit)
  
  cr <- if(!fixed_eff) is_crossed(fit) else FALSE
  
  if(robust & any(cr) || robust & fixed_eff) { 
    
    robust <- FALSE
    message("Robust estimation not available for models with", if(any(cr))" crossed random-" else " only fixed-", "effects.")
  }
  
  res <- if(!robust) { 
    
    a <- coef(summary(fit))
    
    nm <- c("Estimate","SE","t","Df","p-value","Lower","Upper")
    
    nm <- if(fit$test == "t") { nm } else { nm[3] <- "z";
    nm[-4] }
    
    setNames(a, nm)
    
  } else {
    
    a <- as.data.frame(conf_int(fit, vcov = "CR2"))
    
    a$p_Satt <- coef_test(fit, vcov = "CR2")$p_Satt
    
    a <- a[c(1:3,6,4:5)]
    
    setNames(a, c("Estimate","SE","Df","p-value","Lower","Upper"))
  }
  
  u <- get_error_rho(fit)
  cte <- length(u) == 1
  
  d6 <- data.frame(r = if(cte) u else mean(u, na.rm = TRUE), 
                   row.names = paste0("Within Corr. (",if(cte) "constant" else "average",")"))
  
  if(QE){
    qe <- data.frame(Estimate = fit$QE, Df = nobs.rma(fit), 
                     pval = fit$QEp, row.names = "QE") %>%
      dplyr::rename("p-value"="pval") 
    
    res <- bind_rows(res, qe)
  }
  
  
  if(fixed_eff) { 
    
    out <- roundi(dplyr::bind_rows(res, d6), digits = digits)
    out[is.na(out)] <- blank_sign
    
    if(sig){ 
      
      p.values <- as.numeric(out$"p-value")
      
      Signif <- symnum(p.values, corr = FALSE, 
                       na = FALSE, cutpoints = 
                         c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                       symbols = c("***", "**", "*", ".", " ")) 
      
      out <- add_column(out, Sig. = Signif, .after = "p-value")
    }
    
    if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
    if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
    if(!is.null(drop_rows)) out <- out[-drop_rows, ]
    
    out <- dplyr::select(out, -dplyr::all_of(drop_cols))
    
    return(out)
  }
  
  res <- rbind(res, "(RANDOM)" = NA)
  
  Sys.setlocale(category = "LC_ALL", locale = "Greek")
  
  if(fit$withS){
    
    d1 <- data.frame(Sigma = sqrt(fit$sigma2), 
                     row.names = paste0(names(cr), ifelse(cr," (Crossed Ave.)"," (Nested Ave.)"))) 
    
    d1 <- setNames(d1, intToUtf8(963))
  } else { d1 <- NULL}
  
  if(fit$withG){
    
    h <- paste(fit$struct[1], "Corr.")
    is_un <- fit$struct[1] == "UN"
    is_gen <- fit$struct[1] == "GEN"
    is_diag <- fit$struct[1] == "DIAG"
    is_simple <- length(fit$tau2) == 1
    
    rnm <- paste("Level:", tail(fit$g.names,1))
    clnm <- clean_GH_names(fit)
    
    d2 <- data.frame(Tau = sqrt(fit$tau2), 
                     row.names = paste0(if(!is_simple) clnm else fit$g.names[1],
                                        paste0(if(is_diag)" (Uncor. " 
                                               else " (Cor. ",if(!is_simple & !is_gen) paste0(" ", fit$g.names[1]),")")))
    
    d2 <- setNames(d2, intToUtf8(964))
    
    d2 <- rbind(NA, d2)
    rownames(d2)[1] <- rnm
    
    d3 <- data.frame(Rho = fit$rho, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~") 
                     else paste0(h," (",shorten_(clnm, cat_shown),")")) 
    
    d3 <- setNames(d3, intToUtf8(961))
    
  } else { d2 <- NULL; d3 <- NULL}
  
  if(fit$withH){
    
    h <- paste(fit$struct[2], "Corr.")
    is_un <- fit$struct[2] == "UN"
    is_gen <- fit$struct[2] == "GEN"
    is_diag <- fit$struct[2] == "DIAG"
    is_simple <- length(fit$gamma2) == 1
    
    rnm <- paste("Level:", paste0(tail(fit$h.names,1)," "))
    
    clnm <- clean_GH_names(fit, G=FALSE)
    
    d4 <- data.frame(Gamma = sqrt(fit$gamma2), 
                     row.names = paste0(if(!is_simple) clnm else fit$h.names[1],
                                        paste0(if(is_diag)" (Uncor. " 
                                               else " (Cor. ",if(!is_simple) paste0(" ",if(!is_gen)fit$h.names[1]),") "))) 
    
    d4 <- setNames(d4, intToUtf8(933))
    
    d4 <- rbind(NA, d4)
    rownames(d4)[1] <- rnm
    
    d5 <- data.frame(Phi = fit$phi, 
                     row.names = if(is_un || is_gen) apply(combn(clnm,2),2,paste0, collapse = "~ ")
                     else paste0(h," (",shorten_(clnm, cat_shown),") "))
    
    d5 <- setNames(d5, intToUtf8(966))
    
  } else { d4 <- NULL; d5 <- NULL}
  
  out <- roundi(dplyr::bind_rows(res, d1, d2, d3, d4, d5, d6), digits = digits)
  
  out[is.na(out)] <- blank_sign
  
  if(sig){ 
    
    p.values <- as.numeric(out$"p-value")
    
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " ")) 
    
    out <- add_column(out, Sig. = Signif, .after = "p-value")
  }
  
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  
  out <- dplyr::select(out, -dplyr::all_of(drop_cols))
  
  return(out)
}                     
                     
#=================================================================================================================================================                   
                   
roundi <- function(x, digits = 7){
  
  if(!inherits(x, c("data.frame","tibble"))) stop("'x' must be a 'data.frame'.", call. = FALSE)
  
  num <- sapply(x, is.numeric)
  
  x[num] <- lapply(x[num], round, digits)
  
  return(x)
}                    
                   
#=================================================================================================================================================
                   
mc_rma <- function(fit, specs, var = NULL, by = NULL, horiz = TRUE, 
                   adjust = "tukey", compare = FALSE, plot = FALSE, 
                   reverse = TRUE, digits = 3, xlab = "Estimated Effect", 
                   shift_up = NULL, shift_down = NULL, drop_rows = NULL, 
                   drop_cols = 7, full = FALSE, na.rm = TRUE, ...){
  
  
  if(!inherits(fit, "rma.mv")) stop("Model is not 'rma.mv()'.", call. = FALSE)
  if(length(inter_index(fit))==0) full <- TRUE
  
  dat_ <- eval(fit$call$data)
  lm_fit <- lm(fixed_form_rma(fit), data = dat_, na.action = "na.omit")
  lm_fit$call$data <- dat_
  
  is_singular <- anyNA(coef(lm_fit))
  
  fit <- rma2gls(fit)
  
  if(is_singular) { 
    
    fit$coefficients <- replace(lm_fit$coefficients, !is.na(lm_fit$coefficients), fit$coefficients)
    
    fit <- suppressMessages(emmeans::ref_grid(fit))
    fit@nbasis <- suppressMessages(emmeans::ref_grid(lm_fit)@nbasis)
  }
  
  infer <- c(TRUE, TRUE)
  
  ems <- if(is.null(var)){
    
    emmeans(object = fit, specs = specs, infer = infer, ...)
    
  } else {
    
    emtrends(object = fit, specs = specs, var = var, infer = infer, ...)
    
  }
  
  is_pair <- "pairwise" %in% as.character(specs)
  
  if(is_pair){
    
    if(plot) print(plot(ems, by = by, comparisons = compare, horizontal = horiz, adjust = adjust, xlab = xlab)) 
    
    out <- as.data.frame(pairs(ems, reverse = reverse, each = "simple", infer = infer)$emmeans)[c(1:3,7,4,8,5:6)]
    
    names(out) <- c("Hypothesis","Estimate","SE","t","Df","p-value","Lower","Upper")
    
    p.values <- as.numeric(out$"p-value")
    Signif <- symnum(p.values, corr = FALSE, 
                     na = FALSE, cutpoints = 
                       c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " "))
    
    out <- add_column(out, Sig. = Signif, .after = "p-value")
    
    out <- roundi(out, digits = digits)
    
    if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
    if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
    if(!is.null(drop_rows)) out <- out[-drop_rows, ]
    
    if(!full) out <- bet_with(out, Hypothesis)
    
    out <- dplyr::select(out, -dplyr::all_of(drop_cols))
    
    if(na.rm) out <- na.omit(out)
    
    return(out)
  }
  
  else {
    
    ems
  }
}
  
                     
#=================================================================================================================================================
                     
mc_robust_rma <- function(fit, constraints, vcov = "CR2", test = "HTZ", digits = 3,
                          shift_up = NULL, shift_down = NULL, drop_rows = NULL, 
                          clean_names = TRUE, ...){
  
  if(clean_names) obj <- clean_reg_names(fit)
  
  out <- roundi(as.data.frame(Wald_test(obj=obj, constraints=constraints, vcov=vcov, 
                                        test=test, tidy=TRUE, ...))[-c(2,4)], digits)
  
  out <- setNames(out, c("Hypothesis","F","Df1","Df2","p-value"))
  
  p.values <- as.numeric(out$"p-value")
  Signif <- symnum(p.values, corr = FALSE, 
                   na = FALSE, cutpoints = 
                     c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " "))
  
  out <- add_column(out, Sig. = Signif, .after = "p-value")
  if(!is.null(shift_up)) out <- shift_rows(out, shift_up)
  if(!is.null(shift_down)) out <- shift_rows(out, shift_down, up = FALSE)
  if(!is.null(drop_rows)) out <- out[-drop_rows, ]
  return(out)
}
 
#=================================================================================================================================================
                                      
random_GH_form <- function(fit, G = TRUE){
  
  fm <- fit$random
  if(G) fm[[1]] else fm[[2]]
  
}

#=================================================================================================================================================
                                      
clean_GH_names <- function(fit, G = TRUE) {
  
  fmla <- random_left(random_GH_form(fit, G = G))
  vec <- if(G) rownames(fit$G) else rownames(fit$H)
  clean_reg(fmla, vec)
  
}                                     
                       
#=================================================================================================================================================
                       
clean_reg2 <- function(fmla, vec) {
  v1 <- as.character(attr(terms(fmla), "variables"))[-1L] 
  v2 <- setdiff(vec, v1)
  #v1 <- gsub(r"(([\^$.?*|+()[{]))", r"(\\\1)", v1)
  v1 <- gsub("([\\\\^$.?*|+()[\\]{}])", "\\\\\\1", v1, perl = TRUE) # For compatibilty with older R versions
  v3 <- gsub(paste(v1, collapse = "|"), "", v2)
  vec[vec %in% v2] <- v3
  vec[vec=="intrcpt"] <- "Intercept"
  return(vec) 
}                       

#=================================================================================================================================================                     
  
clean_reg <- function(fm, nm) {
  vars <- vapply(attr(terms(fm), "variables"), deparse, "")[-1L]
  subpat <- paste0(gsub("([()])", "\\\\\\1", vars), collapse = "|")
  l <- rapply(strsplit(nm, ":"), sub, how = "list",
              perl = TRUE,
              pattern = sprintf("^(?!(%1$s)$)(%1$s)(.+)$", subpat),
              replacement = "\\3")
  vec <- vapply(l, paste0, "", collapse = ":")
  vec[vec=="intrcpt"] <- "Intercept"
  return(vec)
}                     
                                        
#=================================================================================================================================================                       
set_rownames <- function (object = nm, nm) 
{
  rownames(object) <- nm
  object
}       
                     
#=================================================================================================================================================                                      
#=================================================================================================================================================
#=================================================================================================================================================
                     
round_effects <- function(data, digits, yi_vi_names = c("yi","vi")){
  
  eff <- names(data) %in% yi_vi_names
  data[eff] <- lapply(data[eff], round, digits)
  return(data)
}

#==============

random_rows <- function(x) {
  
  is_char <- is.character(x) 
  is_fact <- is.factor(x)
  
  unique_levels <- if(is_fact) levels(x) else unique(x)
  
  r_min <- if(is_char|is_fact) 0 else min(unique_levels)-1 
  r <- c(r_min, unique_levels)
  n <- sample(r, 1)
  #If n = r[1] select all rows else 
  #select row for corresponding x value
  if(n == r[1]) TRUE else x == n
}

#================

sample_order <- function(data, vars, cat_var_to_remove = "time"){
  g_by <- vars[1]
  vars <- setdiff(vars, c(g_by, cat_var_to_remove))
  
  Fn <- function(data, var, g_by) { 
    data %>% 
      group_by(.data[[g_by]]) %>% 
      filter(random_rows(.data[[var]])) %>% 
      ungroup() %>% as.data.frame() 
  }
  purrr::reduce(vars, Fn, g_by, .init=data) %>%
    group_by(.data[[g_by]]) %>% 
    mutate(across(all_of(vars), 
                  function(x) if(is.numeric(x)) 
                    as.numeric(as.factor(x)) else x)) %>%
                              as.data.frame()
}



#==============

rbetab <- function(n, mu.p, disp){
  
  if(mu.p < 0 || mu.p > 1) { message("Warning: 'mu.p' is 'average probability' of a 'beta dist.' bound between '0' & '1'.") ;
    mu.p[mu.p < 0] <- 0 ;
    mu.p[mu.p > 1] <- 1 }
  
  shape1 <- mu.p * disp
  shape2 <- (1 - mu.p) * disp
  rbeta(n, shape1 = shape1, shape2 = shape2)
}

#==============


sim_rma.mv.dat <- function(..., seed = NULL, mean_sigma = 0.75, sd_df = 1,
                           cat_var_to_remove = NULL,
                           n_high_to_row_remove = 1,
                           drop_levels = NULL, LRR = FALSE, mu.p_sigma = .5,
                           mu.p_disp = 10, digits = 6, ordered_cat_var_to_remove = TRUE
){
  
  set.seed(seed)
  dat <- expand_grid(...)
  
  dots <- rlang::list2(...) 
  
  cat_var_to_remove <- if(!is.null(cat_var_to_remove)){
    cat_var_to_remove <- rlang::ensym(cat_var_to_remove)
  } else cat_var_to_remove
  
  high_name <- names(dots)[1]
  
  
  rows <- sapply(group_split(dat, !!!rlang::syms(high_name)), nrow)
  
 # set.seed(seed)
  
  dat <- if(!LRR){
    mutate(dat, yi = as.vector(unlist(mapply(rnorm, n = rows, 
                                                   mean = rnorm(rows, mean_sigma, .5), 
                                                   sd = rchisq(rows, sd_df)))),
           
           vi = as.vector(unlist(mapply(runif, n = rows))) )
  } else {
    
    mutate(dat, yi = as.vector(unlist(mapply(rbetab, n = rows, 
                                             mu.p = rbetab(rows, mu.p_sigma, mu.p_disp), 
                                             disp = rnbinom(rows, mu.p_disp, mu = mu.p_disp)))),
           
           vi = yi * as.vector(unlist(mapply(runif, n = rows, min = .001, max = .9))) )
  }
  

  dat <- round_effects(dat, digits)

  
  if(!is.null(cat_var_to_remove)){ 
    
    set.seed(seed)
    
    outcome_to_remove <- if(ordered_cat_var_to_remove) sort(unique(dat[[cat_var_to_remove]]))[-(1:2)] else unique(dat[[cat_var_to_remove]])
    study_available <- unique(dat[[high_name]])
    
    
    for(i in outcome_to_remove) {
      
      studies_to_remove <- sample(study_available, size = n_high_to_row_remove)
      study_available <- setdiff(study_available, studies_to_remove)
      #message('Dropping ', cat_var_to_remove, " ", i, " in ", high_name, " ", studies_to_remove)
      dat <- dat %>%
        filter(
          !(    .data[[high_name]] %in% studies_to_remove &
                  if(ordered_cat_var_to_remove) !!cat_var_to_remove  >= i else !!cat_var_to_remove %in% i
          ))
    }
    
  }
  
  #set.seed(seed)
  
  dat %>%
    sample_order(names(dots), cat_var_to_remove) %>%
  select(-tidyselect::all_of(drop_levels)) %>%
    mutate(row_id = row_number())
    
}

#=======

any_nesting <- function(data, ...){
  
input <- is_nestor(data, ...)

nesting. <- unlist(Filter(
  length,
  lapply(
    names(input), 
    function(i) {
      k <- which(input[[i]])
      if (length(k)) paste0(i, '/', names(input[[i]])[k])
    }
  )
))

 if(is.null(nesting.)) nesting. else noquote(nesting.)
}

#=======


hlister <- function(data, highest_level, str_cols){
  
  hlist <- data %>%
  dplyr::group_by(!!!highest_level) %>%
  dplyr::mutate(grp = dplyr::across(tidyselect::all_of(str_cols), ~ {
    tmp <- dplyr::n_distinct(.)
    dplyr::case_when(tmp == 1 ~ 1, tmp == n() ~ 2, tmp > 1 & tmp < n() ~ 3,  TRUE ~ 4)
  }) %>%
    purrr::reduce(stringr::str_c, collapse = "")) %>%
  dplyr::ungroup(.) %>%
  dplyr::group_split(grp, .keep = FALSE)

Filter(NROW, rev(hlist))
}

#=======

disply_highest_level_names <- function(res, sss, struc){

typic <- function(vec) vec[ceiling(length(vec)/2)]

lapply(res, function(i){
  nr <- sapply(split(i, i[[sss]]), nrow);
  study_type <- if(struc == "typical") {typic(as.numeric(names(table(nr))))
  } else if(struc == "simple") {min(as.numeric(names(table(nr))))
  } else {max(as.numeric(names(table(nr))))};
  names(nr)[nr == study_type][1]
  })
}

#=======


meta_tree4 <- function(data = NULL, ..., highest_level_name = NULL, reset = TRUE,
                      structure = c("simple","typical","complex"), 
                      output_highest_level = FALSE,
                      toplab = NULL, cex = 1, main = NULL, rowcount = FALSE, 
                      main_extra_name = FALSE, subset, 
                      n_high_to_row_remove = 15, cat_var_to_remove = "time", 
                      ordered_cat_var_to_remove = TRUE, random_rm = TRUE,
                      drop_levels = NULL, row_id = FALSE, seed = NULL) 
{
  
  
  if(reset){
    graphics.off()
    org.par <- par(no.readonly = TRUE)
    on.exit(par(org.par))
  }
  
  
 if(!is.null(data)){
   
 data <- rm.colrowNA(trim_(data)) %>%
    mutate(row_id = row_number())

  dot_cols <- rlang::ensyms(...)
  str_cols1 <- purrr::map_chr(dot_cols, rlang::as_string)
  
  idx <- str_cols1 %in% names(data)
  if(!all(idx)) stop(toString(dQuote(str_cols1[!idx]))," not found in the 'data'.", call. = FALSE)

  str_cols <- str_cols1[-1]

  highest_level <- dot_cols[[1]]

nesting_suggest <- any_nesting(data, !!!dot_cols)

  ss <- substitute(highest_level)
  sss <- deparse(ss)
  
  
  if(!missing(subset)){
    
    s <- substitute(subset)
    data <- filter(data, eval(s))
  }
  
  data <- data %>%
    dplyr::select(!!! dot_cols)# %>%
    #dplyr::select(-tidyselect::all_of(drop_levels)) 
  
 }  else {
   
dots <- rlang::list2(...)    
   
data <- tidyr::expand_grid(...) 
if(row_id) data <- mutate(data, row_id = row_number())

dotnames <- names(dots)#[!names(dots) %in% drop_levels]
str_cols <- dotnames[-1]

highest_level <- rlang::sym(dotnames[1])

ss <- substitute(highest_level)
sss <- dotnames[1]

if(any(drop_levels %in% sss)) stop("Highest level can't be dropped.", call. = FALSE)

data <- if(random_rm) {
  
  cat_var_to_remove <- if(!is.null(cat_var_to_remove)){
    cat_var_to_remove <- rlang::ensym(cat_var_to_remove)
  } else {cat_var_to_remove}
  
  
  if(!is.null(cat_var_to_remove)){    
    
    outcome_to_remove <- if(ordered_cat_var_to_remove) sort(unique(data[[cat_var_to_remove]]))[-(1:2)] else unique(data[[cat_var_to_remove]])
    study_available <- unique(data[[sss]])
    
    set.seed(seed)
    
    for(i in outcome_to_remove) {
      
      studies_to_remove <- sample(study_available, size = n_high_to_row_remove)
      study_available <- setdiff(study_available, studies_to_remove)
      #message('Dropping ', cat_var_to_remove, " ", i, " in ", high_name, " ", studies_to_remove)
      data <- data %>%
        filter(
          !(    .data[[sss]] %in% studies_to_remove &
                  if(ordered_cat_var_to_remove) !!cat_var_to_remove  >= i else !!cat_var_to_remove %in% i
          ))
    } 
  }
  
  data %>%
    sample_order(dotnames, cat_var_to_remove) %>%
   dplyr::select(-tidyselect::all_of(drop_levels)) 
  
} else { data }

nesting_suggest <- any_nesting(data, !!!rlang::syms(names(data)))

}
  str_cols <- str_cols[!(str_cols %in% drop_levels)]
 
if(is.null(highest_level_name)){
    
    struc <- match.arg(structure) 
    
    res <- hlister(data, highest_level, str_cols)    
    
    main_no. <- sapply(res, function(i) length(unique(i[[sss]])))
    
    nms <- disply_highest_level_names(res, sss, struc)  
    
    list2plot <- lapply(seq_along(res),function(i) filter(res[[i]], eval(ss) == nms[i]))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL),
                      mar = c(2.5, 2.6, 1.8, .5))
 
      }
    
    main <- if(is.null(main)) ifelse(main_no. > 1, pluralify_(sss), sss) else main
    
    main <- paste(main_no., main)
    
    if(main_extra_name) main <- paste0(main, " [",nms,"]")
    
    invisible(lapply(seq_along(list2plot), function(i) data.tree_(list2plot[[i]], main = main[i], toplab, cex, rowcount)))
    
    invisible(if(output_highest_level & length(nesting_suggest) ==0) list(Structure_Types = res, data_used = data) else 
      if(output_highest_level & output_highest_level !=0) list(tibble(`Nesting Suggestion:` = nesting_suggest), Structure_Types = res, data_used = data) else
        if(length(nesting_suggest)!=0) list(tibble(`Nesting Suggestion:` = nesting_suggest), data_used = data)) 
    
    
  } else {
    
    highest_level_name <- trimws(highest_level_name)
    highest_level_names <- unique(data[[sss]])
    
    idx <- highest_level_name %in% highest_level_names 
    
    if(!all(idx)) stop(toString(dQuote(highest_level_name[!idx]))," not found in the ", paste0("'",sss,"'", " column."), call. = FALSE)
    
    list2plot <- lapply(highest_level_name, function(i) filter(data, eval(ss) == i))
    
    LL <- length(list2plot)
    
    if(LL > 1L) { par(mfrow = n2mfrow(LL)) }
    
    invisible(lapply(list2plot, data.tree_, toplab, cex, rowcount))
    
    if(length(nesting_suggest)!=0) tibble(`Nesting Suggestion:` = nesting_suggest)
    
  }
}                   
                     
                     
#=================================================================================================================================================
                     
needzzsf <- c('metafor', 'clubSandwich', 'nlme', 'effects', 'lexicon', 'plotrix', 'rlang', 'fastDummies', 'multcomp','emmeans','tidyverse')      

not.have23 <- needzzsf[!(needzzsf %in% installed.packages()[,"Package"])]
if(length(not.have23)) install.packages(not.have23)

suppressWarnings(
  suppressMessages({ 
    
    invisible(lapply(needzzsf, base::require, character.only = TRUE))
    
  }))
  
options(dplyr.summarise.inform = FALSE)                        
