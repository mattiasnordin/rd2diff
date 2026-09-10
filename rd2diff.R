create_phi <- function(x, c){
  G <- length(x)
  inds <- setdiff(2:G, c)
  phi <- (x[c]-x[c-1]) / (x[inds] - x[inds-1])
  return(phi)
}

create_shifted_matrix <- function(nrows){
  R1 <- matrix(0, nrows, nrows+1)
  if (nrows > 0){
    R1[1:nrows, 1:nrows] <- R1[1:nrows, 1:nrows] + diag(nrows)
    R1[1:nrows, 2:(nrows+1)] <- R1[1:nrows, 2:(nrows+1)] - diag(nrows)
  }
  return(R1)
}

create_R_matrix <- function(x, c){
  G <- length(x)
  phi <- create_phi(x, c)
  R1 <- create_shifted_matrix(c-2)
  R1 <- rbind(R1, matrix(0, G-c, c-1))
  R2 <- create_shifted_matrix(G-c)
  R2 <- rbind(matrix(0, c-2, G-c+1), R2)
  R <- cbind(R1, R2)
  R <- R * phi
  R[, c-1] <- R[, c-1] - 1
  R[, c] <- R[, c] + 1
  return(R)
}

overall_t <- function(tstats){
  Tmin <- min(tstats)
  Tmax <- max(tstats)
  if (sign(Tmin) != sign(Tmax)){
    return(0)
  }
  else {
    return(min(c(abs(Tmin), abs(Tmax))))
  }
}

get_ci <- function(b_est, stde, alphaval){
  z_alpha2 <- qnorm(1-alphaval/2)
  lb <- b_est - stde * z_alpha2
  ub <- b_est + stde * z_alpha2
  return(c(min(lb), max(ub)))
}

parse_vectors <- function(depvar, runvar, controls=NULL){
  return(
    list(
      depvar = depvar,
      runvar = runvar,
      controls = controls
    )
  )
}

parse_formula <- function(formula, data) {
  if (!inherits(formula, "formula")){
    stop("Formula incorrectly specified.", call. = FALSE)
  }

  depvar <- all.vars(formula[[2]])

  if (length(depvar) != 1)
    stop("There can only be one dependent variable.", call. = FALSE)

  rhs <- formula[[3]]

  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    runvar <- all.vars(rhs[[2]])
    controlnames <- if (identical(rhs[[3]], as.name("."))) {
      setdiff(names(data), c(depvar, runvar))
    } else {
      all.vars(rhs[[3]])
    }
    controls <- data[, controlnames, drop=FALSE]
  } else {
    runvar <- all.vars(rhs)
    controls <- NULL
  }
  if (length(runvar) != 1){
    stop("There can only be one running variable.", call. = FALSE)
  }
  return(
    list(
      depvar = data[[depvar]],
      runvar = data[[runvar]],
      controls = controls
    )
  )
}

create_data <- function(depvar=NULL, runvar=NULL, controls=NULL,
                        formula=NULL, data=NULL){
  using_direct <- !is.null(depvar) || !is.null(runvar) || !is.null(controls)
  using_formula <- !is.null(formula)

  if (using_formula && using_direct) {
    stop("'Formula' cannot be combined with 'depvar/runvar/controls'.",
         call. = FALSE)
  }

  if (!using_formula && !using_direct) {
    stop("Either 'formula' with data or 'depvar/runvar/controls' ",
         "must be supplied.", call. = FALSE)
  }

  if (using_formula && !using_direct) {
    if (is.null(data)) {
      stop("'Data' must be supplied when using 'formula'.", call. = FALSE)
    } else {
      out <- parse_formula(formula=formula, data=data)
    }
  }
  if (!using_formula && using_direct) {

    if (inherits(depvar, "formula")) {
      stop("'depvar' should not be a formula. ",
           "Perhaps you wanted to specify the keyword 'formula'?",
           call. = FALSE)
    }
    if (inherits(runvar, "formula")) {
      stop("'runvar' should not be a formula. ",
           "Perhaps you wanted to specify the keyword 'formula'?",
           call. = FALSE)
    }

    if (is.null(depvar) || is.null(runvar)){
      stop("Both dependent variable and running variable must be supplied.",
           call. = FALSE)
    }

    if (!is.character(depvar) && !is.character(runvar)){
      if (!is.null(data)){
        warning(
          "'data' was supplied but is not used because 'depvar' and 'runvar' ",
          "were supplied directly.", call. = FALSE
        )
      }
      if (!is.character(controls)){
        out <- parse_vectors(depvar, runvar, controls=as.data.frame(controls))
      } else {
        stop("'controls' can not be character when 'depvar' and 'runvar'. ",
             "are not", call. = FALSE)
      }
    }

    if (!is.null(data) && is.character(depvar) && is.character(runvar)) {
      if (is.null(controls)){
        out <- parse_vectors(data[[depvar]], data[[runvar]])
      } else {
        if (is.character(controls)){
          out <- parse_vectors(data[[depvar]], data[[runvar]],
                               controls=as.data.frame(data[controls]))
        } else {
          stop("'controls' have to be character when 'depvar' and 'runvar' ",
               "are character.", call. = FALSE)
        }
      }
    }
  }
  return(out)
}

convert_depvar <- function(depvar) {

  if (is.numeric(depvar)){
    return(depvar)
  }

  if (is.logical(depvar)){
    return(as.numeric(depvar))
  }

  stop(
    "'depvar' must be either numeric or logical.",
    call. = FALSE
  )
}

convert_runvar <- function(runvar){
  if (is.numeric(runvar)){
    fac <- factor(runvar, levels = sort(unique(runvar)))
    return(model.matrix(~ fac - 1))
  }
  stop(
    "'runvar' must be numeric.",
    call. = FALSE
  )
}

convert_controls <- function(controls){
  if (is.null(controls)){
    return(NULL)
  } else {
    f <- reformulate(names(controls))
    X_controls <- model.matrix(
      f,
      data = controls
    )[, -1, drop = FALSE]
    return(X_controls)
  }
}

sel_subset <- function(var, ind){
  if (is.null(var)){
    return(NULL)
  } else if (is.null(dim(var))){
    return(var[ind])
  } else {
    return(var[ind, , drop=FALSE])
  }
}

combine_runvar_controls <- function(xvar, cvar){
  if (is.null(cvar)){
    return(xvar)
  }
  C <- cvar
  qr_c <- qr(C)
  C <- C[, qr_c$pivot[seq_len(qr_c$rank)], drop = FALSE]
  X <- cbind(xvar, C)
  qr_X <- qr(X)
  if (ncol(X) > qr_X$rank){
    X <- X[, qr_X$pivot[seq_len(qr_X$rank)], drop = FALSE]
    warning(
      "Some control variables have been dropped because they are ",
      "collinear with levels of the running variable.", call. = FALSE
    )
  }
  nr_runvar_levels <- ncol(xvar)
  if (sum(colnames(xvar) == colnames(X)[1:nr_runvar_levels]) !=
      nr_runvar_levels){
      stop(
        "If this error message is displayed, something went ",
        "unexpectedly wrong. Please send a bug report.",
        call. = FALSE
      )
  }
  return(X)
}

in_bw <- function(runvar, cutoff, h=NULL){
  if (length(h) == 1) {
    in_bw <- (runvar >= cutoff - h) & (runvar <= cutoff + h)
  } else if (length(h) == 2){
    in_bw <- (runvar >= cutoff - h[1]) & (runvar <= cutoff + h[2])
  }
  else if (is.null(h)){
    in_bw = rep(TRUE, length(runvar))
  }
  else{
    stop("'h' option incorrectly specified. Should either be a scalar ",
         "or a vector with two values", call. = FALSE)
  }
  return(
    in_bw
  )
}

trim_data_bw <- function(depvar, runvar, cutoff, controls=NULL, h=NULL){
  if (length(h) == 1) {
    in_bw <- (runvar >= cutoff - h) & (runvar <= cutoff + h)
  } else if (length(h) == 2){
    in_bw <- (runvar >= cutoff - h[1]) & (runvar <= cutoff + h[2])
  }
  else if (is.null(h)){
    in_bw = rep(TRUE, length(runvar))
  }
  else{
    stop("'h' option incorrectly specified. Should either be a scalar ",
         "or a vector with two values", call. = FALSE)
  }
  return(
    list(
      depvar = sel_subset(depvar, in_bw),
      runvar = sel_subset(runvar, in_bw),
      controls = sel_subset(controls, in_bw)
    )
  )
}



listwise_deletion <- function(depvar, runvar, controls=NULL){
  keep_rows <- complete.cases(cbind(depvar, runvar, controls))
  return(
    list(
      depvar = sel_subset(depvar, keep_rows),
      runvar = sel_subset(runvar, keep_rows),
      controls = sel_subset(controls, keep_rows)
    )
  )
}

dd_core <- function(x, c, b, cov_mat, n=NULL, normalize=TRUE, alphaval=0.05){
  G <- length(x)
  if (normalize){
    R <- create_R_matrix(x, c)
  } else {
    R <- create_R_matrix(1:G, c)
  }
  if (is.null(n)){
    n_below <- NULL
    n_above <- NULL
  } else {
    n_below <- sum(n[1:(c-1)])
    n_above <- sum(n[c:G])
  }
  D_g <- R %*% b
  stderr_D_g <- sqrt(diag(R %*% cov_mat %*% t(R)))
  Tstats <- D_g / stderr_D_g
  T_overall <- overall_t(Tstats)
  pval <- 2 * (1-pnorm(T_overall))
  CI_VALs <- get_ci(D_g, stderr_D_g, alphaval)

  out <- list(t.stat.overall=T_overall, p.val=pval, ci.lower=CI_VALs[1],
                 ci.upper=CI_VALs[2], cov.mat=unname(cov_mat),
                 double.differences=D_g,
                 all.t.stats=Tstats, std.error.double.differences=stderr_D_g,
                 sig.level=alphaval, te.lower.bound=min(D_g),
                 te.upper.bound=max(D_g), nr.support.points.below=c-1,
                 nr.support.points.above=G-c+1, nr.obs.below=n_below,
                 nr.obs.above=n_above)
}

print.rd2diffResults <- function(x) {
  if (!is.null(x$call)) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }
  cat("Estimated bound: (", x$te.lower.bound, ", ", x$te.upper.bound, ")\n",
      100*(1-x$sig.level), "% confidence interval: (", x$ci.lower,
      ", ", x$ci.upper, ")\n",
      "T-statistic: ", x$t.stat.overall, "\n",
      "p-value: ", x$p.val, "\n",sep="")
}

rd2diff <- function(depvar=NULL, runvar=NULL, controls=NULL,
                    formula=NULL, data=NULL, cutoff=0, h=NULL,
                    normalize=TRUE, sig.level=0.05,
                    vcov = sandwich::vcovHC, vcov_args = list(type = "HC1")){
  
  vcov_defaults <- list(type = "HC1")
  vcov_args <- modifyList(
    vcov_defaults,
    vcov_args
  )
  
  cr_d <- create_data(depvar=depvar, runvar=runvar, controls=controls,
                      formula=formula, data=data)
  
  if (!is.null(vcov_args$cluster)){
    if (all.equal(vcov_args$cluster, cr_d$runvar)){
      stop("Clustering on the running variable is not allowed.", call. = FALSE)
    }
  }
  
  ind_in_bw <- in_bw(cr_d$runvar, cutoff, h=h)
  cr_d$depvar <- sel_subset(cr_d$depvar, ind_in_bw)
  cr_d$runvar <- sel_subset(cr_d$runvar, ind_in_bw)
  cr_d$controls <- sel_subset(cr_d$controls, ind_in_bw)
  
  cr_d_ld <- listwise_deletion(cr_d$depvar, cr_d$runvar, cr_d$controls)
  ng <- table(cr_d_ld$runvar)
  
  if (length(unique(cr_d$runvar[!is.na(cr_d$runvar)])) !=
      length(unique(cr_d_ld$runvar))){
    stop("There are values of the running variable with no observations",
         call. = FALSE)
  }
  
  min_ng <- min(ng)
  if (min_ng == 1){
    warning("There is a at least one support point with only a single ",
            "observation. It is highly likely that the covariance matrix ",
            "is incorrectly estimated.")
  } else if (min_ng < 20){
    warning("There is a at least one support point with fewer than 20 ",
            "observations. Be aware that the asymptotic normal approximation ",
            "may not be valid.")
  }
  
  xvar_unique <- sort(unique(cr_d_ld$runvar))
  yvar <- convert_depvar(cr_d_ld$depvar)
  xvar <- convert_runvar(cr_d_ld$runvar)
  cvar <- convert_controls(cr_d_ld$controls)
  nr_vals <- length(xvar_unique)
  xc <- min(xvar_unique[xvar_unique >= cutoff])
  g <- rank(xvar_unique)
  c <- g[xvar_unique == xc]
  
  X <- combine_runvar_controls(xvar, cvar)
  ols_model <- lm(yvar ~ X - 1)
  coefs <- coef(ols_model)[1:nr_vals]
  
  covmat <- do.call(vcov, c(list(ols_model), vcov_args))[1:nr_vals, 1:nr_vals]
  
  if (min(startsWith(names(coefs), "Xfac")) != 1){
    stop("If this error message is displayed, something went ",
         "unexpectedly wrong.",
         call. = FALSE
    )
  }
  out <- dd_core(xvar_unique, c, coefs, covmat, ng,
                 normalize=normalize, alphaval=sig.level)
  out$call=match.call()
  return(structure(out, class="rd2diffResults"))
}

rd2diff.aggregate <- function(depvar_means, runvar, cutoff=0, h=NULL,
                              std_dev_depvar_means=NULL,
                              covmat=NULL, ng=NULL, normalize=TRUE, sig.level=.05){
  if ((is.null(std_dev_depvar_means) && is.null(covmat)) |
      (!is.null(std_dev_depvar_means) && !is.null(covmat))){
    stop("One and only one of 'std_dev_depvar_means' and 'covmat' ",
         "must be supplied.", call. = FALSE)
  }
  if (length(runvar) != length(unique(runvar))){
    stop("The running variable should only contain unique values.")
  }
  if (length(runvar) != length(depvar_means)){
    stop("Dependent variable and running variable should be of equal length.")
  }
  if (!is.null(std_dev_depvar_means)){
    if (is.null(ng)){
      stop("'ng' must be supplied if 'std_dev_depvar_means' ",
           "is supplied.", call. = FALSE)
    } else {
      if (length(std_dev_depvar_means) != length(depvar_means)){
        stop("'depvar_means' and 'std_dev_depvar_means' ",
             "must have the same length.", call. = FALSE)
      }
      if (length(std_dev_depvar_means) != length(ng)){
        stop("'std_dev_depvar_means' and 'ng' ",
             "must have the same length.", call. = FALSE)
      }
      covmat <- diag(std_dev_depvar_means^2 / ng)
    }
  } else {
    dim_covmat <- dim(covmat)
    if (dim_covmat[1] != dim_covmat[2]){
      stop("'covmat' must be square.", call. =FALSE)
    } else {
      if (dim_covmat[1] != length(depvar_means)){
        stop("'depvar_means' must be the same length as each dimension",
             " of 'covmat'.", call. =FALSE)
      }
      if (!is.null(ng) && (length(depvar_means) != length(ng))){
        stop("'depvar_means' and 'ng' ",
             "must have the same length.", call. = FALSE)
      }
    }
  }
  xc <- min(runvar[runvar >= cutoff])
  g <- rank(runvar)
  c <- g[runvar == xc]
  out <- dd_core(runvar, c, depvar_means, covmat,
                 n=ng, normalize=normalize, alphaval=sig.level)
  out$call=match.call()
  return(structure(out, class="rd2diffResults"))
}
