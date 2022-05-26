#  File R/EgoStat.node.attr.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################
#' @name nodal_attributes-API
#' @title Helper functions for specifying nodal attribute levels
#'
#' @description These functions are meant to be used in `EgoStat` and other
#'   implementations to provide the user with a way to extract nodal attributes
#'   and select their levels in standardized and flexible ways. They are
#'   intended to parallel [ergm::nodal_attributes-API] of `ergm` package.
#'
#' @param object An argument specifying the nodal attribute to select
#'   or which levels to include.
#' @param df Table of egos or of alters.
#' @param attr A vector of length equal to the number of nodes,
#'   specifying the attribute vector.
#' @param levels Starting set of levels to use; defaults to the sorted
#'   list of unique attributes.
#' @param multiple Handling of multiple attributes or matrix or data
#'   frame output. See the Details section for the specification.
#' @param accept A character vector listing permitted data types for
#'   the output. See the Details section for the specification.
#' @param ... Additional argument to the functions of network or to
#'   the formula's environment.
#'
#' @details The `accept` argument is meant to allow the user to
#'   quickly check whether the output is of an *acceptable* class or
#'   mode. Typically, if a term accepts a character (i.e.,
#'   categorical) attribute, it will also accept a numeric one,
#'   treating each number as a category label. For this reason, the
#'   following outputs are defined:
#' \describe{
#'
#' \item{`"character"`}{Accept any mode or class (since it can
#' beconverted to character).}
#' 
#' \item{`"numeric"`}{Accept real, integer, or logical.}
#' 
#' \item{`"logical"`}{Accept logical.}
#' 
#' \item{`"integer"`}{Accept integer or logical.}
#' 
#' \item{`"natural"`}{Accept a strictly positive integer.}
#' 
#' \item{`"0natural"`}{Accept a nonnegative integer or logical.}
#' 
#' \item{`"nonnegative"`}{Accept a nonnegative number or logical.}
#'
#' \item{`"positive"`}{Accept a strictly positive number or logical.}
#' }
#'
#' \describe{
#' 
#' \item{`"paste"`}{Paste together with dot as the separator.}
#' 
#' \item{`"stop"`}{Fail with an error message.}
#'
#' \item{`"matrix"`}{Construct and/or return a matrix whose rows correspond to vertices.}
#'
#' }
#'
#'
NULL

#' @rdname nodal_attributes-API
#'
#' @description `ergm.ego_get_vattr` extracts and processes the specified
#'   nodal attribute vector. It is strongly recommended that
#'   [check.ErgmTerm()]'s corresponding
#'   `vartype="function,formula,character"` (using the
#'   `ERGM_VATTR_SPEC` constant).
#' 
#' @return `ergm.ego_get_vattr` returns a vector of length equal to the number of nodes giving the
#'   selected attribute function. It may also have an attribute
#'   `"name"`, which controls the suggested name of the attribute
#'   combination.
#'
#' @examples
#' data(florentine)
#' flomego <- as.egor(flomarriage)
#' ergm.ego_get_vattr("priorates", flomego)
#' ergm.ego_get_vattr(~priorates, flomego)
#' ergm.ego_get_vattr(c("wealth","priorates"), flomego)
#' ergm.ego_get_vattr(~priorates>30, flomego)
#' (a <- ergm.ego_get_vattr(~cut(priorates,c(-Inf,0,20,40,60,Inf),label=FALSE)-1, flomego))
#' @export
ergm.ego_get_vattr <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)
  UseMethod("ergm.ego_get_vattr")
}

.handle_multiple <- function(a, multiple){
  name <- attr(a, "name")
  if(!is.list(a)) a <- list(a)
  a <- do.call(cbind, a)
  structure(
    if(ncol(a)>1)
      switch(multiple,
             paste =  apply(a, 1, paste, collapse="."),
             matrix = a,
             stop = ergm_Init_abort("This term does not accept multiple vertex attributes or matrix vertex attribute functions."))
    else c(a),
    name = name)
}

.rightsize_vattr <- function(a, df){
  name <- attr(a, "name")
  rep_len_warn <- function(x, length.out){
    if(length.out%%NVL(nrow(x), length(x))) ergm_Init_warn("Length of vertex attribute vector is not a multiple of network size.")
    if(is.null(nrow(x))) rep_len(x, length.out) else apply(x, 2, rep_len, length.out)
  }
  structure(rep_len_warn(a, nrow(df)), name=name)
}

.check_acceptable <- function(x, accept=c("character", "numeric", "logical", "integer", "natural", "0natural", "nonnegative"), xspec=NULL){
  accept <- match.arg(accept)

  ACCNAME <- list(character = "a character",
                  logical = "a logical",
                  numeric = "a numeric or logical",
                  integer = "an integer or logical",
                  natural = "a natural (positive integer) numeric",
                  `0natural` = "a nonnegative integer or logical",
                  nonnegative = "a nonnegative numeric or logical",
                  positive = "a positive numeric or logical")
  OK <-
    if(accept == "character") TRUE
    else if(!is.numeric(x) && !is.logical(x)) FALSE
    else switch(accept,
                numeric = TRUE,
                logical = all(x %in% c(FALSE, TRUE)),
                integer = all(round(x)==x),
                natural = all(round(x)==x) && x>0,
                `0natural` = all(round(x)==x) && x>=0,
                nonnegative = x>=0,
                positive = x>0)

  if(!OK) ergm_Init_abort("Attribute ", NVL3(xspec, paste0(sQuote(paste(deparse(.),collapse="\n")), " ")), "is not ", ACCNAME[[accept]], " vector as required.")
  if(is.matrix(x) && !is.null(cn <- colnames(x))){
    if(any(cn=="")){
      ergm_Init_warn("Attribute specification ", NVL3(xspec, paste0(sQuote(paste(deparse(.),collapse="\n")), " ")), "is a matrix with some column names set and others not; you may need to set them manually. See example(nodal_attributes) for more information.")
      colnames(x) <- NULL
    }
  }
  x
}

#' @rdname nodal_attributes-API
#' @importFrom purrr "%>%" "map" "pmap_chr"
#' @export
ergm.ego_get_vattr.character <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)
  df <- as_tibble(df)

  missing_attr <- setdiff(object, names(df))
  if(length(missing_attr)){
    ergm_Init_abort(paste.and(sQuote(missing_attr)), " is/are not valid nodal attribute(s).")
  }

  object %>% map(~df[[.]]) %>% set_names(object) %>% .handle_multiple(multiple=multiple) %>%
    .rightsize_vattr(df) %>% structure(name=paste(object, collapse=".")) %>%
    .check_acceptable(accept=accept, xspec=object)
}


#' @rdname nodal_attributes-API
#' @export
ergm.ego_get_vattr.function <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)
  df <- as_tibble(df)

  args <- list()
  for(aname in c("accept", "multiple"))
    if('...' %in% names(formals(object)) || aname %in% names(formals(object)))
      args[[aname]] <- get(aname)
  args <- c(list(df), list(...), args)

  ERRVL(try({
    a <- do.call(object, args)
    while(is(a,'formula')||is(a,'function')) a <- ergm.ego_get_vattr(a, df, accept=accept, multiple=multiple, ...)
    a %>% .rightsize_vattr(df) %>% .handle_multiple(multiple=multiple) %>%
      structure(., name=NVL(attr(.,"name"), strtrim(despace(paste(deparse(body(object)),collapse="\n")),80)))
  }, silent=TRUE),
  ergm_Init_abort(.)) %>%
    .check_acceptable(accept=accept)
}


#' @rdname nodal_attributes-API
#' @importFrom purrr "%>%" map set_names when
#' @importFrom tibble lst
#' @export
ergm.ego_get_vattr.formula <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)
  df <- as_tibble(df)

  a <- names(df)
  vlist <- c(a %>% map(~df[[.]]) %>% set_names(a),
             lst(`.`=df, .df=df, ...))

  e <- ult(object)
  ERRVL(try({
    a <- eval(e, envir=vlist, enclos=environment(object))
    while(is(a,'formula')||is(a,'function')) a <- ergm.ego_get_vattr(a, df, accept=accept, multiple=multiple, ...)
      a %>% .rightsize_vattr(df) %>% .handle_multiple(multiple=multiple) %>%
      structure(., name=NVL(attr(.,"name"), if(length(object)>2) eval_lhs.formula(object) else despace(paste(deparse(e),collapse="\n"))))
  }, silent=TRUE),
  ergm_Init_abort(.)) %>%
    .check_acceptable(accept=accept, xspec=object)
}

#' @rdname nodal_attributes-API
#'
#' @description `ergm.ego_attr_levels` filters the levels of the
#'   attribute.  It is strongly recommended that [check.ErgmTerm()]'s
#'   corresponding
#'   `vartype="function,formula,character,numeric,logical,AsIs,NULL"` (using the
#'   `ERGM_LEVELS_SPEC` constant).
#'
#' @param egor An [`egor`] object.
#' 
#' @return `ergm.ego_attr_levels` returns a vector of levels to use and their order.
#' @examples
#' ergm.ego_attr_levels(NULL, a, flomego)
#' ergm.ego_attr_levels(-1, a, flomego)
#' ergm.ego_attr_levels(1:2, a, flomego)
#' ergm.ego_attr_levels(I(1:2), a, flomego)
#' @export
ergm.ego_attr_levels <- function(object, attr, egor, levels=sort(unique(attr)), ...){
  UseMethod("ergm.ego_attr_levels")
}

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.numeric <- function(object, attr, egor, levels=sort(unique(attr)), ...){
  levels[object]
}

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.logical <- ergm.ego_attr_levels.numeric

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.AsIs <- function(object, attr, egor, levels=sort(unique(attr)), ...){
  object
}

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.character <- ergm.ego_attr_levels.AsIs

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.NULL <- function(object, attr, egor, levels=sort(unique(attr)), ...){
  levels
}

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.matrix <- function(object, attr, egor, levels=sort(unique(attr)), ...){

  # This should get the levels in the right order.
  ol <- levels %>% map(1L) %>% unique
  nol <- length(ol)
  il <- levels %>% map(2L) %>% unique
  nil <- length(il)

  # Construct a matrix indicating where on the levels list does each
  # element go. Then, indexing elements of m with either a logical
  # matrix or a two-column matrix of cell indices will produce a list
  # of level indices selected along with 0s, which can then be
  # dropped.
  ol2c <- match(levels%>%map(1L), ol)
  il2c <- match(levels%>%map(2L), il)
  m <- matrix(0L, nol, nil)
  m[cbind(ol2c,il2c)] <- seq_along(levels)

  sel <- switch(mode(object),
                logical = { # Binary matrix
                  if(any(dim(object)!=c(nol,nil))) ergm_Init_abort("Level combination selection binary matrix should have dimension ", nol, " by ", nil, " but has dimension ", nrow(object), " by ", ncol(object), ".") # Check dimension.
                  if(identical(ol,il)) object <- object | t(object) # Symmetrize, if appropriate.
                  object
                },
                numeric = { # Two-column index matrix
                  if(ncol(object)!=2) ergm_Init_abort("Level combination selection two-column index matrix should have two columns but has ", ncol(object), ".")
                  if(identical(ol,il)) object <- rbind(object, object[,2:1,drop=FALSE]) # Symmetrize, if appropriate.
                  object
                },
                ergm_Init_abort("Level combination selection matrix must be either numeric or logical.")
                )

  sel <- m[sel] %>% keep(`!=`,0L) %>% sort %>% unique
  levels[sel]
}

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.function <- function(object, attr, egor, levels=sort(unique(attr)), ...){
  object <- if('...' %in% names(formals(object))) object(levels, attr, egor, ...)
            else switch(length(formals(object)),
                        object(levels),
                        object(levels, attr),
                        object(levels, attr, egor))
  ergm.ego_attr_levels(object, attr, egor, levels, ...)
}

#' @rdname nodal_attributes-API
#' @export
ergm.ego_attr_levels.formula <- function(object, attr, egor, levels=sort(unique(attr)), ...){
  vlist <- lst(`.`=levels, .levels=levels, .attr=attr, .egor=egor, ...)
  e <- ult(object)
  object <- eval(e, envir=vlist, enclos=environment(object))  
  ergm.ego_attr_levels(object, attr, egor, levels, ...)
}

#' @describeIn nodal_attributes-API
#' A version of [ergm::COLLAPSE_SMALLEST()] that can handle both [`network`] and [`egodata`] objects.
#'
#' @param n,into see [ergm::COLLAPSE_SMALLEST()].
#'
#' @export
COLLAPSE_SMALLEST <- function(object, n, into){
  attr <- object
  function(.x, ...){
    vattr <- if(is.network(.x)) ergm_get_vattr(attr, .x, ...)
             else if(is.data.frame(.x)){
               ergm_Init_warn(paste(sQuote("COLLAPSE_SMALLEST()"), " may behave unpredictably with egocentric data and is not recommended at this time."))
               ergm.ego_get_vattr(attr, .x, ...)
             }else stop("Unrecognised data type. This indicates a bug.")
    lvls <- unique(vattr)
    vattr.codes <- match(vattr,lvls)
    smallest <- which(order(tabulate(vattr.codes), decreasing=FALSE)<=n)
    vattr[vattr.codes %in% smallest] <- into
    vattr
  }
}
