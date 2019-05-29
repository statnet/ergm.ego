#  File R/EgoStat.node.attr.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2019 Statnet Commons
#######################################################################
#' @name node-attr-api
#' @title Helper functions for specifying nodal attribute levels
#'
#' @description These functions are meant to be used in `EgoStat` and
#'   other implementations to provide the user with a way to extract
#'   nodal attributes and select their levels in standardized and
#'   flexible ways described under [`node-attr`]. They are intended to
#'   parallel [node-attr-api] of `ergm` package.
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

#' @rdname node-attr-api
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
#' flomego <- as.egodata(flomarriage)
#' ergm.ego_get_vattr("priorates", flomego$egos)
#' ergm.ego_get_vattr(~priorates, flomego$alters)
#' ergm.ego_get_vattr(c("wealth","priorates"), flomego$egos)
#' ergm.ego_get_vattr(~priorates>30, flomego$alters)
#' (a <- ergm.ego_get_vattr(~cut(priorates,c(-Inf,0,20,40,60,Inf),label=FALSE)-1, flomego$egos))
#' @export
ergm.ego_get_vattr <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)
  UseMethod("ergm.ego_get_vattr")
}

.handle_multiple <- function(a, multiple){
  if(!is.list(a)) a <- list(a)
  a <- do.call(cbind, a)
  if(ncol(a)>1)
    switch(multiple,
           paste =  apply(a, 1, paste, collapse="."),
           matrix = a,
           stop = ergm_Init_abort("This term does not accept multiple vertex attributes or matrix vertex attribute functions."))
  else c(a)
}

.rightsize_vattr <- function(a, df){
  rep_len_warn <- function(x, length.out){
    if(length.out%%NVL(nrow(x), length(x))) ergm_Init_warn("Length of vertex attribute vector is not a multiple of network size.")
    if(is.null(nrow(x))) rep_len(x, length.out) else apply(x, 2, rep_len, length.out)
  }
  rep_len_warn(a, nrow(df))
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
  x
}

#' @rdname node-attr-api
#' @importFrom purrr "%>%" "map" "pmap_chr"
#' @importFrom rlang set_attrs
#' @export
ergm.ego_get_vattr.character <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  missing_attr <- setdiff(object, names(df))
  if(length(missing_attr)){
    ergm_Init_abort(paste.and(sQuote(missing_attr)), " is/are not valid nodal attribute(s).")
  }

  object %>% map(~df[[.]]) %>% set_names(object) %>% .handle_multiple(multiple=multiple) %>%
    .rightsize_vattr(df) %>% set_attrs(name=paste(object, collapse=".")) %>%
    .check_acceptable(accept=accept, xspec=object)
}


#' @rdname node-attr-api
#' @export
ergm.ego_get_vattr.function <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  ERRVL(try(object(df, ...) %>%
            .rightsize_vattr(df) %>% .handle_multiple(multiple=multiple) %>%
            set_attrs(name=strtrim(despace(paste(deparse(body(object)),collapse="\n")),80)),
            silent=TRUE),
        ergm_Init_abort(.)) %>%
    .check_acceptable(accept=accept)
}


#' @rdname node-attr-api
#' @importFrom purrr "%>%" map set_names when
#' @importFrom tibble lst
#' @export
ergm.ego_get_vattr.formula <- function(object, df, accept="character", multiple=if(accept=="character") "paste" else "stop", ...){
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  a <- names(df)
  vlist <- c(a %>% map(~df[[.]]) %>% set_names(a),
             lst(`.`=df, .df=df, ...))

  e <- object[[length(object)]]
  ERRVL(try({
    eval(e, envir=vlist, enclos=environment(object)) %>%
      .rightsize_vattr(df) %>% .handle_multiple(multiple=multiple) %>%
      set_attrs(name=if(length(object)>2) eval_lhs.formula(object) else despace(paste(deparse(e),collapse="\n")))
  }, silent=TRUE),
  ergm_Init_abort(.)) %>%
    .check_acceptable(accept=accept, xspec=object)
}

#' @rdname node-attr-api
#'
#' @description `ergm.ego_attr_levels` filters the levels of the
#'   attribute.  It is strongly recommended that [check.ErgmTerm()]'s
#'   corresponding
#'   `vartype="function,formula,character,numeric,logical,AsIs,NULL"` (using the
#'   `ERGM_LEVELS_SPEC` constant).
#'
#' @param egodata An [`egodata`] object.
#' 
#' @return `ergm.ego_attr_levels` returns a vector of levels to use and their order.
#' @examples
#' ergm.ego_attr_levels(NULL, a, flomego$egos)
#' ergm.ego_attr_levels(-1, a, flomego$egos)
#' ergm.ego_attr_levels(1:2, a, flomego$egos)
#' ergm.ego_attr_levels(I(1:2), a, flomego$egos)
#' @export
ergm.ego_attr_levels <- function(object, attr, egodata, levels=sort(unique(attr)), ...){
  UseMethod("ergm.ego_attr_levels")
}

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.numeric <- function(object, attr, egodata, levels=sort(unique(attr)), ...){
  levels[object]
}

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.logical <- ergm.ego_attr_levels.numeric

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.AsIs <- function(object, attr, egodata, levels=sort(unique(attr)), ...){
  object
}

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.character <- ergm.ego_attr_levels.AsIs

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.NULL <- function(object, attr, egodata, levels=sort(unique(attr)), ...){
  levels
}

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.function <- function(object, attr, egodata, levels=sort(unique(attr)), ...){
  object <- if('...' %in% names(formals(object))) object(levels, attr, egodata, ...)
            else switch(length(formals(object)),
                        object(levels),
                        object(levels, attr),
                        object(levels, attr, egodata))
  ergm.ego_attr_levels(object, attr, egodata, levels, ...)
}

#' @rdname node-attr-api
#' @export
ergm.ego_attr_levels.formula <- function(object, attr, egodata, levels=sort(unique(attr)), ...){
  vlist <- lst(`.`=levels, .levels=levels, .attr=attr, .egodata=egodata, ...)
  e <- object[[length(object)]]
  object <- eval(e, envir=vlist, enclos=environment(object))  
  ergm.ego_attr_levels(object, attr, egodata, levels, ...)
}

