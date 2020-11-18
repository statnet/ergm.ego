#  File R/locator.R in package statnet.common, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

locate_function_cache <- local({
  cache <- list()
  watchlist <- character(0) # Packages being watched for unloading.
  pkglist <- character(0) # Current list of packages.
  # Reset the cache and update the list of watched packages.
  reset <- function(...){
    pkglist <<- .packages()
    new <- setdiff(pkglist, watchlist)
    for(pkg in new){
      setHook(packageEvent(pkg, "detach"), reset)
      setHook(packageEvent(pkg, "onUnload"), reset)
    }
    watchlist <<- c(watchlist, new)
    cache <<- list()
  }
  # Check if new namespaces have been added.
  checknew <- function(){
    if(!setequal(.packages(), pkglist)) reset()
  }
  function(name, env=NULL){
    checknew()
    if(is.null(env)){
      cache[[name]]
    }else{
      cache[[name]] <<- env
    }
  }
})

locate_function <- function(name, env = globalenv(), ...){
  if(is.call(name)) name <- name[[1]]
  name <- as.character(name)
  
  # Try the given environment...
  if(!is.null(obj<-get0(name, mode='function', envir=env))){
    env <- environment(obj)
    envname <- environmentName(env)
    # Check that environment name is not blank or globalenv(), and
    # that the detected environment actually contains the object.
    if(! NVL(envname,"") %in% c("", "R_GlobalEnv") && exists(name, mode='function', envir=env, inherits=FALSE)) return(call(":::",as.name(envname),as.name(name)))
    else return(as.name(name))
  }

  # Try the cache...
  envname <- locate_function_cache(name)
  if(!is.null(envname)) return(call(":::",as.name(envname),as.name(name)))

  # Use getAnywhere()...
  #' @importFrom utils getAnywhere
  m <- getAnywhere(name)
  if(length(m$objs)){
    ## Prioritise visible over not:
    if(any(m$visible)){
      m <- lapply(m[-1], "[", m$visible)
    }
    if(length(m$objs)>1) warning("Name ",name," matched by multiple objects; using the first one on the list.", ...)
    envname <- environmentName(environment(m$objs[[1]]))
    locate_function_cache(name, envname)
    return(call(":::",as.name(envname),as.name(name)))
  }
  NULL
}

locate_prefixed_function <- function(name, prefix, errname, env = globalenv(), ..., call.=FALSE){
  if(is.call(name)) name <- name[[1]]
  name <- as.character(name)
  fname <- paste(prefix,name,sep=".")
  f <- locate_function(fname, env, ...)
  if(is.null(f) && !is.null(errname)) stop(errname,' ', sQuote(name), " function ", sQuote(fname), " not found.", ..., call.=call.)
  else f
}
