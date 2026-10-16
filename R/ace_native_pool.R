# Analysis-owned native scratch. No cache enters a prepared tree or result.
.fast_ace_native_pool <- function() {
  budget <- getOption("fastphylosig.ace_cache_bytes",64*1024^2)
  if (!is.numeric(budget) || length(budget)!=1L || !is.finite(budget) || budget<0)
    stop("fastphylosig.ace_cache_bytes must be a finite non-negative number.",call.=FALSE)
  pool <- new.env(parent=emptyenv())
  pool$budget <- budget;pool$pid <- Sys.getpid();pool$clock <- 0;pool$next_id <- 0
  pool$bytes <- 0
  pool$entries <- new.env(parent=emptyenv())
  reg.finalizer(pool,.fast_ace_pool_close,onexit=TRUE)
  pool
}

.fast_ace_pool_close <- function(pool) {
  if(is.null(pool))return(invisible(NULL))
  keys <- ls(pool$entries,all.names=TRUE)
  for(key in keys)try(fast_ace_fit_release_cpp(pool$entries[[key]]$pointer),silent=TRUE)
  if(length(keys))rm(list=keys,envir=pool$entries)
  pool$bytes <- 0
  invisible(NULL)
}

.fast_ace_pool_trim <- function(pool,protect=NULL) {
  if(pool$bytes<=pool$budget)return(invisible(NULL))
  keys <- ls(pool$entries,all.names=TRUE)
  if(!length(keys))return(invisible(NULL))
  entries <- mget(keys,pool$entries,inherits=FALSE)
  for(i in order(vapply(entries,`[[`,numeric(1),"used"))) {
    if(pool$bytes<=pool$budget)break
    key <- keys[i];if(identical(key,protect))next
    fast_ace_fit_release_cpp(entries[[i]]$pointer)
    rm(list=key,envir=pool$entries);pool$bytes <- pool$bytes-entries[[i]]$bytes
  }
  invisible(NULL)
}

.fast_ace_pool_acquire <- function(workspace,states,threads) {
  pool <- workspace$native_pool
  key <- paste(workspace$native_key,threads,sep=":")
  entry <- pool$entries[[key]]
  if(is.null(entry)) {
    ptr <- fast_ace_fit_prepare_cpp(workspace$edge,workspace$edge_length,states,
      workspace$rate$rate_index,threads)
  } else {
    ptr <- entry$pointer;fast_ace_fit_update_cpp(ptr,states)
  }
  pool$clock <- pool$clock+1
  bytes <- fast_ace_fit_bytes_cpp(ptr)
  pool$bytes <- pool$bytes+bytes-if(is.null(entry))0 else entry$bytes
  pool$entries[[key]] <- list(pointer=ptr,used=pool$clock,bytes=bytes)
  .fast_ace_pool_trim(pool,protect=key)
  list(pointer=ptr,key=key)
}

.fast_ace_pool_account <- function(pool,key) {
  entry <- pool$entries[[key]]
  if(is.null(entry))return(invisible(NULL))
  bytes <- fast_ace_fit_bytes_cpp(entry$pointer)
  pool$bytes <- pool$bytes+bytes-entry$bytes
  entry$bytes <- bytes
  pool$entries[[key]] <- entry
  .fast_ace_pool_trim(pool)
}
