# Independent R matrix pruning and the existing ACE downward estimator. The
# fixture rate made LAPACK's general eigenvectors singular on Windows R 4.6.1.
# Other platforms may successfully decompose the same matrix: validate the
# answer rather than requiring a platform-dependent decomposition failure.
.er_recovery_oracle <- function(edge,lengths,states,k,q,marginal=FALSE) {
  n <- length(states); lik <- matrix(0,2*n-1L,k); log_scale <- 0
  for(i in seq_len(n)) if(is.na(states[i])) lik[i,] <- 1 else lik[i,states[i]] <- 1
  transition <- function(t) {
    exponent <- if(t==0) 0 else -k*q*t
    diag(exp(exponent),k)+matrix(-expm1(exponent)/k,k,k)
  }
  for(e in seq(1L,nrow(edge),by=2L)) {
    left <- as.numeric(transition(lengths[e]) %*% lik[edge[e,2L],])
    right <- as.numeric(transition(lengths[e+1L]) %*% lik[edge[e+1L,2L],])
    value <- left*right; scale <- sum(value)
    lik[edge[e,1L],] <- value/scale; log_scale <- log_scale+log(scale)
  }
  ancestors <- lik[n+seq_len(n-1L),,drop=FALSE]
  if(!marginal) for(e in seq(nrow(edge)-1L,1L,by=-2L)) {
    parent <- edge[e,1L]-n
    for(r in e:(e+1L)) {
      child <- edge[r,2L]-n
      if(child>0) {
        P <- transition(lengths[r])
        ancestors[child,] <- as.numeric((ancestors[parent,]/
          as.numeric(ancestors[child,] %*% P)) %*% P)*ancestors[child,]
      }
    }
    ancestors <- ancestors/rowSums(ancestors)
  }
  list(deviance=-2*log_scale,lik.anc=ancestors)
}

test_that("ER recovery agrees with independent pruning and ancestral updates", {
  f <- readRDS(test_path("fixtures","ace-er-singular.rds"))
  k <- nrow(f$rate_index)
  for(states in list(f$states,replace(f$states,c(3L,17L),NA_integer_))) {
    for(lengths in list(f$edge_length,replace(f$edge_length,1L,0))) {
      ptr <- fast_ace_fit_prepare_cpp(f$edge,lengths,states,f$rate_index,1L)
      set.seed(391); seed <- .Random.seed
      for(marginal in c(FALSE,TRUE)) {
        expected <- .er_recovery_oracle(f$edge,lengths,states,k,f$par,marginal)
        value <- fast_ace_fit_liks_cpp(ptr,f$par,marginal)
        expect_equal(value$deviance,expected$deviance,tolerance=1e-9)
        expect_equal(value$lik.anc,expected$lik.anc,tolerance=1e-9)
        expect_true(all(is.finite(value$lik.anc)))
        expect_equal(rowSums(value$lik.anc),rep(1,nrow(value$lik.anc)),tolerance=1e-12)
        for(reuse in c(FALSE,TRUE)) {
          fresh <- fast_ace_discrete_liks_cpp(f$edge,lengths,states,f$rate_index,
            f$par,marginal,reuse)
          expect_equal(fresh$deviance,value$deviance,tolerance=1e-10)
          expect_equal(fresh$lik.anc,value$lik.anc,tolerance=1e-10)
        }
      }
      # Returning to an ordinary eigensystem must not reuse analytic coefficients.
      expect_identical(fast_ace_fit_deviance_cpp(ptr,.1),
        fast_ace_discrete_deviance_cpp(f$edge,lengths,states,f$rate_index,.1,TRUE))
      expect_identical(.Random.seed,seed)
      expect_true(fast_ace_fit_recoveries_cpp(ptr) %in% c(0,2))
      fast_ace_fit_update_cpp(ptr,states)
      expect_identical(fast_ace_fit_recoveries_cpp(ptr),0)
      expect_identical(.Random.seed,seed)
      fast_ace_fit_release_cpp(ptr)
      expect_error(fast_ace_fit_recoveries_cpp(ptr),"released")
    }
  }
})

test_that("recovered fits report once and pooled state updates reset counts", {
  f <- readRDS(test_path("fixtures","ace-er-singular.rds"))
  k <- nrow(f$rate_index); n <- length(f$states)
  pool <- .fast_ace_native_pool(); on.exit(.fast_ace_pool_close(pool))
  w <- list(edge=f$edge,edge_length=f$edge_length,nb_tip=n,nb_node=n-1L,
    rate=.fast_ace_rate_index("ER",k),native_pool=pool,native_key=1L,fit_threads=1L)
  run <- function(workspace) {
    messages <- character(); set.seed(761); seed <- .Random.seed
    fit <- withCallingHandlers(.fast_ace_fit_workspace(f$states,letters[seq_len(k)],
      workspace,ip=f$par,CI=TRUE,estimate_se=FALSE),warning=function(c) {
        messages <<- c(messages,conditionMessage(c)); invokeRestart("muffleWarning")
      })
    expect_identical(.Random.seed,seed)
    recovered <- grepl("analytic ER transitions",messages,fixed=TRUE)
    expect_lte(sum(recovered),1L)
    expect_identical(any(recovered),!is.null(fit$numerical_recovery))
    if(any(recovered)) expect_gt(fit$numerical_recovery$evaluations,0)
    list(fit=fit,messages=messages)
  }
  first <- run(w); expect_identical(run(w),first)
  w$native_pool <- NULL; w$native_key <- NULL
  expect_identical(run(w),first)
})
