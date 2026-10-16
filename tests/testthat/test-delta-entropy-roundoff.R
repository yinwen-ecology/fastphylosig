.entropy_roundoff_raw <- function(p,type) {
  k <- ncol(p); out <- numeric(nrow(p))
  for(i in seq_len(nrow(p))) {
    total <- 0
    for(j in seq_len(k)) {
      value <- p[i,j]
      if(type==1L) {
        if(value>1/k)value <- value/(1-k)-1/(1-k)
        total <- total+value
      } else if(type==2L) {
        if(value>0)total <- total-value*log(value)
      } else total <- total+value*value
    }
    out[i] <- if(type==2L)total/log(k) else if(type==3L)((1-total)*k)/(k-1) else total
  }
  out
}

.entropy_roundoff_jitter <- function(x,k) {
  tol <- min(1e-12,8*k*.Machine$double.eps)
  x[x<0 & x>=-tol] <- 0; x[x>1 & x<=1+tol] <- 1
  if(any(x==0)) x[x==0] <- x[x==0]+runif(1)/10000
  if(any(x==1)) x[x==1] <- x[x==1]-runif(1)/10000
  x
}

test_that("real Shannon roundoff input continues with the existing jitter", {
  old_rng <- RNGkind()
  on.exit(do.call(RNGkind,as.list(old_rng)),add=TRUE)
  p <- readRDS(test_path("fixtures","delta-entropy-roundoff.rds"))
  raw <- .entropy_roundoff_raw(p,2L)
  for(uniform in c("Mersenne-Twister","L'Ecuyer-CMRG"))
    for(normal in c("Inversion","Box-Muller","Kinderman-Ramage")) {
      RNGkind(uniform,normal,"Rejection"); set.seed(617)
      expected <- .entropy_roundoff_jitter(raw,ncol(p)); seed <- .Random.seed
      set.seed(617); actual <- suppressWarnings(delta_entropy_cpp(p,2L))
      expect_equal(as.numeric(actual),expected,tolerance=1e-14)
      expect_identical(.Random.seed,seed)
      expect_true(all(actual>0 & actual<1))
      set.seed(617)
      fit <- suppressWarnings(delta_mcmc_cpp(p,.1,.5,500L,5L,30L,2L,TRUE))
      seed <- .Random.seed
      expect_true(is.finite(fit$delta) && fit$delta>0)
      if(!is.null(fit$entropy_recovery)) {
        expect_gt(fit$entropy_recovery$above,0)
        expect_lte(fit$entropy_recovery$maximum_correction,fit$entropy_recovery$tolerance)
      }
      set.seed(617); x <- .entropy_roundoff_jitter(raw,ncol(p))
      reference <- .delta_posterior_reference(x,.1,.5,500L,5L,30L)
      expect_equal(fit$alpha_chain,reference$alpha,tolerance=1e-8)
      expect_equal(fit$beta_chain,reference$beta,tolerance=1e-8)
      expect_equal(fit$delta,mean(reference$beta/reference$alpha),tolerance=1e-9)
      expect_identical(.Random.seed,seed)
    }
})

test_that("roundoff repair shares endpoint draws and reports only recovered rows", {
  # Gini arithmetic guarantees a few-ulp upper excursion independently of log.
  p <- rbind(c(.5,.5)-.Machine$double.eps,c(.5,.5),c(1,0),c(.8,.2))
  raw <- .entropy_roundoff_raw(p,3L)
  expect_gt(raw[1L],1)
  set.seed(711); expected <- .entropy_roundoff_jitter(raw,2L); seed <- .Random.seed
  set.seed(711)
  expect_warning(actual <- delta_entropy_cpp(p,3L),"roundoff in 1 row")
  expect_equal(as.numeric(actual),expected,tolerance=1e-15)
  expect_identical(.Random.seed,seed)
  expect_identical(actual[1L],actual[2L])
  set.seed(712)
  expect_warning(fit <- delta_mcmc_cpp(p,.1,.5,80L,2L,10L,3L),"roundoff in 1 row")
  expect_identical(fit$entropy_recovery$above,1L)
  expect_identical(fit$entropy_recovery$below,0L)
  expect_equal(fit$entropy_recovery$maximum_correction,raw[1L]-1,tolerance=1e-16)
  expect_identical(fit$entropy_recovery$tolerance,16*.Machine$double.eps)
  set.seed(713)
  healthy <- delta_mcmc_cpp(p[-1L,,drop=FALSE],.1,.5,80L,2L,10L,3L)
  expect_null(healthy$entropy_recovery)
})

test_that("material or invalid-row excursions are errors without jitter draws", {
  eps <- .Machine$double.eps
  inputs <- list(c(.5,.5)-1e-8,c(-eps,1+eps),c(1,1),c(-.1,1.1))
  for(value in inputs) {
    p <- rbind(c(1,0),value)
    set.seed(716); seed <- .Random.seed
    expect_error(delta_mcmc_cpp(p,.1,.5,80L,2L,10L,3L),"cannot be repaired")
    expect_identical(.Random.seed,seed)
  }
  # Small entropy excursions must also have a normalized input row.
  p <- matrix(rep(.5-2*.Machine$double.eps,2),nrow=1L)
  raw <- .entropy_roundoff_raw(p,3L)
  expect_lt(raw[1L]-1,16*eps)
  # A more extreme row cannot be repaired even if entropy is finite.
  expect_error(delta_entropy_cpp(matrix(c(.1,.1),1L),3L),"cannot be repaired")
})
