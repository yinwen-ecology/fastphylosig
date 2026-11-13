# Shared inputs and timing normalization for the frozen public K contract.
.compat_k_input <- function() {
  list(tree = ape::read.tree(text = "((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9);"),
       x = cbind(first = c(a=1.2,b=-0.7,c=2.3,d=0.4),
                 second = c(a=-2,b=.3,c=1.1,d=3.5)),
       permutations = rbind(c(1L,2L,3L,4L), c(2L,1L,4L,3L),
                            c(3L,4L,1L,2L), c(4L,3L,2L,1L)))
}
.compat_k_normalize <- function(fit) {
  # Preserve timing schema/status/units, replacing only measured durations.
  timing <- attr(fit, "timing")
  timing$total_elapsed[] <- 0
  timing$stages[] <- 0
  attr(fit, "timing") <- timing
  workflow <- attr(fit, "workflow")
  workflow$total_elapsed[] <- 0
  attr(fit, "workflow") <- workflow
  if (!is.data.frame(fit) && !is.null(fit$timing)) {
    fit$timing$total_elapsed[] <- 0
    fit$timing$stages[] <- 0
  }
  fit
}
.compat_k_capture <- function() {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE))
      rm(".Random.seed", envir=.GlobalEnv)
  })
  # Pin the complete R 4.6.1 RNG configuration used by this frozen fixture.
  # R-devel adds a binomial-generator code to .Random.seed[1]. K uses no
  # binomial draws, but restoring the historical configuration keeps the
  # full seed comparison exact rather than discarding any state metadata.
  withCallingHandlers(RNGversion("4.6.1"), warning = function(w) {
    if (identical(conditionMessage(w), "Buggy BTPE algorithm used for rbinom()"))
      invokeRestart("muffleWarning")
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  input <- .compat_k_input()
  out <- list()
  for (vector in c(FALSE, TRUE)) for (explicit in c(FALSE, TRUE)) {
    set.seed(424243L)
    before <- .Random.seed
    fit <- suppressMessages(fastphylosig::fast_k(
      input$tree, if (vector) input$x[,1L] else input$x,
      test=TRUE, nsim=if (explicit) 4L else 31L,
      permutations=if (explicit) input$permutations else NULL,
      return_sim=TRUE, verbose=FALSE, progress=FALSE, ncores=1L))
    out[[paste(vector, explicit)]] <- list(
      result=.compat_k_normalize(fit), before=before, after=.Random.seed)
  }
  out
}
