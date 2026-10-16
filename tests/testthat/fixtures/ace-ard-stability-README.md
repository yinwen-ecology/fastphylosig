Four 64-tip, five-state ARD failures from the frozen independent calibration
`docs/signal-heterogeneity/ard-depth-calibration-performance-20261015`.
Keys 5:42 and 5:51 have ill-conditioned eigensystems and incorrect old native
objectives; keys 5:20 and 5:68 have roundoff outside probability boundaries.
The stored objective is from independent `expm::expm` transitions and base-R
scaled pruning at the stored rates, not a new fitted objective. Inputs are
regression fixtures, excluded from subsequent independent calibration.
