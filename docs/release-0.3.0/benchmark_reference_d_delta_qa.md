# D/Delta full-run raw audit and post-processing recovery

## Run integrity
The original benchmark runner exited 1 after timing completed: its final pair-order aggregate referenced first_implementation, which was absent from the pair projection. Raw timed rows were saved and left unchanged. The previous pairs CSV was snapshotted before repair. No benchmark calls were repeated.
Executed runner SHA-256, from the matching immediately preceding smoke run: 440cb42d4c18ebc35863cc65e0ba771afb1cc9e918e4f9888333c4e1ab31349f
Corrected runner SHA-256, adding the omitted pair field: a70fd6c4ab980f5d00838e3587ef48f04f661bdcae6d6a95037bb03ad6ee7831
Postprocessor SHA-256: 552e4f44064c48f028d7048ea0584f58fa98636392dfc85a759ab96ac45a7af6
Raw runs SHA-256, verified unchanged: e46a67a58f3e8519b1f50ae8ed31e88005053e947949517ec647e6fe99400318
Pre-repair pair snapshot SHA-256: 96164c7e9e26ed0b720f7ca402c9b1997ba7a153a14c2e92daa803fb7c94bda9
The metadata records the original process exit code; it does not claim the benchmark runner exited 0.

## Design and provenance
Package fastphylosig 0.3.0 ; tarball SHA-256 909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff ; DLL SHA-256 1bd40080f1ea537b3c73c838dc98ff84940f1064f4e1c4dbfb2dd5063d999bdb
D: fast_d(test=TRUE, nsim=199) versus caper::phylo.d(permut=199); all fast calls report 199/199 random and 199/199 Brownian simulations.
Delta reference source SHA-256: 0dc1e987c13f0e972b0e7c3b4829b4d0fdc41ff0de3856bae4803f9ce0115175 ; source uses ape::ace(model=ARD), nentropy() and two emcmc chains.
Fast Delta actual controls and requested/successful total iterations were audited from returned fields: 10000, thin 10, burn 100, lambda0 0.1, proposal_sd 0.5, LSE, ARD, two-chain contract, total 20000.
Each timed interval measured one public call with Sys.time. Tree/trait generation, prepare_tree, comparative.data, matching, seeding and gc were excluded; Delta ACE/reconstruction and MCMC remained inside the call.
The four panels compare runtime for different configured analyses: D uses a permutation test, Delta uses two-chain MCMC, and K/lambda use test=FALSE point estimates. They do not imply equal statistical work.

## Raw audit
Raw calls: 480 (240 fast, 240 reference); pairs: 240 ; cells: 24; paired replicates/cell: 10.
Finite successful calls: 480 / 480 ; errors: 0
Observed actual species match: 240 / 240 ; all counts equal n_tip; NA removals: 0.
Tree, trait and reference-input hash matches: 240 / 240 ; 240 / 240 ; 240 / 240
First-call order balance by cell:
~~~text
 method scenario n_tip first_implementation comparison_id
      D moderate    50         fastphylosig             5
      D moderate    50            reference             5
      D moderate   100         fastphylosig             5
      D moderate   100            reference             5
      D moderate   500         fastphylosig             5
      D moderate   500            reference             5
      D moderate  2000         fastphylosig             5
      D moderate  2000            reference             5
      D     none    50         fastphylosig             5
      D     none    50            reference             5
      D     none   100         fastphylosig             5
      D     none   100            reference             5
      D     none   500         fastphylosig             5
      D     none   500            reference             5
      D     none  2000         fastphylosig             5
      D     none  2000            reference             5
      D   strong    50         fastphylosig             5
      D   strong    50            reference             5
      D   strong   100         fastphylosig             5
      D   strong   100            reference             5
      D   strong   500         fastphylosig             5
      D   strong   500            reference             5
      D   strong  2000         fastphylosig             5
      D   strong  2000            reference             5
  Delta moderate    50         fastphylosig             5
  Delta moderate    50            reference             5
  Delta moderate   100         fastphylosig             5
  Delta moderate   100            reference             5
  Delta moderate   500         fastphylosig             5
  Delta moderate   500            reference             5
  Delta moderate  2000         fastphylosig             5
  Delta moderate  2000            reference             5
  Delta     none    50         fastphylosig             5
  Delta     none    50            reference             5
  Delta     none   100         fastphylosig             5
  Delta     none   100            reference             5
  Delta     none   500         fastphylosig             5
  Delta     none   500            reference             5
  Delta     none  2000         fastphylosig             5
  Delta     none  2000            reference             5
  Delta   strong    50         fastphylosig             5
  Delta   strong    50            reference             5
  Delta   strong   100         fastphylosig             5
  Delta   strong   100            reference             5
  Delta   strong   500         fastphylosig             5
  Delta   strong   500            reference             5
  Delta   strong  2000         fastphylosig             5
  Delta   strong  2000            reference             5
~~~
Statuses (reference APIs do not expose a status field):
~~~text
, , reported_status = not_reported_by_reference, success = TRUE

       implementation
method  fastphylosig reference
  D                0       120
  Delta            0       120

, , reported_status = ok, success = TRUE

       implementation
method  fastphylosig reference
  D              120         0
  Delta          120         0

~~~

## Warning and Delta diagnostic audit
Warning-bearing calls: 113 /480; D=0; Delta reference=112/120; fast Delta=1/120. All warnings remain in raw rows; no retries or row removals.
Warning-call counts by method and implementation:
~~~text
 method implementation warning_calls
  Delta   fastphylosig             1
  Delta      reference           112
~~~
All warning-component frequencies (components may co-occur):
~~~text

                                                                                  NaNs produced 
                                                                                             96 
                                                         Inf replaced by maximum positive value 
                                                                                             70 
                                                          imaginary parts discarded in coercion 
                                                                                             38 
                                          model fit suspicious: gradients apparently non-finite 
                                                                                             11 
transition-rate optimization did not converge: iteration limit reached without convergence (10) 
                                                                                              1 
~~~
Reference Delta warning-component frequencies:
~~~text

                                        NaNs produced 
                                                   96 
               Inf replaced by maximum positive value 
                                                   70 
                imaginary parts discarded in coercion 
                                                   38 
model fit suspicious: gradients apparently non-finite 
                                                   11 
~~~
Reference Delta warning components by species count and scenario:
~~~text
, , scenario = moderate

                                                       n_tip
warning_component                                       50 100 500 2000
  Inf replaced by maximum positive value                 6   8   9    7
  NaNs produced                                          7   8   9    7
  imaginary parts discarded in coercion                  1   2   0    0
  model fit suspicious: gradients apparently non-finite  0   1   0    0

, , scenario = none

                                                       n_tip
warning_component                                       50 100 500 2000
  Inf replaced by maximum positive value                 1   0   0    0
  NaNs produced                                          6   6   6    8
  imaginary parts discarded in coercion                  6   9   9   10
  model fit suspicious: gradients apparently non-finite  2   4   2    1

, , scenario = strong

                                                       n_tip
warning_component                                       50 100 500 2000
  Inf replaced by maximum positive value                 9  10  10   10
  NaNs produced                                          9  10  10   10
  imaginary parts discarded in coercion                  1   0   0    0
  model fit suspicious: gradients apparently non-finite  1   0   0    0

~~~
The pinned reference delta() returns only a numeric scalar and discards the ACE convergence object; result_fields are <numeric>. Call-level warning capture cannot distinguish an optimizer trial warning from an accepted/final-fit warning. Finite estimates do not certify convergence.
Fast Delta chain diagnostics:
~~~text
 calls diagnostics_available ESS_alpha_min ESS_beta_min Rhat_alpha_max
   120                   120      52.40738     52.37619       1.060496
 Rhat_beta_max MCSE_Delta_max
      1.057219     0.06444899
~~~
One fast Delta warning row:
~~~text
 scenario n_tip replicate
     none  2000         6
                                                                                         warning
 transition-rate optimization did not converge: iteration limit reached without convergence (10)
  estimate ESS_alpha ESS_beta Rhat_alpha Rhat_beta   MCSE_Delta
 0.1296409  601.6414 458.0392  0.9993196  1.000341 0.0001259713
~~~
Delta elapsed times describe calls returning finite values only. They do not validate equal accuracy, effective chains, or final-fit convergence.

## Observed runtime
Ratio is median(reference seconds) / median(fast seconds), descriptive only:
~~~text
 method scenario n_tip n_pairs fast_median_seconds reference_median_seconds
      D moderate    50      10         0.007730126               0.04118252
      D moderate   100      10         0.010280013               0.06031001
      D moderate   500      10         0.032377005               0.68169451
      D moderate  2000      10         0.121750474              26.56591403
      D     none    50      10         0.007972479               0.04201508
      D     none   100      10         0.010466576               0.06171238
      D     none   500      10         0.033290982               0.66942251
      D     none  2000      10         0.122554421              26.77241051
      D   strong    50      10         0.007849932               0.04141653
      D   strong   100      10         0.010539055               0.06103456
      D   strong   500      10         0.032403469               0.68839347
      D   strong  2000      10         0.118285060              26.50673604
  Delta moderate    50      10         0.029536963               0.89077008
  Delta moderate   100      10         0.035795927               1.15015960
  Delta moderate   500      10         0.115482330               3.70165396
  Delta moderate  2000      10         0.594941020              14.21268451
  Delta     none    50      10         0.034136534               0.96909451
  Delta     none   100      10         0.042075992               1.37401056
  Delta     none   500      10         0.179730535               4.74966455
  Delta     none  2000      10         0.755549550              18.43380415
  Delta   strong    50      10         0.028530002               0.87821698
  Delta   strong   100      10         0.035148501               1.15239549
  Delta   strong   500      10         0.088249087               3.32639503
  Delta   strong  2000      10         0.335479021              11.26666796
 reference_over_fast_ratio_of_medians paired_ratio_median warning_calls_fast
                             5.327535            5.350725                  0
                             5.866725            5.885837                  0
                            21.054897           20.977776                  0
                           218.199676          218.990880                  0
                             5.270014            5.274250                  0
                             5.896139            5.896300                  0
                            20.108223           20.109586                  0
                           218.453241          218.616245                  0
                             5.276036            5.260622                  0
                             5.791275            5.786766                  0
                            21.244437           20.432458                  0
                           224.092003          231.485745                  0
                            30.157809           30.057670                  0
                            32.131019           32.925588                  0
                            32.053856           33.093828                  0
                            23.889233           24.180337                  0
                            28.388779           28.468340                  0
                            32.655453           30.216100                  0
                            26.426587           25.928093                  0
                            24.397876           22.608851                  1
                            30.782227           30.400661                  0
                            32.786476           32.758031                  0
                            37.693251           38.642811                  0
                            33.583823           33.445713                  0
 warning_calls_reference
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       0
                       8
                       9
                       9
                       7
                       9
                      10
                      10
                      10
                      10
                      10
                      10
                      10
~~~
D ratio range: 5.27 to 224.1 x.
Delta ratio range: 23.89 to 37.69 x.
The Delta figure panel must disclose reference warnings and the absence of a convergence certification.
