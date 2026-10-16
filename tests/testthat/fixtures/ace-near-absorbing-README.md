Development regression from the failed 0.3.3.9010 independent calibration:
case8/run36, second captured bounded ARD retry, 64 tips and three states.
The original optimizer return and selected parameters/workspace/states are saved.

The independent deviance at the saved selected parameters is
98.5522241033134354379055884762753643383395791796919788744937200729738612029982932715231952.
It was calculated from the exact IEEE double inputs with Python Decimal,
Taylor matrix exponential with scaling/squaring, and independent pruning;
60 and 90 digit runs agree. Default expm, Ward77 and PadeRBS agree as well.
The unguarded Windows spectral deviance was 98.552217306507657, exceeding the
predefined 1e-6 tolerance. No threshold was changed or sample removed.

Reproducible inputs, scripts and full failure diagnostics are preserved in
docs/signal-heterogeneity/ard-local-recovery-20261015/oracle-exception-* and
oracle-decimal.py. This case is development evidence for 9011, not fresh
independent validation of the repair.
