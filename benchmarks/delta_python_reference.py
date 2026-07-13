"""Benchmark wrapper for the original Delta Python implementation.

The Delta functions below are a direct, formatting-only transcription of
diogo-s-ribeiro/delta-statistic/Delta-Python/delta_functs.py.  The command
line wrapper is local benchmark glue: it loads an ancestral-probability matrix,
runs one first call to include Numba compilation, then runs warmed calls.
"""

import argparse
import csv
import math
import time

import numpy as np
from numba import float64, int64, njit
from scipy.stats import entropy


@njit(float64(float64, float64, float64[::1], float64, float64))
def mhalpha(a, b, x, l0, se):
    a1 = np.exp(np.random.normal(np.log(a), se, 1))[0]
    lp_a = np.exp(
        (
            len(x) * (math.lgamma(a1 + b) - math.lgamma(a1))
            - a1 * (l0 - np.sum(np.log(x)))
        )
        - (
            len(x) * (math.lgamma(a + b) - math.lgamma(a))
            - a * (l0 - np.sum(np.log(x)))
        )
    )
    r = min(1, lp_a)
    while np.isnan(lp_a) == True:
        a1 = np.exp(np.random.normal(np.log(a), se, 1))[0]
        lp_a = np.exp(
            (
                len(x) * (math.lgamma(a1 + b) - math.lgamma(a1))
                - a1 * (l0 - np.sum(np.log(x)))
            )
            - (
                len(x) * (math.lgamma(a + b) - math.lgamma(a))
                - a * (l0 - np.sum(np.log(x)))
            )
        )
        r = min(1, lp_a)
    if np.random.uniform(0, 1) < r:
        return a1
    return a


@njit(float64(float64, float64, float64[::1], float64, float64))
def mhbeta(a, b, x, l0, se):
    b1 = np.exp(np.random.normal(np.log(b), se, 1))[0]
    lp_b = np.exp(
        (
            len(x) * (math.lgamma(a + b1) - math.lgamma(b1))
            - b1 * (l0 - np.sum(np.log(1 - x)))
        )
        - (
            len(x) * (math.lgamma(a + b) - math.lgamma(b))
            - b * (l0 - np.sum(np.log(1 - x)))
        )
    )
    r = min(1, lp_b)
    while np.isnan(lp_b) == True:
        b1 = np.exp(np.random.normal(np.log(b), se, 1))[0]
        lp_b = np.exp(
            (
                len(x) * (math.lgamma(a + b1) - math.lgamma(b1))
                - b1 * (l0 - np.sum(np.log(1 - x)))
            )
            - (
                len(x) * (math.lgamma(a + b) - math.lgamma(b))
                - b * (l0 - np.sum(np.log(1 - x)))
            )
        )
        r = min(1, lp_b)
    if np.random.uniform(0, 1) < r:
        return b1
    return b


@njit(
    float64[:, ::1](
        float64, float64, float64[::1], float64, float64, int64, int64, int64
    )
)
def emcmc(alpha, beta, x, l0, se, sim, thin, burn):
    n_size = np.linspace(burn, sim, int((sim - burn) / thin + 1))
    usim = np.round(n_size, 0, np.empty_like(n_size))
    gibbs = []
    p = 0
    for i in range(sim + 1):
        alpha = mhalpha(alpha, beta, x, l0, se)
        beta = mhbeta(alpha, beta, x, l0, se)
        if i == usim[p]:
            gibbs.append((alpha, beta))
            p += 1
    gibbs = np.asarray(gibbs)
    return gibbs


def entropy_type(prob, ent_type):
    if ent_type == "LSE":
        k = np.shape(prob)[1]
        prob = np.asarray(
            np.where(prob <= (1 / k), prob, prob / (1 - k) - 1 / (1 - k))
        )
        tent = np.sum(prob, 1)
        tent = np.asarray(
            np.where(tent != 0, tent, tent + np.random.uniform(0, 1, 1) / 10000)
        )
        tent = np.asarray(
            np.where(tent != 1, tent, tent - np.random.uniform(0, 1, 1) / 10000)
        )
        return tent
    elif ent_type == "SE":
        k = np.shape(prob)[1]
        tent = entropy(prob, base=k, axis=1)
        tent = np.asarray(
            np.where(tent != 0, tent, tent + np.random.uniform(0, 1, 1) / 10000)
        )
        tent = np.asarray(
            np.where(tent != 1, tent, tent - np.random.uniform(0, 1, 1) / 10000)
        )
        return tent
    else:
        k = np.shape(prob)[1]
        tent = ((1 - np.sum(prob**2, axis=1)) * k) / (k - 1)
        tent = np.asarray(
            np.where(tent != 0, tent, tent + np.random.uniform(0, 1, 1) / 10000)
        )
        tent = np.asarray(
            np.where(tent != 1, tent, tent - np.random.uniform(0, 1, 1) / 10000)
        )
        return tent


def delta(x, lambda0, se, sim, thin, burn, ent_type):
    mc1 = emcmc(
        np.random.exponential(),
        np.random.exponential(),
        entropy_type(x, ent_type),
        lambda0,
        se,
        sim,
        thin,
        burn,
    )
    mc2 = emcmc(
        np.random.exponential(),
        np.random.exponential(),
        entropy_type(x, ent_type),
        lambda0,
        se,
        sim,
        thin,
        burn,
    )
    mchain = np.concatenate((mc1, mc2), axis=0)
    deltaA = (np.mean(mchain[:, 1])) / (np.mean(mchain[:, 0]))
    return deltaA


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prob", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--entropy-out", required=True)
    parser.add_argument("--lambda0", type=float, default=0.1)
    parser.add_argument("--se", type=float, default=0.5)
    parser.add_argument("--sim", type=int, default=10000)
    parser.add_argument("--thin", type=int, default=10)
    parser.add_argument("--burn", type=int, default=100)
    parser.add_argument("--entropy", default="LSE")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    np.random.seed(args.seed)
    prob = np.loadtxt(args.prob, delimiter=",", dtype=np.float64)
    prob = np.ascontiguousarray(prob)

    ent = entropy_type(prob.copy(), args.entropy)
    np.savetxt(args.entropy_out, ent, delimiter=",")

    t0 = time.perf_counter()
    first_delta = delta(
        prob.copy(), args.lambda0, args.se, args.sim, args.thin, args.burn,
        args.entropy
    )
    first_time = time.perf_counter() - t0

    rows = []
    for i in range(args.repeats):
        t0 = time.perf_counter()
        value = delta(
            prob.copy(), args.lambda0, args.se, args.sim, args.thin, args.burn,
            args.entropy
        )
        rows.append(
            {
                "repeat": i + 1,
                "delta": value,
                "time_s": time.perf_counter() - t0,
                "first_delta": first_delta,
                "first_time_s": first_time,
            }
        )

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "repeat", "delta", "time_s", "first_delta", "first_time_s"
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
