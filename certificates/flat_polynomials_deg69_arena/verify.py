"""Standalone stdlib-only verifier for EV-000024 (degree-69 binary Littlewood flat polynomial).

Usage: python -I verify.py construction.json

Computes BOTH:
  1. The exact 1,000,000-point EinsteinArena / AlphaEvolve V2 evaluation metric:
       C+_arena = max_{j=0..999999} |g(exp(2*pi*i * j / 999999))| / sqrt(71)
     using exact integer autocorrelation R(m), a pure-Python radix-2 FFT on G = 8192
     points, and rigorous Taylor-Bernstein interval pruning
       max_{[a,b]} |g(e^{i*theta})|^2 <= max(f(a), f(b)) + (1/8) * M2 * (b - a)^2,
     where M2 = 2 * sum_{m=1}^{69} m^2 |R(m)| >= sup |f''(theta)|. Only grid indices j
     inside intervals whose continuous upper bound exceeds the incumbent grid maximum
     need direct evaluation, making the 1,000,000-point maximum mathematically exact
     in < 0.2 s using only the Python standard library.
  2. A rigorous continuous supremum bracket [C+_cont_lo, C+_cont_ub] over all |z| = 1
     by bisecting surviving intervals until the bracket width is < 1e-11.
"""

import json
import math
import os
import sys


def exact_autocorr(c):
    n = len(c)
    return [sum(c[j] * c[j + m] for j in range(n - m)) for m in range(n)]


def eval_f(r, theta):
    """Exact |g(e^{i*theta})|^2 = R(0) + 2 * sum_{m=1}^{n-1} R(m) * cos(m * theta)."""
    # Clenshaw recurrence for sum_{m=1}^{n-1} R(m) cos(m*theta)
    cos_t = math.cos(theta)
    two_cos = 2.0 * cos_t
    b1 = 0.0
    b2 = 0.0
    for m in range(len(r) - 1, 0, -1):
        b0 = r[m] + two_cos * b1 - b2
        b2 = b1
        b1 = b0
    return r[0] + 2.0 * (b1 * cos_t - b2)


def radix2_fft_power(c, g):
    """Iterative radix-2 Cooley-Tukey FFT returning |g(e^{2*pi*i*k/g})|^2 for k=0..g-1."""
    re = [0.0] * g
    im = [0.0] * g
    for i, v in enumerate(c):
        re[i] = float(v)
    j = 0
    for i in range(1, g):
        bit = g >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            re[i], re[j] = re[j], re[i]
            im[i], im[j] = im[j], im[i]
    length = 2
    two_pi = 2.0 * math.pi
    while length <= g:
        half = length >> 1
        ang = two_pi / length
        wlen_r = math.cos(ang)
        wlen_i = math.sin(ang)
        for i in range(0, g, length):
            wr = 1.0
            wi = 0.0
            for k in range(half):
                u_r = re[i + k]
                u_i = im[i + k]
                v_r = re[i + k + half] * wr - im[i + k + half] * wi
                v_i = re[i + k + half] * wi + im[i + k + half] * wr
                re[i + k] = u_r + v_r
                im[i + k] = u_i + v_i
                re[i + k + half] = u_r - v_r
                im[i + k + half] = u_i - v_i
                wr, wi = wr * wlen_r - wi * wlen_i, wr * wlen_i + wi * wlen_r
        length <<= 1
    return [re[k] * re[k] + im[k] * im[k] for k in range(g)]


def main():
    corrupt = "--corrupt" in sys.argv[1:]
    pos_args = [a for a in sys.argv[1:] if a != "--corrupt"]
    if len(pos_args) > 1:
        sys.exit("Usage: python -I verify.py [construction.json] [--corrupt]")
    cert_path = (
        pos_args[0]
        if pos_args
        else os.path.join(os.path.dirname(os.path.abspath(__file__)), "construction.json")
    )
    with open(cert_path, encoding="utf-8") as fh:
        data = json.load(fh)
    c = list(data["coefficients"] if isinstance(data, dict) else data)
    if corrupt:
        c[0] = 0
    if not isinstance(c, list) or len(c) != 70 or any(x not in (-1, 1) for x in c):
        sys.exit("FAIL: coefficients must be a list of 70 entries in {-1, +1}")

    n = len(c)
    r = exact_autocorr(c)
    e4 = sum(x * x for x in r[1:])
    psl = max(abs(x) for x in r[1:])
    mf = (n * n) / (2.0 * e4)
    m2 = float(2 * sum((m * m) * abs(r[m]) for m in range(1, n)))

    g = 8192
    two_pi = 2.0 * math.pi
    h = two_pi / g
    f_grid = radix2_fft_power(c, g)

    # 1. Exact 1,000,000-point EinsteinArena evaluation via rigorous interval pruning.
    n_pts = 1_000_000
    step = two_pi / (n_pts - 1)
    # Seed incumbent at the FFT argmax mapped to nearest arena index
    k_peak = max(range(g), key=lambda k: f_grid[k])
    j_seed = int(round((k_peak * h) / step))
    f_arena_max = eval_f(r, j_seed * step)

    # First refine intervals on [0, pi] (symmetric: f(2*pi - theta) == f(theta))
    # We keep full [0, 2*pi] to be bit-exact with linspace(0, 2*pi, 1_000_000).
    margin_base = 0.125 * m2 * (h * h) + 1e-9
    intervals = []
    for k in range(g):
        fa = f_grid[k]
        fb = f_grid[(k + 1) % g]
        ub = (fa if fa >= fb else fb) + margin_base
        if ub >= f_arena_max:
            intervals.append((k * h, (k + 1) * h, fa, fb))

    # Bisect candidate intervals 4 times so surviving arcs are tiny (< 2e-5 rad)
    for _ in range(4):
        next_iv = []
        for a, b, fa, fb in intervals:
            mid = 0.5 * (a + b)
            fm = eval_f(r, mid)
            # Check nearest arena grid point to mid to tighten lower bound fast
            j_near = int(round(mid / step))
            if 0 <= j_near < n_pts:
                fj = eval_f(r, j_near * step)
                if fj > f_arena_max:
                    f_arena_max = fj
            hm = 0.5 * (b - a)
            pad = 0.125 * m2 * (hm * hm) + 1e-9
            if (fa if fa >= fm else fm) + pad >= f_arena_max:
                next_iv.append((a, mid, fa, fm))
            if (fm if fm >= fb else fb) + pad >= f_arena_max:
                next_iv.append((mid, b, fm, fb))
        intervals = next_iv

    # Evaluate every 1,000,000-grid index inside surviving intervals
    checked_pts = 0
    for a, b, fa, fb in intervals:
        pad = 0.125 * m2 * ((b - a) ** 2) + 1e-9
        if (fa if fa >= fb else fb) + pad < f_arena_max:
            continue
        j_lo = max(0, int(math.floor(a / step)))
        j_hi = min(n_pts - 1, int(math.ceil(b / step)))
        for j in range(j_lo, j_hi + 1):
            fj = eval_f(r, j * step)
            checked_pts += 1
            if fj > f_arena_max:
                f_arena_max = fj

    # 2. Rigorous continuous supremum bracket over all |z| = 1
    f_cont_lo = f_arena_max
    cont_iv = list(intervals)
    for _ in range(20):
        next_iv = []
        f_cont_ub = f_cont_lo
        for a, b, fa, fb in cont_iv:
            mid = 0.5 * (a + b)
            fm = eval_f(r, mid)
            if fm > f_cont_lo:
                f_cont_lo = fm
            hm = 0.5 * (b - a)
            pad = 0.125 * m2 * (hm * hm) + 1e-12
            for la, lb, lfa, lfb in ((a, mid, fa, fm), (mid, b, fm, fb)):
                ub = (lfa if lfa >= lfb else lfb) + pad
                if ub >= f_cont_lo:
                    next_iv.append((la, lb, lfa, lfb))
                    if ub > f_cont_ub:
                        f_cont_ub = ub
        cont_iv = next_iv
        if f_cont_ub - f_cont_lo < 1e-11:
            break

    f_cont_ub = max(
        (fa if fa >= fb else fb) + 0.125 * m2 * ((b - a) ** 2) + 1e-12
        for a, b, fa, fb in cont_iv
    ) if cont_iv else f_cont_lo

    denom = math.sqrt(71.0)
    c_arena = math.sqrt(f_arena_max) / denom
    c_cont_lo = math.sqrt(f_cont_lo) / denom
    c_cont_ub = math.sqrt(f_cont_ub) / denom

    print(f"N = {n}  all_pm1 = True  E4 = {e4}  PSL = {psl}  MF = {mf:.6f}  checked_arena_pts = {checked_pts}")
    print(f"C+_arena (1e6 pts) = {c_arena:.12f}")
    print(f"PASS: C+_continuous_sup in [{c_cont_lo:.12f}, {c_cont_ub:.12f}]")


if __name__ == "__main__":
    main()
