#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Glitchy quantum period finding: simulator for the block-shift model ("Model I").

A periodic sequence with period r and length N consists of M = ceil(N/r) copies of the
period.  Glitch j sits on the boundary before copy m_j (0 < m_1 <= ... <= m_k < M) and
moves every later copy by +1 (insertion) or -1 (deletion).  The quantum period-finding
measurement distribution is then

    P_k(y) = |z(y)|^2 / (M N),   z(y) = sum_j e^{2 pi i y s_j / N} sum_{m in block j} e^{2 pi i m r y / N},

where s_j is the cumulative shift of block j.  This is exact when r divides N (see
glitchy_period_finding_report.pdf, Definition 1).  Probabilities are reported relative to the ideal
probability 1/r of each good outcome, so "1" means "as good as with no glitch".

Usage (uv installs numpy and matplotlib automatically):

    uv run glitchy.py plot  --N 3000 --r 20 --k 2              # random placement
    uv run glitchy.py plot  --N 3000 --r 50 --blocks 6 54      # given placement
    uv run glitchy.py worst --N 3000 --r 50 --k 2              # exhaustive worst case
    uv run glitchy.py worst --N 1500 --r 50 --k 4 --samples 2000
    uv run glitchy.py check                                    # every number in glitchy_period_finding_report.pdf
"""
import argparse
import itertools
import math

import numpy as np


# ---------------------------------------------------------------------------
# Model I
# ---------------------------------------------------------------------------
def shifts(signs):
    """Cumulative shifts s_1 = 0, s_{j+1} = s_j + sign_j, one per block."""
    return np.concatenate([[0], np.cumsum(signs)]).astype(int)


def distribution(N, r, blocks, signs, y=None):
    """P_k(y) for the outcomes y (default: all of 0..N-1)."""
    M = math.ceil(N / r)
    b = [0, *blocks, M]
    s = shifts(signs)
    y = np.arange(N) if y is None else np.asarray(y)
    z = np.zeros(y.shape, complex)
    for j in range(len(b) - 1):
        m = np.arange(b[j], b[j + 1])
        z += np.exp(2j * np.pi * np.multiply.outer(y, m) * r / N).sum(-1) * np.exp(2j * np.pi * y * s[j] / N)
    return np.abs(z) ** 2 / (M * N)


def good_outcomes(N, r):
    """The integer nearest c N / r for c = 0..r-1 (equal to c M when r divides N)."""
    return np.rint(np.arange(r) * N / r).astype(int) % N


def ratios(N, r, blocks, signs, radius):
    """R[c, h + radius] = P_k(good_c + h) / (1/r) for |h| <= radius."""
    y = (good_outcomes(N, r)[:, None] + np.arange(-radius, radius + 1)) % N
    return r * distribution(N, r, blocks, signs, y)


def worst_ratio(N, r, blocks, signs, radius=None):
    """min over c of the probability within `radius` of good outcome c, relative to ideal."""
    radius = len(blocks) if radius is None else radius
    return ratios(N, r, blocks, signs, radius).sum(1).min()


def worst_case_search(N, r, k, signs, radius=None, samples=None, seed=None):
    """Worst ratio for every placement of k glitches (or `samples` random ones).

    Returns a list of (worst ratio, blocks), sorted from worst to best.
    """
    radius = k if radius is None else radius
    M = math.ceil(N / r)
    y = ((good_outcomes(N, r)[:, None] + np.arange(-radius, radius + 1)) % N).ravel()
    # prefix sums over copies: block sum over [a, b) = C[:, b] - C[:, a]
    C = np.zeros((len(y), M + 1), complex)
    C[:, 1:] = np.cumsum(np.exp(2j * np.pi * np.outer(y, np.arange(M)) * r / N), axis=1)
    phase = np.exp(2j * np.pi * np.outer(y, shifts(signs)) / N)

    if samples is None:
        placements = itertools.combinations(range(1, M), k)
    else:
        rng = np.random.default_rng(seed)
        unique = set()
        while len(unique) < min(samples, math.comb(M - 1, k)):
            unique.add(tuple(sorted(rng.choice(np.arange(1, M), k, replace=False).tolist())))
        placements = iter(sorted(unique))

    results = []
    while chunk := list(itertools.islice(placements, 512)):
        B = np.array([[0, *p, M] for p in chunk])
        z = np.einsum("ypj,yj->py", C[:, B[:, 1:]] - C[:, B[:, :-1]], phase)
        worst = (r * np.abs(z) ** 2 / (M * N)).reshape(len(chunk), r, -1).sum(2).min(1)
        results += zip(worst.tolist(), [list(map(int, p)) for p in chunk])
    return sorted(results)


# ---------------------------------------------------------------------------
# Continuum limit (M -> infinity): glitch j at u_j = m_j / M in (0, 1)
# ---------------------------------------------------------------------------
def continuum_coeffs(u, signs, alpha, radius):
    """a_h = int_0^1 phi(t) e^{2 pi i h t} dt for |h| <= radius, phi = e^{2 pi i alpha s_j} on piece j."""
    u = np.concatenate([[0.0], u, [1.0]])
    phi = np.exp(2j * np.pi * alpha * shifts(signs))
    a = []
    for h in range(-radius, radius + 1):
        if h == 0:
            a.append(np.sum(phi * np.diff(u)))
        else:
            e = np.exp(2j * np.pi * h * u)
            a.append(np.sum(phi * np.diff(e)) / (2j * np.pi * h))
    return np.array(a)


def continuum_ratio(u, signs, alpha, radius):
    """Probability within `radius` of the good outcome alpha = c/r, relative to ideal."""
    return float(np.sum(np.abs(continuum_coeffs(u, signs, alpha, radius)) ** 2))


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def plot(N, r, blocks, signs, radius, out):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.ticker import MaxNLocator

    blue, orange, muted, ink = "#2a78d6", "#eb6834", "#898781", "#52514e"
    ramp = LinearSegmentedColormap.from_list("blue", ["#fcfcfb", "#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])

    R = ratios(N, r, blocks, signs, radius)
    c = np.arange(r)
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 7), height_ratios=[1.3, 1],
                                      layout="constrained")
    for ax in (top, bottom):
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color(muted)
        ax.tick_params(colors=ink, labelsize=9)

    top.axhline(1, color=muted, lw=1)
    top.axhline(0.5, color=muted, lw=1, ls="--")
    top.text(r - 0.5, 1.01, "no glitch", color=ink, fontsize=8, ha="right", va="bottom")
    top.text(r - 0.5, 0.51, "conjectured bound 1/2", color=ink, fontsize=8, ha="right", va="bottom")
    top.plot(c, R.sum(1), "o-", color=blue, ms=4, lw=1.5, label=f"within {radius} of the good outcome")
    top.plot(c, R[:, radius], "o", color=orange, ms=4, label="exactly at the good outcome")
    top.set_ylim(0, 1.1)
    top.set_ylabel("probability / ideal", color=ink)
    top.grid(axis="y", color="#e1e0d9", lw=0.6)
    top.legend(frameon=False, fontsize=9, ncol=2, loc="lower left", bbox_to_anchor=(0, 1))
    sign_text = "".join("+" if s > 0 else "-" for s in signs) or "none"
    fig.suptitle(f"N = {N}, r = {r}, {len(blocks)} glitches at copies {list(blocks)}, signs {sign_text}",
                 fontsize=10, color=ink, x=0.01, ha="left")

    im = bottom.imshow(R.T, aspect="auto", origin="lower", cmap=ramp, vmin=0, interpolation="nearest",
                       extent=[-0.5, r - 0.5, -radius - 0.5, radius + 0.5])
    bottom.set_ylabel("offset h from good outcome", color=ink)
    bottom.yaxis.set_major_locator(MaxNLocator(integer=True))
    bottom.set_xlabel("good outcome index c   (outcome nearest cN/r)", color=ink)
    cb = fig.colorbar(im, ax=bottom, orientation="horizontal", location="bottom", aspect=50, shrink=0.5)
    cb.set_label("probability / ideal", color=ink, fontsize=9)
    cb.outline.set_visible(False)
    fig.align_ylabels()
    fig.savefig(out, dpi=150)
    print(f"worst ratio within {radius}: {R.sum(1).min():.4f} (at c = {int(R.sum(1).argmin())}); saved {out}")


# ---------------------------------------------------------------------------
# Reproduce every number quoted in glitchy_period_finding_report.pdf
# ---------------------------------------------------------------------------
def check_literal_two_glitches():
    """Check the residue identity and total-mass bounds on small literal edit cases."""
    def one_edit(sequence, values):
        for pos in range(len(sequence) + 1):
            for value in values:
                yield sequence[:pos] + (value,) + sequence[pos:]
        for pos in range(len(sequence)):
            yield sequence[:pos] + sequence[pos + 1:]

    for N, r in [(1, 1), (2, 1), (2, 2), (6, 1), (6, 2), (6, 3), (18, 3), (24, 4)]:
        # Two spare originals allow two deletions before restricting to length N.
        # Encode inserted label w as -w-1 to track surviving insertions, even
        # when they match an original value. Include two possible new labels.
        original = tuple(x % r for x in range(N + 2))
        values = range(-r - 2, 0)
        sequences = {
            second[:N]
            for first in one_edit(original, values)
            for second in one_edit(first, values)
        }
        good = good_outcomes(N, r)
        near = np.unique((good[:, None] + np.arange(-2, 3)) % N)
        lower = r / (3 * r + 2)
        minimum = 1.0
        for sequence in sequences:
            encoded = np.asarray(sequence)
            inserted = int(np.count_nonzero(encoded < 0))
            entries = np.where(encoded < 0, -encoded - 1, encoded)
            classes = entries[None, :] == np.arange(r + 2)[:, None]
            probability = (np.abs(np.fft.fft(classes, axis=1)) ** 2).sum(0) / N**2
            counts = np.zeros((r + 2, r), dtype=int)
            np.add.at(counts, (entries, np.arange(N) % r), 1)
            exact_mass = float(probability[good].sum())
            residue_mass = r * np.sum(counts**2) / N**2
            assert 0 <= inserted <= 2
            assert np.count_nonzero(counts) <= 3 * r + inserted
            assert abs(probability.sum() - 1) < 1e-12
            assert abs(exact_mass - residue_mass) < 1e-12
            assert exact_mass >= r / (3 * r + inserted) - 1e-12
            assert exact_mass >= (N - inserted) ** 2 / (3 * N**2) - 1e-12
            assert probability[near].sum() >= exact_mass - 1e-12
            minimum = min(minimum, exact_mass)
        print(f"N={N} r={r}: {len(sequences)} literal sequences checked; "
              f"minimum exact-ideal mass {minimum:.6f}, bound {lower:.6f}")


def check():
    def section(title):
        print(f"\n== {title}")

    section("Table 2: 2019 logged worst cases, recomputed, then re-searched over all placements")
    runs = [  # (log, N, r, k, sign, logged worst blocks, logged worst probability)
        ("Plot2", 3000, 50, 2, +1, [6, 54], 0.015445582305137994),
        ("Plot3", 3000, 100, 3, -1, [3, 6, 27], 0.006135891098008545),
        ("Plot4", 1500, 50, 4, +1, [3, 24, 25, 28], 0.014458930695021374),
    ]
    for name, N, r, k, sign, blocks, logged in runs:
        signs = [sign] * k
        M = N // r
        print(f"{name}: N={N} r={r} k={k} sum P={distribution(N, r, blocks, signs).sum():.6f} "
              f"logged worst {logged * r:.4f} at {blocks}, recomputed {worst_ratio(N, r, blocks, signs):.4f}")
        res = worst_case_search(N, r, k, signs)
        worst, blocks = res[0]
        c = int(ratios(N, r, blocks, signs, k).sum(1).argmin())
        cont = min(continuum_ratio(np.array(blocks) / M, signs, cc / r, k) for cc in range(r))
        print(f"    all {len(res)} placements: worst {worst:.4f} at {blocks}, attained at c/r = {c / r:.2f}, "
              f"continuum {cont:.4f}")

    section("Sign patterns (N=3000, r=50): worst ratio over all placements")
    for signs in ([1, 1], [1, -1], [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1]):
        worst, blocks = worst_case_search(3000, 50, len(signs), signs)[0]
        c = int(ratios(3000, 50, blocks, signs, len(signs)).sum(1).argmin())
        print(f"signs {signs}: {worst:.4f} at {blocks}, attained at c/r = {c / 50:.2f}")

    section("One-glitch theorem")
    D = lambda u: 2 * u * (1 - u) - 4 * np.sin(np.pi * u) ** 2 / np.pi**2
    for u, alpha in [(0.1, 0.5), (0.2, 0.3), (0.37, 0.5)]:
        closed = 1 - (1 - math.cos(2 * math.pi * alpha)) * D(u)
        print(f"u={u} alpha={alpha}: closed form {closed:.6f}, direct {continuum_ratio([u], [1], alpha, 1):.6f}")
    uu = np.linspace(0, 1, 2_000_001)
    g = 1 - 2 * D(uu)
    half = uu[uu <= 0.5]
    print(f"min D = {D(uu).min():.1e}; min g = {g.min():.4f} at u = {uu[g.argmin()]:.4f}")
    print(f"Jordan bound min = {(1 - 4 * half + (4 + 32 / np.pi**2) * half**2).min():.6f}; "
          f"8/(8+pi^2) = {8 / (8 + np.pi**2):.6f}")
    for N, r in [(3000, 50), (20000, 20)]:
        print(f"Model I, N={N} r={r}: worst over all placements = {worst_case_search(N, r, 1, [1])[0][0]:.4f}")

    section("Section 8.1: packed family (k equal-sign glitches spaced evenly in [0, delta), alpha = 1/2)")
    for k in [2, 4, 6, 10, 20, 50, 100, 200]:
        vals = [(continuum_ratio(np.linspace(0, d, k + 1)[1:], [1] * k, 0.5, k), d)
                for d in np.linspace(0.05, 0.95, 181)]
        best = min(vals)
        print(f"k={k:4d}: min S = {best[0]:.4f} at delta = {best[1]:.3f}")

    section("Radius lower bound: equally spaced glitches, k odd, alpha = 1/2: S(L) for L = 0..k")
    for k in [3, 5, 9]:
        u = np.arange(1, k + 1) / (k + 1)
        print(f"k={k}: " + " ".join(f"{continuum_ratio(u, [1] * k, 0.5, L):.3f}" for L in range(k + 1)))

    section("Two literal glitches: residue identity and total mass at exact ideal outcomes")
    check_literal_two_glitches()


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(q):
        q.add_argument("--N", type=int, required=True, help="domain size")
        q.add_argument("--r", type=int, required=True, help="period")
        q.add_argument("--k", type=int, default=2, help="number of glitches")
        q.add_argument("--signs", type=int, nargs="+", default=[1],
                       help="+1 insertion / -1 deletion; one value applies to every glitch")
        q.add_argument("--radius", type=int, help="neighbourhood radius (default k)")
        q.add_argument("--seed", type=int, help="random seed")

    q = sub.add_parser("plot", help="plot the distribution around every good outcome")
    common(q)
    q.add_argument("--blocks", type=int, nargs="+", help="glitch copies m_1..m_k (default: random)")
    q.add_argument("--out", default="glitchy_plot.png")

    q = sub.add_parser("worst", help="search placements for the worst neighbourhood ratio")
    common(q)
    q.add_argument("--samples", type=int, help="random placements instead of all of them")
    q.add_argument("--top", type=int, default=10)

    sub.add_parser("check", help="reproduce every number in glitchy_period_finding_report.pdf")

    a = p.parse_args()
    if a.cmd == "check":
        return check()

    M = math.ceil(a.N / a.r)
    if a.cmd == "plot" and a.blocks:
        a.k = len(a.blocks)
    signs = a.signs * a.k if len(a.signs) == 1 else a.signs
    if len(signs) != a.k:
        p.error(f"--signs needs 1 or {a.k} values")
    radius = a.k if a.radius is None else a.radius

    if a.cmd == "plot":
        blocks = a.blocks or sorted(np.random.default_rng(a.seed).choice(np.arange(1, M), a.k, replace=False).tolist())
        plot(a.N, a.r, blocks, signs, radius, a.out)
    else:
        res = worst_case_search(a.N, a.r, a.k, signs, radius, a.samples, a.seed)
        print(f"{len(res)} placements of {a.k} glitches in {M} copies; worst ratio within {radius}:")
        for worst, blocks in res[: a.top]:
            print(f"  {worst:.4f}  at copies {blocks}")


if __name__ == "__main__":
    main()
