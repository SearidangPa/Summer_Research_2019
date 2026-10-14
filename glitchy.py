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
    uv run glitchy.py recover --N 3001 --r 50 --k 3 --samples 100  # chance of recovering r
    uv run glitchy.py check                                   # every number in glitchy_period_finding_report.pdf
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
    """Legacy block-shift formula (exact Model I when r divides N).

    For nondivisible N this uses ceil(N/r), not the actual truncated literal
    sequence; the nondivisible proof checks use an exact literal-prefix FFT.
    """
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
# Period recovery on literal prefixes (the Section 5 proposal)
# ---------------------------------------------------------------------------
def literal_prefix(N, r, positions, signs):
    """Apply edits in order to the tagged sequence 0, 1, 2, ... and return N symbols.

    An insertion at `pos` puts a new symbol r + j before the current entry `pos`;
    a deletion removes the current entry `pos`.  Original tag t has symbol t mod r.
    """
    sequence = list(range(N + len(signs)))
    for j, (pos, sign) in enumerate(zip(positions, signs)):
        if sign > 0:
            sequence.insert(pos, -(r + j) - 1)
        else:
            sequence.pop(pos)
    tags = np.array(sequence[:N])
    return np.where(tags >= 0, tags % r, -tags - 1)


def literal_probability(entries):
    """Exact Fourier measurement distribution of a finite symbol sequence."""
    _, inverse = np.unique(entries, return_inverse=True)
    classes = inverse[None, :] == np.arange(inverse.max() + 1)[:, None]
    return (np.abs(np.fft.fft(classes, axis=1)) ** 2).sum(0) / len(entries) ** 2


def convergent_denominators(a, N, limit):
    """Denominators q <= limit of the continued-fraction convergents of a/N."""
    qs, (q_prev, q), (num, den) = [], (1, 0), (a, N)
    while den:
        t, num, den = num // den, den, num % den
        q_prev, q = q, t * q + q_prev
        if q > limit:
            break
        qs.append(q)
    return qs


def recovery(N, r, entries, radius, tau):
    """Exact success probability of one run of the Section 5 procedure.

    Measure y; collect convergent denominators q <= sqrt(N) of y'/N for |y' - y| <= radius;
    output the smallest q whose shift test f(x) = f(x+q) passes on more than a 1 - tau
    fraction of the prefix.  Returns (mass within radius of any peak, mass within radius
    of an informative peak gcd(c, r) = 1, probability that the output is r).
    """
    probability = literal_probability(entries)
    good = good_outcomes(N, r)
    offsets = np.arange(-radius, radius + 1)
    informative = [c for c in range(r) if math.gcd(c, r) == 1]
    near = np.unique((good[:, None] + offsets) % N)
    near_informative = np.unique((good[informative][:, None] + offsets) % N)

    limit = math.isqrt(N)
    candidates = [convergent_denominators(y, N, limit) for y in range(N)]
    passes = {q: np.mean(entries[:-q] == entries[q:]) > 1 - tau for q in range(1, limit + 1)}
    success = 0.0
    for y in range(N):
        qs = sorted({q for h in offsets for q in candidates[(y + h) % N]})
        if next((q for q in qs if passes[q]), None) == r:
            success += probability[y]
    return float(probability[near].sum()), float(probability[near_informative].sum()), float(success)


def recovery_report(N, r, placements, radius, tau):
    """Print recovery numbers for each placement [(positions, signs), ...] and a summary."""
    phi = sum(math.gcd(c, r) == 1 for c in range(r))
    _, _, ideal = recovery(N, r, literal_prefix(N, r, [], []), radius, tau)
    print(f"no edits: success {ideal:.4f};  phi(r)/r = {phi / r:.4f};  "
          f"conjectured floor (1/2) phi(r)/r = {phi / (2 * r):.4f}")
    rows = [(*recovery(N, r, literal_prefix(N, r, pos, signs), radius, tau), pos, signs)
            for pos, signs in placements]
    for mass, informative, success, pos, signs in rows[:20]:
        sign_text = "".join("+" if s > 0 else "-" for s in signs)
        print(f"  edits {sign_text} at {pos}: window mass {mass:.4f}, informative {informative:.4f}, "
              f"success {success:.4f}")
    if len(rows) > 1:
        worst = min(rows, key=lambda row: row[2])
        print(f"{len(rows)} placements: success min {worst[2]:.4f} (at {worst[3]}), "
              f"mean {np.mean([row[2] for row in rows]):.4f}; "
              f"window mass min {min(row[0] for row in rows):.4f}")
    return rows


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
def assert_literal_bounds(sequence, N, r, k):
    """Check one tagged prefix against direct Fourier measurement probabilities.

    Nonnegative tags are original indices. Negative tag -w-1 denotes an inserted
    value w, even when w equals an original periodic value.
    """
    encoded = np.asarray(sequence)
    assert len(encoded) == N and N % r == 0
    original = encoded >= 0
    inserted = int(np.count_nonzero(~original))
    entries = np.where(original, encoded % r, -encoded - 1)
    labels, inverse = np.unique(entries, return_inverse=True)
    classes = inverse[None, :] == np.arange(len(labels))[:, None]
    probability = (np.abs(np.fft.fft(classes, axis=1)) ** 2).sum(0) / N**2
    counts = np.zeros((len(labels), r), dtype=int)
    retained = np.zeros_like(counts)
    residues = np.arange(N) % r
    np.add.at(counts, (inverse, residues), 1)
    np.add.at(retained, (inverse[original], residues[original]), 1)
    displacements = np.arange(N)[original] - encoded[original]
    good = good_outcomes(N, r)
    near = np.unique((good[:, None] + np.arange(-k, k + 1)) % N)
    exact_mass = float(probability[good].sum())
    residue_mass = r * np.sum(counts**2) / N**2
    support_bound = r / ((k + 1) * r + inserted)
    retained_bound = (N - inserted) ** 2 / ((k + 1) * N**2)
    assert 0 <= inserted <= min(k, N)
    assert len(np.unique(displacements)) <= k + 1
    assert np.count_nonzero(retained) <= (k + 1) * r
    assert np.count_nonzero(counts) <= (k + 1) * r + inserted
    assert abs(probability.sum() - 1) < 1e-12
    assert abs(exact_mass - residue_mass) < 1e-12
    assert exact_mass >= max(support_bound, retained_bound) - 1e-12
    assert support_bound >= 1 / (2 * k + 1) - 1e-12
    assert probability[near].sum() >= exact_mass - 1e-12
    return exact_mass


def check_literal_glitches():
    """Exhaustive small edits, seeded general-k cases, and the worked example."""
    def one_edit(sequence, values):
        for pos in range(len(sequence) + 1):
            for value in values:
                yield sequence[:pos] + (value,) + sequence[pos:]
        for pos in range(len(sequence)):
            yield sequence[:pos] + sequence[pos + 1:]

    cases = [(N, r, 2) for N, r in
             [(1, 1), (2, 1), (2, 2), (6, 1), (6, 2), (6, 3), (18, 3), (24, 4)]]
    cases += [(1, 1, 3), (4, 2, 3), (6, 3, 3)]
    for N, r, k in cases:
        # Spare originals make this exactly the infinite-sequence prefix model,
        # even if every operation deletes an entry in the prefix.
        states = {tuple(range(N + k))}
        for step in range(k):
            states = {
                edited[:N] if step == k - 1 else edited
                for state in states
                for edited in one_edit(state, range(-r - 2, 0))
            }
        minimum = min(assert_literal_bounds(state, N, r, k) for state in states)
        print(f"N={N} r={r} k={k}: {len(states)} exhaustive tagged prefixes; "
              f"minimum exact mass {minimum:.6f}, uniform bound {1 / (2 * k + 1):.6f}")

    rng = np.random.default_rng(2019)
    checked = 0
    for N, r in [(1, 1), (2, 2), (6, 3), (12, 4), (30, 5), (72, 8), (60, 15)]:
        for k in [0, 1, 2, 3, 4, 8, 16, 40, 80]:
            for trial in range(40):
                sequence = list(range(N + k))
                A = D = 0
                for j in range(k):
                    insert = (trial % 4 == 0 or
                              (trial % 4 == 2 and j % 2 == 0) or
                              (trial % 4 == 3 and bool(rng.integers(2))))
                    pos = int(rng.integers(len(sequence) + int(insert)))
                    if insert:
                        # Old values, repeated new values, and distinct new values.
                        value = [int(rng.integers(r)), r, r + j][j % 3]
                        sequence.insert(pos, -value - 1)
                        A += 1
                    else:
                        sequence.pop(pos)
                        D += 1
                prefix = sequence[:N]
                for x, t in enumerate(prefix):
                    if t >= 0:
                        assert -D <= x - t <= A
                assert sum(t < 0 for t in prefix) <= A
                assert_literal_bounds(prefix, N, r, k)
                checked += 1
    # Entire prefix inserted; an insertion subsequently deleted; edits after N.
    for sequence, N, r, k in [
        ([-4, -4, -4], 3, 3, 3),
        (list(range(12)), 12, 3, 2),
        ([-1, -2, -3, -4], 4, 2, 8),
    ]:
        assert_literal_bounds(sequence, N, r, k)
    print(f"{checked} seeded general-k cases checked, including k > N and k >= r.")

    example = list(range(18))
    example.insert(3, -4)  # d before original 3
    example.insert(7, -5)  # e before original 6
    assert example.pop(11) == 9
    mass = assert_literal_bounds(example[:15], 15, 3, 3)
    assert abs(mass - 1 / 3) < 1e-12
    print("Worked example: N=15 r=3 k=3 I=2; exact mass 1/3 >= 3/14 >= 1/7.")

    for k in [0, 1, 2, 3, 5, 10]:
        r, M = k + 2, 2 * (k + 1)
        N, blocks = r * M, list(range(2, M, 2))
        for signs in [[1] * k, [-1] * k, [(-1) ** j for j in range(k)]]:
            lengths = np.diff([0, *blocks, M]) / M
            grouped = np.bincount(shifts(signs) % r, weights=lengths, minlength=r)
            exact_mass = distribution(N, r, blocks, signs, good_outcomes(N, r)).sum()
            assert abs(exact_mass - np.sum(grouped**2)) < 1e-12
            assert exact_mass >= 1 / (k + 1) - 1e-12
            assert abs(distribution(N, r, blocks, signs).sum() - 1) < 1e-12
            if len(set(signs)) <= 1:
                assert abs(exact_mass - 1 / (k + 1)) < 1e-12
    print("Model I: grouped-length identity and equal-block equality cases checked.")


def assert_neighbourhood_bound(sequence, N, r, k):
    """Check the division-free Fejer proof using an actual literal prefix.

    Tagged entries have the same encoding as assert_literal_bounds. These checks
    deliberately do not use distribution(), whose nondivisible block-shift
    formula is not the literal-prefix measurement distribution.
    """
    encoded = np.asarray(sequence)
    assert len(encoded) == N and 2 <= r <= N and k >= 1
    inserted = int(np.count_nonzero(encoded < 0))
    assert 0 <= inserted <= min(k, N)
    entries = np.where(encoded >= 0, encoded % r, -encoded - 1)
    labels, inverse = np.unique(entries, return_inverse=True)
    classes = inverse[None, :] == np.arange(len(labels))[:, None]
    probability = (np.abs(np.fft.fft(classes, axis=1)) ** 2).sum(0) / N**2
    good = good_outcomes(N, r)
    near = np.unique((good[:, None] + np.arange(-k, k + 1)) % N)
    mass = float(probability[near].sum())
    assert abs(probability.sum() - 1) < 1e-12
    assert mass >= 133 / 358 - 1e-12
    assert mass >= 1 / (2 * k + 1) - 1e-12

    t = 2 * N / (r * (2 * k + 1))
    if t <= 2:
        assert len(near) == N
        assert abs(mass - 1) < 1e-12
        return mass

    T = 2 if t <= 6 else math.ceil(3 * t / 4)
    ell = np.arange(1, T)
    matching = np.array([np.mean(entries == np.roll(entries, -int(h)))
                         for h in ell * r])
    assert np.all(matching >= 1 - ((k + 1) * ell * r + inserted) / N - 1e-12)
    assert np.all(matching >= 1 - (4 * ell + 1) / (3 * t) - 1e-12)
    phase = r * np.arange(N) / N
    moments = np.cos(2 * np.pi * np.outer(ell, phase)) @ probability
    assert np.allclose(moments, matching, atol=1e-11, rtol=0)

    # Build the kernel as a squared geometric sum, independently of its moments.
    kernel = np.abs(np.exp(2j * np.pi * np.outer(np.arange(T), phase)).sum(0)) ** 2 / T
    expectation = float(probability @ kernel)
    assert abs(expectation - (1 + 2 * np.dot(1 - ell / T, matching))) < 1e-10
    lower = 1 - 4 * (T**2 - 1) / (9 * t * T) - (T - 1) / (3 * t * T)
    assert expectation / T >= lower - 1e-11
    outside = np.ones(N, dtype=bool)
    outside[near] = False
    distance = np.abs((phase + 0.5) % 1 - 0.5)
    assert np.all(distance[outside] >= 1 / t - 1e-12)

    if t <= 6:
        q = math.cos(np.pi / t) ** 2
        assert t * math.sin(np.pi / t) ** 2 >= 4 / 3 - 1e-12
        assert lower >= 1 - 5 / (6 * t) - 1e-12
        assert (lower - q) / (1 - q) >= 3 / 8 - 1e-12
    else:
        q = 1 / (T**2 * math.sin(np.pi / t) ** 2)
        assert lower >= 29 / 54 - 1e-12
        assert q <= 64 / 243 + 1e-12
        assert (lower - q) / (1 - q) >= 133 / 358 - 1e-12
    assert np.all(kernel[outside] / T <= q + 1e-11)
    assert expectation / T <= mass + (1 - mass) * q + 1e-11
    assert mass >= (lower - q) / (1 - q) - 1e-11
    return mass


def check_nondividing_glitches():
    """Exhaustive and seeded literal edits: both Fejer regimes and full coverage."""
    def one_edit(sequence, r):
        for pos in range(len(sequence) + 1):
            for value in range(r + 1):  # Every old value and one new value.
                yield sequence[:pos] + (-value - 1,) + sequence[pos:]
        for pos in range(len(sequence)):
            yield sequence[:pos] + sequence[pos + 1:]

    for N, r, k in [(7, 2, 1), (9, 2, 1), (11, 3, 1), (31, 2, 1),
                    (9, 2, 2), (13, 2, 2), (17, 3, 2)]:
        states = {tuple(range(N + k))}
        for step in range(k):
            states = {
                edited[:N] if step == k - 1 else edited
                for state in states
                for edited in one_edit(state, r)
            }
        minimum = min(assert_neighbourhood_bound(state, N, r, k) for state in states)
        print(f"N={N} r={r} k={k}: {len(states)} exhaustive prefixes; "
              f"minimum neighbourhood mass {minimum:.6f} >= {133 / 358:.6f}")

    rng = np.random.default_rng(2019)
    checked = 0
    for N, r in [(7, 2), (12, 2), (31, 2), (37, 3), (101, 7),
                 (211, 17), (257, 107), (1009, 17), (2003, 31)]:
        for k in [1, 2, 3, 4, 8, 16, 40, 80]:
            for trial in range(20):
                sequence = list(range(N + k))
                for j in range(k):
                    insert = (trial % 4 == 0 or
                              (trial % 4 == 2 and j % 2 == 0) or
                              (trial % 4 == 3 and bool(rng.integers(2))))
                    pos = int(rng.integers(len(sequence) + int(insert)))
                    if insert:
                        value = [int(rng.integers(r)), r, r + j][(j + trial) % 3]
                        sequence.insert(pos, -value - 1)
                    else:
                        sequence.pop(pos)
                assert_neighbourhood_bound(sequence[:N], N, r, k)
                checked += 1
    print(f"{checked} seeded cases checked: nondivisible and divisible controls, "
          "mixed signs, half-integer peaks, old/new symbols, and overlapping windows.")

    # The poster's nondivisible parameters, now with literal edits and an exact FFT.
    N, r, k = 11000, 107, 3
    sequence = list(range(N + k))
    sequence.insert(3 * r + 7, -r - 1)
    sequence.insert(6 * r + 2, -r - 2)
    sequence.pop(9 * r + 1)
    mass = assert_neighbourhood_bound(sequence[:N], N, r, k)
    print(f"Literal poster-sized example: N={N} r={r} k={k}; mass {mass:.6f} >= 1/3.")

    # No-edit leakage rules out blindly extending the divisible k=0 lower bound 1.
    entries = np.arange(5) % 2
    probability = (np.abs(np.fft.fft(entries[None, :] == np.arange(2)[:, None],
                                    axis=1)) ** 2).sum(0) / 25
    mass = float(probability[good_outcomes(5, 2)].sum())
    assert mass < 1 - 1e-12
    print(f"No edits, N=5 r=2: mass at rounded peaks {mass:.6f} < 1 (expected leakage).")


def check_recovery():
    """Exact success probabilities of the Section 5 procedure.

    The placements are the worst found by a local search over positions and signs
    (for k = 1, every third position with both signs). They are not proved worst cases.
    """
    from fractions import Fraction

    for N in [7, 12, 100, 360]:
        for a in range(N):
            x, terms = Fraction(a, N), []
            while True:
                terms.append(math.floor(x))
                if x == terms[-1]:
                    break
                x = 1 / (x - terms[-1])
            denominators = []
            for i in range(1, len(terms) + 1):
                f = Fraction(terms[i - 1])
                for t in reversed(terms[:i - 1]):
                    f = t + 1 / f
                denominators.append(f.denominator)
            expected = list(itertools.takewhile(lambda q: q <= math.isqrt(N), denominators))
            assert convergent_denominators(a, N, math.isqrt(N)) == expected
    print("Continued-fraction denominators agree with exact Fraction expansions.")

    # Without edits and with r | N, exactly the outcomes cN/r with gcd(c, r) = 1 succeed.
    for N, r in [(360, 12), (3000, 50)]:
        phi = sum(math.gcd(c, r) == 1 for c in range(r))
        assert abs(recovery(N, r, literal_prefix(N, r, [], []), 0, 0.25)[2] - phi / r) < 1e-12
    print("No edits, r | N, radius 0: success = phi(r)/r exactly.")

    for N, r, positions, signs in [
        (3000, 50, [2448], [1]),
        (3000, 50, [250, 2750], [-1, -1]),
        (3000, 50, [299, 600, 900], [1, -1, 1]),
        (3000, 50, [107, 308, 558, 2808], [1, 1, -1, 1]),
        (3001, 50, [501], [-1]),
        (3001, 50, [350, 2650], [1, 1]),
        (3001, 50, [600, 300, 850], [1, -1, -1]),
        (3001, 50, [203, 452, 700, 2752], [-1, 1, -1, 1]),
        (2500, 47, [471], [-1]),
        (2500, 47, [328, 610], [-1, 1]),
        (2500, 47, [282, 518, 753], [1, -1, 1]),
        (2500, 47, [1677, 1866, 2288, 2100], [1, -1, -1, 1]),
    ]:
        k = len(signs)
        phi = sum(math.gcd(c, r) == 1 for c in range(r))
        entries = literal_prefix(N, r, positions, signs)
        ideal = recovery(N, r, literal_prefix(N, r, [], []), k, 0.25)[2]
        mass, _, success = recovery(N, r, entries, k, 0.25)
        assert mass >= 133 / 358 - 1e-12
        assert success >= phi / (2 * r)
        for tau in [0.05, 0.5]:
            assert abs(recovery(N, r, entries, k, tau)[2] - success) < 1e-12
        sign_text = "".join("+" if s > 0 else "-" for s in signs)
        print(f"N={N} r={r} k={k} {sign_text} at {positions}: success {success:.4f}, "
              f"no edits {ideal:.4f}, ratio {success / ideal:.3f}  (phi(r)/r = {phi / r:.4f})")


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

    section("General literal glitches: residue identity and total mass at exact ideal outcomes")
    check_literal_glitches()

    section("Nondividing periods: matching moments, Fejer kernel and constant neighbourhood mass")
    check_nondividing_glitches()

    section("Section 5: exact success probability of the recovery procedure (numerical evidence only)")
    check_recovery()


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

    q = sub.add_parser("recover", help="exact success probability of the Section 5 recovery procedure")
    common(q)
    q.add_argument("--positions", type=int, nargs="+",
                   help="literal edit positions, applied in order (default: random)")
    q.add_argument("--samples", type=int, default=1, help="random placements to evaluate")
    q.add_argument("--tau", type=float, default=0.25, help="shift test passes above 1 - tau")

    sub.add_parser("check", help="reproduce every number in glitchy_period_finding_report.pdf")

    a = p.parse_args()
    if a.cmd == "check":
        return check()

    M = math.ceil(a.N / a.r)
    if a.cmd == "plot" and a.blocks:
        a.k = len(a.blocks)
    if a.cmd == "recover" and a.positions:
        a.k = len(a.positions)
    signs = a.signs * a.k if len(a.signs) == 1 else a.signs
    if len(signs) != a.k:
        p.error(f"--signs needs 1 or {a.k} values")
    radius = a.k if a.radius is None else a.radius

    if a.cmd == "plot":
        blocks = a.blocks or sorted(np.random.default_rng(a.seed).choice(np.arange(1, M), a.k, replace=False).tolist())
        plot(a.N, a.r, blocks, signs, radius, a.out)
    elif a.cmd == "recover":
        if a.r * a.r > a.N:
            p.error("continued fractions need r^2 <= N")
        rng = np.random.default_rng(a.seed)
        placements = ([(a.positions, signs)] if a.positions else
                      [(rng.integers(a.N, size=a.k).tolist(), signs) for _ in range(a.samples)])
        recovery_report(a.N, a.r, placements, radius, a.tau)
    else:
        res = worst_case_search(a.N, a.r, a.k, signs, radius, a.samples, a.seed)
        print(f"{len(res)} placements of {a.k} glitches in {M} copies; worst ratio within {radius}:")
        for worst, blocks in res[: a.top]:
            print(f"  {worst:.4f}  at copies {blocks}")


if __name__ == "__main__":
    main()
