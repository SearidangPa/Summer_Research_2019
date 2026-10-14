#!/usr/bin/env python3
"""Build the exact, data-derived counting illustration used in Section 4.1.

Standard library only. Run from any directory; the SVG is written beside this file.
"""
from collections import Counter
from fractions import Fraction
from html import escape
from pathlib import Path

N, r, k = 15, 3, 3
sequence = [(t, "abc"[t % r]) for t in range(N + k)]
sequence.insert(3, (-1, "d"))  # the j-th insertion gets tag -j
sequence.insert(7, (-2, "e"))
assert sequence.pop(11) == (9, "a")
sequence = sequence[:N]
counts = Counter((symbol, x % r) for x, (_, symbol) in enumerate(sequence))
squares = sum(n * n for n in counts.values())
assert squares == 25 and Fraction(r * squares, N * N) == Fraction(1, 3)
assert sum(t < 0 for t, _ in sequence) == 2

parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="690" viewBox="0 0 1000 690">',
         '<title>Three literal glitches: from original indices to residue counts</title>',
         '<desc>Two insertions and one deletion in a period-three sequence. '
         'The count table has total 15 and sum of squares 25. '
         'The probability on ideal outcomes 0, 5, and 10 is one third.</desc>',
         '<rect width="1000" height="690" fill="white"/>']
ink, blue, green, muted = "#24384b", "#1f4e79", "#276442", "#5d6c79"


def rect(x, y, w, h, fill, stroke="none", radius=0):
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" '
                 f'fill="{fill}" stroke="{stroke}"/>')


def text(x, y, content, size=18, color=ink, bold=False, anchor="start"):
    parts.append(f'<text x="{x}" y="{y}" font-family="Arial, sans-serif" '
                 f'font-size="{size}" fill="{color}" text-anchor="{anchor}"'
                 f' font-weight="{"bold" if bold else "normal"}">{escape(str(content))}</text>')


text(24, 30, "1. Follow the original entries", 23, blue, True)
text(24, 57, "r = 3, N = 15, k = 3.  Insert d before original 3; insert e before original 6; delete original 9.", 17)
labels = [("Position x", 99), ("Symbol", 138), ("Tag τ(x)", 177), ("s = x − τ(x)", 216)]
for label, y in labels:
    text(24, y, label, 17)
palette = {0: "#e6eff8", 1: "#e4f2e9", 2: "#fff1d6", None: "#f0e7f8"}
for x, (t, symbol) in enumerate(sequence):
    left, width = 194 + x * 51, 46
    shift = None if t < 0 else x - t
    rect(left, 76, width, 151, palette[shift], radius=4)
    for entry, (_, y) in zip([x, symbol, f"−{-t}" if t < 0 else t,
                              "—" if t < 0 else shift], labels):
        text(left + width / 2, y, entry, 19, anchor="middle", bold=(y == 138))
text(194, 253, "The same displacement may recur: the deletion brings the last block from 2 back to 1.", 17, muted)
rect(24, 273, 952, 1, "#dbe3eb")

text(24, 308, "2. Put entries into (symbol, residue) boxes", 22, blue, True)
text(24, 335, "Each column collects positions with the same x mod 3.", 17, muted)
left, top, cw, ch = 123, 379, 75, 44
text(36, 367, "symbol u", 16, muted)
for b in range(r):
    text(left + b * cw + cw / 2, 367, f"b = {b}", 17, anchor="middle", bold=True)
for row, symbol in enumerate("abcde"):
    text(65, top + row * ch + 29, symbol, 20, bold=True, anchor="middle")
    for b in range(r):
        n = counts[symbol, b]
        fill = ["#fafbfd", "#e5eef7", "#c4d9ed", "#94b9dc"][n]
        if symbol in "de" and n:
            fill = "#eadbf5"
        rect(left + b * cw, top + row * ch, cw - 3, ch - 3, fill, "#d4dde6", 3)
        text(left + b * cw + (cw - 3) / 2, top + row * ch + 28, n, 20,
             color=muted if not n else ink, anchor="middle")
text(24, 628, "Total entries: 15.   Sum of squared counts: 25.", 17, bold=True)
text(24, 653, "d and e are inserted symbols; each adds one box.", 17, muted)

rect(397, 353, 579, 305, "#f3f7fb", radius=8)
text(420, 388, "3. Read off the exact Fourier mass", 22, blue, True)
text(420, 426, "Ideal outcomes: G = {0, 5, 10}", 20)
text(420, 468, "P(G) = (r / N²) × sum of squared counts", 20)
text(420, 505, "         = (3 / 225) × 25 = 1/3", 22, green, True)
text(420, 550, "The theorem guarantees at least", 18)
text(420, 584, "r / ((k + 1)r + I) = 3/14", 22, blue, True)
text(420, 628, "This also exceeds the universal bound 1/7.", 18, muted)
parts.append("</svg>")
out = Path(__file__).with_name("general_k_example.svg")
out.write_text("\n".join(parts) + "\n")
print(f"{out}: exact mass = 1/3, support bound = 3/14, universal bound = 1/7")
