// Quantum Glitchy Period Finding: status report and open problems.
// Compile with:  typst compile glitchy_period_finding_report.typ
// Every number quoted below is reproduced by `uv run glitchy.py check` (repository root).

#set document(
  title: [Quantum Glitchy Period Finding: Status Report and Open Problems],
  author: ("Searidang Pa",),
)
#set page(paper: "us-letter", margin: (x: 1in, y: 1in), numbering: "1")
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show heading.where(level: 1): set block(above: 1.6em, below: 0.9em)
#show figure.caption: set text(size: 9pt)
#show link: set text(fill: rgb("#1f4e79"))
#show ref: it => {
  let el = it.element
  if el != none and el.func() == math.equation {
    link(el.location(), [Eq.~#numbering(el.numbering, ..counter(math.equation).at(el.location()))])
  } else { it }
}

// ---------------------------------------------------------------------------
// Status pills
// ---------------------------------------------------------------------------
#let pill(body, col) = box(
  inset: (x: 4pt, y: 1.5pt),
  outset: (y: 1pt),
  radius: 5pt,
  fill: col.lighten(88%),
  stroke: 0.4pt + col,
  text(size: 7.5pt, weight: "semibold", fill: col.darken(25%), body),
)
#let from2019 = pill([from 2019 work], rgb("#7a5200"))
#let claim2019 = pill([2019 claim, proof not on file], rgb("#7a5200"))
#let new = pill([new in this note], rgb("#1f4e79"))
#let proved = pill([proof below], rgb("#1f6f43"))
#let leanchecked = pill([Lean checked], rgb("#31647a"))
#let sketch = pill([proof sketch only], rgb("#b35900"))
#let numerical = pill([numerical evidence only], rgb("#5b3f99"))
#let conjecture = pill([conjecture], rgb("#a11d1d"))

// ---------------------------------------------------------------------------
// Theorem-like environments (share one counter)
// ---------------------------------------------------------------------------
#show figure.where(kind: "thm"): set block(breakable: true)
#show figure.where(kind: "thm"): it => {
  set align(left)
  set par(justify: false)
  block(
    width: 100%,
    inset: (x: 10pt, y: 8pt),
    radius: 2pt,
    fill: luma(248),
    stroke: (left: 2pt + rgb("#1f4e79")),
  )[#strong[#it.supplement #context it.counter.display(it.numbering)]#it.body]
}
#let env(supplement) = (name: none, tags: (), body) => figure(
  kind: "thm",
  supplement: supplement,
  numbering: "1",
  outlined: false,
)[#if name != none [ #strong[(#name).]] else [#strong[.]] #tags.join(h(3pt)) #h(2pt) #body]
#let definition = env([Definition])
#let lemma = env([Lemma])
#let proposition = env([Proposition])
#let theorem = env([Theorem])
#let corollary = env([Corollary])
#let conj = env([Conjecture])

#let proof(title: [Proof], body) = block(
  width: 100%,
  breakable: true,
  inset: (left: 10pt, y: 2pt),
  stroke: (left: 0.5pt + luma(170)),
)[#set par(justify: false); _#title._ #body #h(1fr) $square$]

#let remark(body) = block(width: 100%, breakable: true)[#set par(justify: false); _Remark._ #body]

// ---------------------------------------------------------------------------
// Title
// ---------------------------------------------------------------------------
#align(center)[
  #text(size: 18pt, weight: "bold")[Quantum Glitchy Period Finding]
  #v(-4pt)
  #text(size: 13pt)[Status report and open problems]
  #v(6pt)
  Searidang Pa \
  #text(size: 9.5pt)[based on summer 2019 research with Shelby Kimmel, Department of Computer Science, Middlebury College] \
  #text(size: 9pt, fill: luma(90))[Compiled October 2026 from the 2019 poster, simulation code and logs in this repository]
]

#v(6pt)
#block(inset: (x: 18pt))[
  #set text(size: 9.5pt)
  *Abstract.* Quantum period finding assumes an exactly periodic oracle. We study a periodic
  sequence corrupted by $k$ insertions or deletions ("glitches"). Each glitch shifts every later
  entry, so the corrupted function disagrees with every periodic function on a constant fraction of
  its domain. This note records the state of the 2019 project: the idealized model used throughout,
  the exact measurement distribution, the one-glitch theorem (probability near a good outcome at
  least $0.447$ times ideal), the multi-glitch conjecture, and the numerical evidence. The 2019
  materials (poster, code, logs) contain no written proofs. Material marked #new was reconstructed
  or added while compiling this note. It includes an identification of the one-glitch constant as
  $8 \/ (8 + pi^2)$, a Parseval reformulation of the multi-glitch problem, a crude bound with
  quadratic radius, and a family of configurations showing that the conjectured constant $1\/2$
  would be asymptotically sharp. Separately, an exact residue-counting argument proves total
  probability at least $1\/(2k+1)$ at the ideal outcomes themselves after any $k$ literal
  glitches, with a stronger bound approaching $1\/(k+1)$ when $k=o(N)$. A step-by-step proof
  and a worked example explain this total-mass guarantee; it does not resolve the stronger
  per-outcome conjecture or establish period recovery. A separate Fejér-kernel argument proves
  total probability at least $1\/3$ within distance $k$ of the rounded peaks for every $k>=1$,
  even when $r divides.not N$. A list of open problems closes the note.
]

#v(4pt)
#block(inset: (x: 18pt))[
  #set text(size: 9pt)
  *How to read the status labels.*
  #from2019 the statement or formula appears in the 2019 poster, code or logs.
  #claim2019 stated on the 2019 poster; no proof is in the repository.
  #new added while compiling this note.
  #proved a complete argument is given here.
  #leanchecked the general-$k$ arguments in @sec:total and @sec:nondividing have checked Lean proofs.
  #sketch the argument is outlined and the details still need to be checked.
  #numerical supported by computation only.
  #conjecture open.
  Statements marked #new were checked numerically with
  `uv run glitchy.py check` (@app:repro). They were prepared with AI assistance
  and have not been independently refereed.
]

#outline(depth: 1, indent: auto)

= Summary of the current state <sec:summary>

- *Model.* All 2019 analysis and simulation use one idealized _block-shift model_
  (@def:model). In it every glitch sits on a period boundary and shifts all later copies of the
  period by $plus.minus 1$. Its measurement distribution has an exact closed form (@lem:exact).
- *One glitch.* The poster states that for $r divides N$ the probability of an outcome within distance $1$
  of a good outcome is more than $0.447$ times the ideal probability. In the continuum limit we
  reconstruct the argument (@thm:one). The neighbourhood ratio is at least
  $1 - 4u(1-u) + (8 \/ pi^2) sin^2(pi u)$, and Jordan's inequality turns this into
  $8\/(8+pi^2) = 0.4477...$, matching the poster's constant. The sharp value is $0.6400$.
- *Several glitches.* The poster's "gruesome expression" is exactly the continuum neighbourhood
  mass (@sec:poster-formula). The 2019 target was to prove that it is at least $1\/2$
  (@conj:half). Exhaustive simulations for $k = 2, 3, 4$ give worst cases $0.772$, $0.614$,
  $0.689$ (@tab:worst). The 2019 grid search for $k = 4$ had found only $0.723$. The worst good
  outcome is $c = r\/2$ for $k = 1, 3$ but not for $k = 2, 4$.
- *Added in this note.* (a) The neighbourhood mass is $1$ minus the high-frequency energy of a
  unimodular step function, with an explicit formula in terms of jumps (@prop:parseval).
  (b) Radius $L = ceil(4(k+1)^2 \/ pi^2)$ already gives mass $>= 1\/2$ (@cor:quadratic), so the
  open problem is to reduce the radius to $O(k)$. (c) Error bounds for passing between the
  discrete and continuum models (@lem:continuum, @prop:transfer). (d) A packed family whose mass at
  radius $k$ decreases towards $1\/2$, reaching $0.521$ at $k = 200$ (@sec:packed). (e) Any bound
  uniform in $c$ needs radius at least $(k+1)\/2$ (@prop:radius).
- *Any number of glitches, total mass.* For $k$ literal insertions/deletions at arbitrary
  positions, with arbitrary inserted symbols, the total probability at the exact ideal
  outcomes is at least $r\/((k+1)r+I) >= 1\/(2k+1)$, where $I <= k$ is the number of
  surviving inserted entries (@thm:total). A second bound is $(1-I\/N)^2\/(k+1)$.
  Normalized Model I gives $1\/(k+1)$ (@prop:model-total). These also bound the probability
  within distance $k$, but are averaged statements, not lower bounds near each good outcome.
- *Nondivisible periods, total mass.* For every $k>=1$, literal edits have total probability
  at least $1\/3$ within distance $k$ of the rounded peaks, whether or not $r divides N$
  (@thm:nondividing). This is again a union-of-windows guarantee, not a per-outcome bound.
- *Not done.* No multi-glitch theorem uniform in $c$ at radius $O(k)$. No period-recovery theorem
  (@sec:recovery). No per-outcome analysis of $r divides.not N$. No public follow-up paper or manuscript exists as of October 2026
  (@sec:open, item 8).

= Standard period finding <sec:background>

Let $f : ZZ_N -> W$ have period $r$, with $r < sqrt(N)$. For most of this note we also assume
#footnote[The poster's figure 5 treats $r divides.not N$; nothing is proved there.]
$ r divides N, quad M := N\/r, quad f(x) = w_(x mod r) "with" w_0, ..., w_(r-1) "distinct". $ <eq:standing>
The algorithm prepares $N^(-1\/2) sum_x lr(|x chevron.r) lr(|f(x) chevron.r)$, applies the quantum Fourier
transform to the first register and measures. Outcome $y$ appears with probability
$ P_f (y) = 1/N^2 sum_(u in W) abs(sum_(x : f(x) = u) e^(2 pi i x y \/ N))^2 . $ <eq:pf>
For the ideal function, value $w_a$ is taken on ${a + m r : 0 <= m < M}$, and the inner sum is
$e^(2 pi i a y \/ N) sum_(m=0)^(M-1) e^(2 pi i m y \/ M) = M e^(2 pi i a y \/ N) bb(1)[M divides y]$. Hence
$ P_0 (y) = 1/r dot bb(1)[M divides y], quad G := {c M : 0 <= c < r}. $ <eq:ideal>
From a good outcome $y = c M$ we get $y\/N = c\/r$. Reducing the fraction gives $r \/ gcd(c, r)$,
and $Pr[gcd(c, r) = 1] = phi(r)\/r = Omega(1 \/ log log r)$. When $r divides.not N$, the peaks sit at the
integers nearest $c N \/ r$. Continued fractions then recover $c\/r$ because $r^2 <= N$. This is
where the assumption $r < sqrt(N)$ is used.

= The glitch model <sec:model>

== Informal picture

An _insertion_ adds one symbol to the periodic sequence and shifts everything after it by one
place. A _deletion_ removes one symbol and shifts everything after it back by one place.
@fig:glitch shows both. With period $3$, `abcabcabc…` becomes `abcadbcabc…` after an insertion and
`abcacabc…` after a deletion. After the glitch the sequence is still locally periodic, but its
phase has moved. A single glitch therefore changes $f$ on about half of the domain, so
perturbation bounds that count changed entries are useless (research prompt, Approach 2).

#figure(
  grid(
    columns: 2,
    gutter: 8pt,
    image("poster_figures/single_insertion.png", width: 100%),
    image("poster_figures/single_deletion.png", width: 100%),
  ),
  caption: [Poster figures 1--2. A periodic function before and after one insertion (left pair)
    and one deletion (right pair).],
) <fig:glitch>

== The block-shift model used in 2019

#definition(name: [Block-shift model, "Model I"], tags: (from2019,))[
  Assume @eq:standing. A _$k$-glitch configuration_ consists of
  - block boundaries $0 = m_0 <= m_1 <= dots.h <= m_k <= m_(k+1) = M$. Glitch $j$ sits between copy
    $m_j - 1$ and copy $m_j$ of the period;
  - signs $epsilon_1, ..., epsilon_k in {+1, -1}$, where $+1$ is an insertion and $-1$ a deletion.
  Put $s_1 = 0$, $s_(j+1) = s_j + epsilon_j$ (cumulative shift) and
  $B_j = {m_(j-1), ..., m_j - 1}$ for $j = 1, ..., k+1$. Copy $m in B_j$ of the period starts at
  $m r + s_j$, so value $w_a$ occupies $X_a = {a + m r + s_j : 1 <= j <= k+1, m in B_j}$, read
  modulo $N$. The Model I distribution is
  $ P_k (y)
    = 1/N^2 sum_(a=0)^(r-1) abs(sum_(x in X_a) e^(2 pi i x y \/ N))^2
    = r/N^2 abs(sum_(j=1)^(k+1) e^(2 pi i y s_j \/ N) sum_(m in B_j) e^(2 pi i m y \/ M))^2 . $ <eq:modelI>
  The second equality holds because $a$ enters only through the unimodular factor
  $e^(2 pi i a y \/ N)$.
] <def:model>

#remark[
  *Correspondence with the code.* @eq:modelI is `distribution` in `glitchy.py`, which is
  `probability_distr_with_insertion_deletion` from the 2019 scripts. There `ListK`
  $= (m_0, ..., m_(k+1))$, `Shift` $= (s_1, ..., s_(k+1))$, and `abs(z)**2/(mb*N)` equals
  $r abs(z)^2 \/ N^2$. Every 2019 run used equal signs, so $s_j = plus.minus (j - 1)$: the code picked
  "all insertions" or "all deletions" at random. For $k < r$ the points of each $X_a$ are distinct
  modulo $N$, so $sum_y P_k (y) = 1$, which `glitchy.py check` confirms.
]

#remark[
  *What Model I idealizes.* Take the _literal model_: an infinite $r$-periodic sequence, $k$ edits
  at arbitrary positions, arbitrary inserted symbols, and keep the first $N$ symbols. Model I
  differs from it in three ways.
  (i) A glitch strictly inside a period, between residues $b - 1$ and $b$, moves only the classes
  $a >= b$ of that copy. Model I moves the whole copy.
  (ii) The inserted symbols are dropped.
  (iii) Model I keeps exactly $M$ copies of every value and wraps positions outside $[0, N)$
  around. The literal model truncates instead.
  Each effect changes at most $O(k)$ points of each class. @prop:transfer bounds the resulting error.
  The 2019 materials do not settle several questions raised in the research prompt. The bounds
  are meant to hold for every placement: the poster's conjecture is worst-case, and the code
  enumerated placements exhaustively or at random. Inserted values play no role in Model I.
  Adjacent glitches are allowed ($m_j = m_(j+1)$ gives an empty block), but two glitches inside
  one period are not modelled. $k$ counts operations, and the net displacement is
  $abs(s_(k+1)) <= k$. Mixed sign patterns are allowed by @eq:modelI. The 2019 runs did not simulate them; see
  @sec:signs.
]

= The distribution near good outcomes <sec:near>

#lemma(name: [Exact form near good outcomes], tags: (from2019, proved))[
  Write $y = c M + h$ with $0 <= c < r$ and $abs(h) < M\/2$, and let $omega = e^(2 pi i \/ r)$. Then
  $ R_c (h) := (P_k (c M + h)) / (P_0 (c M))
    = 1/M^2 abs(sum_(j=1)^(k+1) omega^(c s_j) e^(2 pi i h s_j \/ N) D_j (h))^2,
    quad D_j (h) := sum_(m in B_j) e^(2 pi i m h \/ M), $ <eq:ratio>
  and $abs(D_j (h)) = abs(sin(pi h L_j \/ M) \/ sin(pi h \/ M))$, where $L_j = m_j - m_(j-1)$.
] <lem:exact>
#proof[
  Substitute $y = c M + h$ into @eq:modelI. Then $e^(2 pi i m y \/ M) = e^(2 pi i m h \/ M)$ because
  $m c in ZZ$, and $e^(2 pi i y s_j \/ N) = omega^(c s_j) e^(2 pi i h s_j \/ N)$. Divide by
  $P_0 (c M) = 1\/r$ and use $r^2 \/ N^2 = 1\/M^2$.
]

For a radius $L >= 0$ define the _neighbourhood ratio_ and its worst case:
$ S_c (L) := sum_(abs(h) <= L) R_c (h) = (Pr[Y "within" L "of" c M]) / (Pr_"ideal" [Y = c M]),
  quad S^* (L) := min_(0 <= c < r) S_c (L). $ <eq:nbhd>
The 2019 choice is $L = k$. The quantity $S_k = sum_(y in G_k) P_k (y)$ in the research prompt
equals $(1\/r) sum_c S_c (k) >= S^* (k)$ when $k < M\/2$. Except for the explicitly averaged
result below, the subsequent neighbourhood bounds are per $c$ and therefore also bound $S_k$.

== General $k$ literal glitches: a total-mass bound <sec:total>

*Answer in brief.* After any $k$ literal glitches, the probability of landing in the _set_ of
ideal outcomes is at least $1\/(2k+1)$. If $k$ is small compared with $N$, a second bound is
close to $1\/(k+1)$. Both are exact finite-$N$ guarantees. They require neither the block-shift
idealization nor a continuum limit.

The key is to count entries by two labels: their *output symbol* $u$ and their *input residue*
$b = x mod r$, where $x$ is the entry's position. Call each pair $(u,b)$ a _box_. Averaging over all ideal Fourier outcomes turns
interfering amplitudes into a sum of squared box counts. Glitches limit the number of occupied
boxes, and spreading $N$ entries among only a few boxes forces that sum of squares to be large.

#block(breakable: false)[
#theorem(name: [General literal glitches, total mass], tags: (new, proved, leanchecked))[
  Let $N, r$ be positive integers with $r divides N$, and start from the infinite periodic
  sequence with symbols $w_0, ..., w_(r-1)$, which need not be distinct. Perform $A$
  single-symbol insertions and $D$ single-symbol deletions, $k = A + D >= 0$ edits in all, in
  any order, at arbitrary positions, and with arbitrary inserted symbols. Keep the first $N$
  entries of the resulting sequence, obtaining $f_k$. Let $I$ be the number of surviving
  inserted entries in this prefix, so $0 <= I <= min(A,N)$.

  For the measurement distribution $P_(f_k)$, define
  $ G := {c N\/r : 0 <= c < r}, quad
    G_k := {y in ZZ_N : min_(g in G) d_N (y,g) <= k}, $
  where $d_N (y,g) := min_(m in ZZ) abs(y-g-m N)$ is circular distance modulo $N$. Then
  $ Pr[Y in G_k] >= Pr[Y in G] >= r/((k+1)r + I)
    >= r/((k+1)r + k) >= 1/(2k+1). $ <eq:total>
  Independently,
  $ Pr[Y in G] >= (N-I)^2/((k+1)N^2)
    = (1-I\/N)^2/(k+1). $ <eq:retained>
  The larger of these two lower bounds may be used. There is no assumption that $k < r$,
  $k < N$, or $r < sqrt(N)$. Overlap between neighbourhoods in $G_k$ causes no problem:
  the proof already bounds the smaller set $G$.
] <thm:total>
]

#figure(
  image("report_figures/general_k_proof-1.svg", width: 70%),
  caption: [Proof map. A box means one pair (output symbol, input residue modulo $r$).
    The count argument controls the _total_ probability on $G$, not its division among the
    individual good outcomes. Mermaid source: #text(font: "DejaVu Sans Mono", size: 8pt)[report_figures/general_k_proof.md].],
) <fig:total-proof>

#pagebreak(weak: true)
#proof[
  *Step 1: Express the Fourier probability using box counts.*

  For any function $f : ZZ_N -> W$, let
  $ n_(u b) := abs({x in {0, ..., N-1} : f(x)=u, x equiv b thick (mod r)}). $
  These counts are nonnegative and $sum_(u,b) n_(u b)=N$. At an ideal outcome $y=c N\/r$,
  the phase depends on $x$ only through its residue $b$:
  $ e^(2 pi i x y\/N) = e^(2 pi i c x\/r) = e^(2 pi i c b\/r). $
  Thus @eq:pf becomes
  $ P_f (c N\/r) = 1/N^2 sum_u abs(sum_(b=0)^(r-1) n_(u b) e^(2 pi i c b\/r))^2. $
  Sum over $c$. The identity
  $ sum_(c=0)^(r-1) e^(2 pi i c(b-b')\/r)
    = cases(r & "if" b=b', 0 & "if" b != b') $
  removes every cross-term between different residues. Therefore
  $ Pr[Y in G]
    &= 1/N^2 sum_u sum_(b=0)^(r-1) sum_(b'=0)^(r-1) n_(u b) n_(u b')
       sum_(c=0)^(r-1) e^(2 pi i c(b-b')\/r) \
    &= r/N^2 sum_u sum_(b=0)^(r-1) n_(u b)^2. $ <eq:residue-mass>

  In words: the Fourier calculation has become an ordinary counting problem. Two entries in
  the same box always reinforce one another in this sum; different boxes do not interfere
  after averaging over $c$.

  *Step 2: Bound the number of occupied boxes.*

  Give every entry a _tag_ recording where it came from. An original entry's tag is its
  original index $t >= 0$; the $j$-th insertion gets the tag $-j$. Write $tau(x)$ for the tag
  at position $x$ of the edited sequence. Its symbol is
  $ f_k (x) = cases(w_(tau(x) mod r) & "if" tau(x) >= 0, v(tau(x)) & "if" tau(x) < 0) $
  where $v$ assigns the arbitrary inserted symbols. Insertions and deletions preserve the
  order of all surviving original entries.

  Consider an original entry that survives at position $x$, so that $t = tau(x) >= 0$. Among
  the original $t$ entries preceding it, let $d(t)$ have been deleted. Let $i(t)$ surviving
  inserted entries precede it in the final sequence. Counting the entries before it gives
  $ x = t-d(t)+i(t), quad s := x-t = i(t)-d(t). $
  Since $0 <= i(t) <= I <= A$ and $0 <= d(t) <= D$,
  $ -D <= s <= A. $ <eq:shift-range>
  This interval contains only $A+D+1=k+1$ integers. Hence there are at most $k+1$ possible
  displacements, regardless of the order or locations of the edits.

  An original entry with symbol $w_a$ has $t equiv a thick (mod r)$. Its final residue is
  $x equiv a+s thick (mod r)$. For each of the $r$ original symbols, at most $k+1$ residues are
  therefore possible. So retained original entries occupy at most $(k+1)r$ boxes.
  Each of the $I$ surviving inserted entries can add at most one further box. In total,
  $ q := abs({(u,b):n_(u b)>0}) <= (k+1)r+I =: Q. $ <eq:box-count>

  This reasoning includes adjacent edits, several edits inside one period, an insertion
  subsequently deleted, and inserted symbols equal to original symbols. Reusing an existing
  box can only reduce the count. Restricting to the first $N$ entries also cannot add boxes.

  *Step 3: Few boxes force a large sum of squares.*

  Cauchy--Schwarz, applied just to the $q$ nonempty boxes, gives
  $ N^2 = (sum_(u,b) n_(u b))^2 <= q sum_(u,b) n_(u b)^2
    <= Q sum_(u,b) n_(u b)^2. $
  Equivalently, the sum of squares is smallest when entries are spread as evenly as possible.
  Substitute into @eq:residue-mass:
  $ Pr[Y in G] >= r/N^2 dot N^2/Q = r/((k+1)r+I). $
  Since $I <= k$ and $r >= 1$,
  $ r/((k+1)r+I) = 1/(k+1+I\/r) >= 1/(2k+1). $
  Finally, $G subset.eq G_k$, proving @eq:total. The same conclusion holds for ordinary,
  non-wrapped distance-$k$ neighbourhoods, since they too contain $G$.

  *Step 4: A second bound, ignoring inserted entries.*

  Let $m_(u b)$ count only retained original entries in each box. Then
  $0 <= m_(u b) <= n_(u b)$, their sum is $N-I$, and they occupy at most $(k+1)r$ boxes.
  A second application of Cauchy--Schwarz yields
  $ sum_(u,b) n_(u b)^2 >= sum_(u,b) m_(u b)^2
    >= (N-I)^2/((k+1)r). $
  This is also valid when $I=N$, when its right-hand side is zero.
  Multiplying by $r\/N^2$ in @eq:residue-mass proves @eq:retained.
]

#block(breakable: false)[
*Worked example: three mixed glitches.* Start with #raw("abcabcabc…"), insert #raw("d") just before
original index $3$, insert #raw("e") just before original index $6$, and delete original index $9$.
Keep $N=15$ entries. Here $r=3$, $k=3$, $A=2$, $D=1$, and $I=2$; all positions and indices
start at zero. In the edited sequence the three edits happen at positions $3$, $7$, and $11$,
and the inserted entries get the tags $-1$ and $-2$.
@fig:residue-example tracks the surviving originals and builds the count table.

#figure(
  image("report_figures/general_k_example.svg", width: 100%),
  caption: [A literal three-edit example, not a block-shift approximation. Inserted entries
    have negative tags and no displacement. The deletion removes original index $9$ and changes
    the later displacement $s$ from $2$ back to $1$. The table counts each final entry exactly
    once.],
) <fig:residue-example>
]

The table has $sum n_(u b)=15$ and $sum n_(u b)^2=25$, so the identity gives
$ Pr[Y in {0,5,10}] = 3/15^2 dot 25 = 1/3
  >= 3/((3+1)3+2) = 3/14 >= 1/7. $
The middle quantity is the theorem's box-count bound; the last is its simpler universal bound.
Repeated displacement values do not need new boxes, so an actual configuration can do better.

*Useful special cases.*
- *No edits:* $k=I=0$, and the bound is $1$, as it should be.
- *Deletions only (or no inserted entry survives in the prefix):* $I=0$ gives
  $Pr[Y in G] >= 1\/(k+1)$.
- *Few edits compared with the sample length:* if $k <= N$, then @eq:retained gives
  $Pr[Y in G] >= (1-k\/N)^2\/(k+1)$. In particular, for $k=o(N)$ the guarantee is
  $(1-o(1))\/(k+1)$. The condition $k <= N$ is needed for this substitution, but not for the
  theorem itself.
- *Two glitches:* setting $k=2$ recovers the earlier result
  $Pr[Y in G_2] >= Pr[Y in G] >= r\/(3r+I) >= r\/(3r+2) >= 1\/5$,
  together with $Pr[Y in G] >= (N-I)^2\/(3N^2)$.

#proposition(name: [General glitches in normalized Model I], tags: (new, proved, leanchecked))[
  For any normalized Model I configuration (in particular, any configuration with $k<r$),
  $ Pr_"Model I" [Y in G] >= 1/(k+1). $ <eq:model-total>
  This bound is attained at the exact set $G$ when the $k+1$ blocks have equal lengths and
  their shifts are distinct modulo $r$.
] <prop:model-total>
#proof[
  Let $ell_j=(m_j-m_(j-1))\/M$ and group together lengths with the same shift residue:
  $ p_b := sum_(j:s_j equiv b thick (mod r)) ell_j. $
  The $p_b$ are nonnegative, sum to $1$, and at most $k+1$ are nonzero. At $h=0$,
  @eq:ratio and root-of-unity orthogonality give
  $ Pr_"Model I" [Y in G]
    = 1/r sum_(c=0)^(r-1) abs(sum_j ell_j e^(2 pi i c s_j\/r))^2
    = sum_(b=0)^(r-1) p_b^2 >= 1/(k+1). $
  For equality, take all insertions, $r>k$, and $M$ divisible by $k+1$, with
  $ell_j=1\/(k+1)$. The shifts $0,1,...,k$ are distinct modulo $r$, and the $k+1$ nonzero
  grouped lengths are all $1\/(k+1)$.
]

#remark[
  *What is proved, and what remains open.* The two probabilities below are different:
  - *Total mass:* @thm:total guarantees probability $Omega(1\/(k+1))$ on the union of good
    outcomes, even before adding their neighbours. Thus the weaker total-mass target in the
    research prompt is proved for every $k$, including literal mixed-sign edits.
  - *Per-outcome mass:* @conj:half asks for a constant lower bound on
    $S_alpha^"cont"(k)$ for _every_ $alpha$. Our proof sums over $c$ to cancel cross-terms;
    it does not bound the summands separately.

  Even exact ideal outcomes can have very different probabilities. For example, in Model I
  with one insertion, two equal blocks, and even $r>1$, the amplitude at $c=r\/2$ is
  $(1\/2)-(1\/2)=0$, although the total mass on $G$ is $1\/2$. This example concerns radius
  zero, not a failure of the radius-$k$ conjecture.

  The equality example in @prop:model-total shows that no $k$-independent positive lower
  bound is possible on the _exact_ set $G$ in Model I. It does not rule out a constant bound
  on the larger set $G_k$. Finally, total mass includes uninformative outcomes such as $c=0$:
  neither theorem supplies the mass near coprime $c$, so period recovery remains a separate
  problem (@sec:recovery).
]

*Machine-checked companion.* The proofs for this subsection are in #raw("formal/"), using
Lean 4 and mathlib v4.19.0. The formalization starts with actual insertion/deletion operations
on an infinite tagged sequence, proves the displacement and box-count bounds, proves character
orthogonality, and connects the counts to the complex-exponential Fourier measurement formula.
It does not assume the key counting or Fourier identities as hypotheses.

The checked results include both bounds in @thm:total, their neighbourhood and finite-size
corollaries, and the Model I bound and sharpness in @prop:model-total. The figure's actual edit
history, count table, and mass $1\/3$ are also checked, not just its final fraction. There are
no #raw("sorry") placeholders or added axioms; the axiom audit reports only Lean's standard
logical foundations.

The formal alphabet is finite, which suffices because a finite period and finitely many edits
use only finitely many symbols. The Fourier measurement formula is taken as the definition;
a quantum circuit implementation is not verified. Lean also covers @sec:nondividing below;
it does not certify the other continuum results, the transfer sketch, or the open
per-outcome and recovery claims. The Lean names follow this note's symbols. See
#raw("formal/README.md") for the theorem-by-theorem map and a table matching each symbol to
its Lean name, and @app:repro for commands.

== Nondivisible periods: a constant neighbourhood bound <sec:nondividing>

The exact box-count identity requires $r divides N$; substituting rounded peak locations into
that identity is not valid. Nevertheless, the _total neighbourhood mass_ has a stronger bound
that does not require divisibility. The following argument uses the actual literal edited
prefix, not the legacy nondivisible Model I formula with $M = ceil(N\/r)$.

#theorem(name: [Literal glitches without a divisibility assumption], tags: (new, proved, leanchecked))[
  Let $N >= 1$ and $r >= 2$. Start from an infinite $r$-periodic sequence, perform $A$
  insertions and $D$ deletions at arbitrary positions, $k = A + D >= 1$ edits in all, and keep
  the first $N$ entries, obtaining $f_k$. Inserted symbols are arbitrary. Let $g_c$ be any nearest integer to $c N\/r$, read modulo
  $N$, and put
  $ G_k := {y in ZZ_N : min_(0 <= c < r) d_N (y,g_c) <= k}. $
  For the Fourier measurement distribution @eq:pf,
  $ Pr[Y in G_k] >= 133/358 > 1/3 >= 1/(2k+1). $ <eq:nondividing>
  This holds whether or not $r divides N$, for all edit placements and mixed signs.
  There is no additional assumption on $k\/N$, $k r\/N$, or $r\/N$, and overlapping windows
  are counted only once. The constant $133\/358$ is not claimed to be sharp.
] <thm:nondividing>

#proof[
  *Step 1: Relate Fourier moments to equal-symbol pairs.*

  For any integer $h >= 1$, define the circular matching fraction
  $ C_h := 1/N sum_(x=0)^(N-1) bb(1)[f_k (x)=f_k ((x+h) mod N)]. $
  Expanding @eq:pf and summing the $N$ roots of unity gives
  $ sum_(y=0)^(N-1) P_(f_k)(y) cos(2 pi h y\/N) = C_h. $ <eq:matching-moment>

  Use the tags $tau$ from Step 2 of @sec:total. Call $x$ a _good start_ for the lag $h$ if
  positions $x$ and $x+h$ both hold original entries whose tags differ by exactly $h$, that is,
  $tau(x) >= 0$ and $tau(x+h) = tau(x)+h$. When $r divides h$ the two tags are congruent
  modulo $r$, so a good start has $f_k (x) = f_k (x+h)$.

  For every length $L$, at most $k h + A$ of the starts $x < L$ are not good. This follows by
  induction on the edit history. In the unedited sequence every start is good. An insertion
  destroys at most $h$ good pairs, those that straddle it, and adds one start. A deletion
  destroys at most $h+1$ good pairs, those that contain or straddle the deleted entry, and
  removes one start. So each insertion adds at most $h+1$ starts that are not good, each
  deletion at most $h$, and $A(h+1) + D h = k h + A$.

  Take $L = N-h$, so that no pair wraps around, and count the $h$ wrapped starts as failures.
  For every multiple $h$ of $r$ this gives
  $ C_h >= 1 - ((k+1)h + A)/N. $ <eq:matching-bound>
  For $h >= N$ the right-hand side is nonpositive, so the inequality is automatic.

  *Step 2: Separate the rounded good windows from their complement.*

  If $N < (2k+1)r$, the rounded windows cover every integer outcome. Indeed, an integer
  is within $N\/(2r)$ of some real peak and within $N\/(2r)+1\/2 < k+1$ of its rounded
  peak, hence within $k$. Thus this case is immediate. Otherwise $N >= (2k+1)r$; put
  $ t := (2N)/(r(2k+1)) >= 2, quad p := Pr[Y in G_k], quad
    theta_y := min_(m in ZZ) abs(r y\/N - m) in [0,1/2]. $
  If $y$ is outside $G_k$, then its integer distance from every $g_c$ is at least $k+1$.
  Since $d_N (g_c,c N\/r) <= 1\/2$, its distance from every real peak is at least $k+1\/2$.
  Consequently $theta_y >= 1\/t$.

  For every integer $ell >= 1$, @eq:matching-bound, $r>=2$, and $A<=k$ imply
  $ C_(ell r) >= 1 - (2(k+1)ell + k)/((2k+1)t)
    >= 1 - (4ell+1)/(3t). $ <eq:normalised-matching>
  The last inequality is equivalent to $(k-1)(2ell-1) >= 0$.

  *Step 3: Average a nonnegative Fourier kernel.*

  For an integer $T>=2$, the Fejér kernel is
  $ F_T (theta) := 1/T abs(sum_(j=0)^(T-1) e^(2 pi i j theta))^2
    = 1 + 2 sum_(ell=1)^(T-1) (1-ell\/T) cos(2 pi ell theta). $
  It satisfies $0 <= F_T <= T$. From @eq:matching-moment and @eq:normalised-matching,
  $ E_T := sum_y P_(f_k)(y) F_T (r y\/N), quad
    E_T/T >= 1 - (4(T^2-1))/(9t T) - (T-1)/(3t T). $ <eq:fejer-expectation>
  Here $sum_(ell=1)^(T-1) ell(1-ell\/T)=(T^2-1)\/6$ and
  $sum_(ell=1)^(T-1)(1-ell\/T)=(T-1)\/2$.

  If $F_T\/T <= q<1$ outside $G_k$, then
  $ E_T/T <= p+(1-p)q, quad p >= (E_T\/T-q)/(1-q). $ <eq:fejer-mass>

  *Step 4: Choose the kernel size, using two elementary estimates.*

  *Short kernel, $2<=t<=6$.* Take $T=2$. Then $F_2 (theta)\/2=cos^2(pi theta)$, so we may take
  $q=cos^2(pi\/t)$. @eq:fejer-expectation gives $E_2\/2 >= 1-5\/(6t)$, hence
  $ p >= 1 - 5/(6t sin^2(pi\/t)) >= 3/8. $
  To check the last step, $t sin^2(pi\/t)>=4\/3$ throughout this range:
  for $2<=t<=3$ use $sin^2(pi\/t)>=3\/4$; for $3<=t<=4$ use
  $sin^2(pi\/t)>=1\/2$; for $4<=t<=6$, concavity of sine on $[0,pi\/4]$ gives
  $sin(pi\/t)>=2 sqrt(2)\/t$, so $t sin^2(pi\/t)>=8\/t>=4\/3$.

  *Long kernel, $t>6$.* Take $T=ceil(3t\/4)$. Using $T<=3t\/4+1$ in @eq:fejer-expectation gives
  $ E_T/T >= 2/3 - 7/(9t) >= 29/54. $
  Outside $G_k$, the geometric-sum formula gives
  $ F_T (theta_y)\/T <= 1/(T^2 sin^2(pi theta_y))
    <= 1/(T^2 sin^2(pi\/t)) <= 64/243. $
  Indeed, concavity of sine on $[0,pi\/3]$ gives
  $sin(pi\/t)>=(3 sqrt(3)\/2)\/t$, and $T>=3t\/4$.
  Now @eq:fejer-mass yields
  $ p >= (29\/54-64\/243)/(1-64\/243) = 133/358 > 1/3. $
  Since $3\/8>133\/358$, both cases prove the stated bound.
]

#remark[
  *Scope.* This proves the original $1\/(2k+1)$ _total-mass_ guarantee for $r divides.not N$
  and $k>=1$, without an extra finite-size condition. It does not extend the box-count
  identity, bound the mass at just the rounded peaks, or bound each peak separately.
  In particular, it still does not establish period recovery or the informative-outcome
  conjecture. The $k=0$ claim with lower bound $1$ cannot be carried over to the nearest
  integers: an unedited nondivisible period already has Fourier leakage.

  *Machine-checked companion.* #raw("formal/Glitchy/Nondividing.lean") proves this theorem
  end to end for the same literal edit histories and complex-exponential measurement
  formula. Matching moments, pair counts, rounded-window separation, the squared geometric
  sum, and both analytic regimes are proved, not supplied as hypotheses. Every rounding-tie
  choice is allowed, and overlapping windows are counted as a set. The axiom audit uses
  only Lean's standard logical foundations, with no placeholders or added axioms.

  The Lean proof follows Steps 1--4 above with the same symbols: the good-start induction, the
  coverage case $N<(2k+1)r$, and the short and long kernels. Deterministic numerical checks
  remain in #raw("glitchy.py check"). See #raw("formal/README.md") for the declaration map.
]

== Continuum limit

Let $u_j = m_j \/ M in [0, 1]$, $ell_j = u_j - u_(j-1)$ and $alpha = c\/r$. Define the _phase function_
and its Fourier coefficients
$ phi_alpha (u) &= e^(2 pi i alpha s_j) quad "for" u in [u_(j-1), u_j), \
  a_h &:= integral_0^1 phi_alpha (u) e^(2 pi i h u) dif u
  = sum_(j=1)^(k+1) e^(2 pi i alpha s_j) e^(i pi h (u_(j-1) + u_j)) (sin(pi h ell_j)) / (pi h), $ <eq:ah>
where the $h = 0$ term reads $ell_j$. The continuum neighbourhood mass is
$S_alpha^"cont" (L) = sum_(abs(h) <= L) abs(a_h)^2$.

#lemma(name: [Continuum approximation], tags: (new, proved))[
  For all $c$ and $h$, with $epsilon_h := pi abs(h) (1 + 2k\/r) \/ M$,
  $ abs(sqrt(R_c (h)) - abs(a_h)) <= epsilon_h, quad
    abs(S_c (L) - S_(c\/r)^"cont" (L)) <= sum_(abs(h) <= L) (2 epsilon_h + epsilon_h^2)
    = O(L^2 (1 + k\/r) \/ M) . $
] <lem:continuum>
#proof[
  $sqrt(R_c (h))$ is the modulus of $A := M^(-1) sum_j omega^(c s_j) e^(2 pi i h s_j \/ N) D_j (h)$.
  The sum $M^(-1) D_j (h)$ is a left Riemann sum, with step $1\/M$, for
  $integral_(u_(j-1))^(u_j) e^(2 pi i h u) dif u$. The integrand has derivative of size at most $2 pi abs(h)$,
  so each cell contributes error at most $pi abs(h) \/ M^2$, and the total over all blocks is at most
  $pi abs(h) \/ M$. In addition, $abs(e^(2 pi i h s_j \/ N) - 1) <= 2 pi abs(h) k \/ N = 2 pi abs(h) k \/ (r M)$,
  and $sum_j abs(D_j) \/ M <= 1$. So $abs(A - a_h) <= epsilon_h$. Finally,
  $abs(abs(A)^2 - abs(a_h)^2) <= epsilon_h (2 + epsilon_h)$ because $abs(a_h) <= 1$.
]

== The poster's "gruesome expression" <sec:poster-formula>

The poster states the multi-glitch problem with $l = k$ glitches at $0 < x_1 < dots.h < x_l < t$
and $0 < alpha < 1$. Write $t = 1$, $x_j = u_j$, and take all insertions, so $s_j = j - 1$. Its
quantities are the continuum Fourier data of @eq:ah:
$ A_0 = abs(a_0)^2, quad A_(h j) = (sin(pi h ell_j)) / (pi h), quad S_h = abs(a_h)^2 + abs(a_(-h))^2, quad
  S = A_0 + sum_(h=1)^l S_h = S_alpha^"cont" (l). $
Expanding the squares gives the poster's formulas:
$ A_0 &= sum_(j=1)^(l+1) ell_j^2 + 2 sum_(1 <= j' < j <= l+1) ell_j ell_(j') cos(2 pi alpha (j - j')), \
  S_h &= 2 sum_(j=1)^(l+1) A_(h j)^2
       + 4 sum_(1 <= j' < j <= l+1) A_(h j) A_(h j') cos(2 pi alpha (j - j')) cos(pi h (u_(j-1) + u_j - u_(j'-1) - u_(j'))) . $
The second line uses $cos(theta + psi) + cos(theta - psi) = 2 cos theta cos psi$. Three
transcription slips on the poster are corrected here. The first sum in $A_0$ runs over
$j = 1, ..., l+1$, not $0, ..., l$. The last term in the second cosine is $x_(j')$, not $x_j$. The
prefactor $t$ in $A_(h j)$ is consistent only with the normalization $t = 1$. For mixed signs,
replace $j - j'$ by $s_j - s_(j')$. All-deletion configurations give the same $S$, because $S$ is
even in $alpha$.

#conj(name: [Poster, 2019], tags: (conjecture, from2019))[
  For every $k >= 1$, every $0 < u_1 < dots.h < u_k < 1$, all-insertion or all-deletion signs, and
  every $alpha in [0, 1)$, we have $S_alpha^"cont" (k) >= 1\/2$. By @lem:continuum this gives
  $S^* (k) >= 1\/2 - O(k^2 (1 + k\/r) \/ M)$ in Model I.
] <conj:half>

The poster interprets this as "even in the worst case, the cost to the runtime would scale like
the number of glitches". The $2k + 1$ candidates around each outcome must be tried, but only $O(1)$
more samples are needed (@sec:recovery). The target $S_k >= c\/(k+1)$ in the research prompt is
weaker than @conj:half: @sec:total proves it in normalized Model I with constant $1$ and
in the literal model with constant $1\/2$, since $1\/(2k+1) >= (1\/2)\/(k+1)$.

= A Parseval reformulation <sec:parseval>

#proposition(name: [Energy and jump formula], tags: (new, proved))[
  Let $phi_j := e^(2 pi i alpha s_j)$ be the value of $phi_alpha$ on piece $j$. Set $u_0 := 0$ and
  define the jumps
  $ Delta_j := phi_j - phi_(j+1) quad (1 <= j <= k), quad quad Delta_0 := phi_(k+1) - phi_1 quad ("wrap-around at" u = 0 equiv 1). $
  + $sum_(h in ZZ) abs(a_h)^2 = 1$, so $S_alpha^"cont" (L) = 1 - sum_(abs(h) > L) abs(a_h)^2$.
  + For $h != 0$, $a_h = (2 pi i h)^(-1) sum_(j=0)^k Delta_j e^(2 pi i h u_j)$.
  + $abs(Delta_j) = 2 abs(sin pi alpha)$ for every glitch $j >= 1$, and $abs(Delta_0) = 2 abs(sin(pi alpha s_(k+1)))$.
  Consequently
  $ 1 - S_alpha^"cont" (L) = 1/(4 pi^2) sum_(abs(h) > L) 1/h^2 abs(sum_(j=0)^k Delta_j e^(2 pi i h u_j))^2 . $ <eq:jump>
] <prop:parseval>
#proof[
  (1) Parseval's identity for $phi_alpha in L^2 [0, 1)$ with $abs(phi_alpha) = 1$. (2) Integrate
  each piece, $integral_(u_(j-1))^(u_j) e^(2 pi i h u) dif u = (e^(2 pi i h u_j) - e^(2 pi i h u_(j-1))) \/ (2 pi i h)$,
  and collect the coefficient of each $e^(2 pi i h u_j)$. Use $e^(2 pi i h dot 1) = 1$ to merge the
  endpoint terms into $Delta_0$. (3) $Delta_j = phi_j (1 - e^(2 pi i alpha epsilon_j))$ with
  $epsilon_j = plus.minus 1$, and $s_1 = 0$.
]

So the multi-glitch problem is purely about step functions. $S$ is the fraction of the energy
of a unimodular step function with $k + 1$ jumps that lies at frequencies $abs(h) <= L$. The $k$
glitch jumps have equal size, and there is one wrap-around jump. All destructive interference is
in the exponential sum $sum_j Delta_j e^(2 pi i h u_j)$ in @eq:jump.

#corollary(name: [Quadratic radius suffices], tags: (new, proved))[
  $S_alpha^"cont" (L) >= 1 - 2(k+1)^2 \/ (pi^2 L)$ for every configuration and every $alpha$. In
  particular, $S_alpha^"cont" (L) >= 1\/2$ for $L >= 4(k+1)^2 \/ pi^2$.
] <cor:quadratic>
#proof[
  In @eq:jump use $abs(sum_j Delta_j e^(2 pi i h u_j)) <= sum_j abs(Delta_j) <= 2(k+1)$ and
  $sum_(abs(h) > L) h^(-2) <= 2\/L$.
]

This bound gives nothing at $L = k$, because $2(k+1)^2 \/ (pi^2 k) > 1$ for every $k >= 1$. The
crude estimate ignores cancellation between jumps. Exploiting that cancellation, for example
with large-sieve-type inequalities for well-separated jump positions combined with merging of
clustered jumps, is the obvious route to radius $O(k)$.

= One glitch <sec:one>

#theorem(name: [One glitch], tags: (claim2019, new, proved))[
  Let $k = 1$, with the glitch at $u in (0, 1)$, either sign, and any $alpha in [0, 1)$. In the continuum
  limit,
  $ S_alpha^"cont" (1) = 1 - (1 - cos 2 pi alpha) D(u), quad D(u) := 2u(1-u) - 4/pi^2 sin^2(pi u) >= 0 . $ <eq:onefull>
  Hence
  $ S_alpha^"cont" (1) >= g(u) := 1 - 4u(1-u) + 8/pi^2 sin^2(pi u), $ <eq:onebound>
  with equality when $alpha = 1\/2$, and $g(u) >= 8\/(8 + pi^2) = 0.44768...$ for all $u$.
  The sharp constant is $min_u g(u) = 0.6400$, attained at $u approx 0.198$ and its mirror
  $u approx 0.802$.
] <thm:one>
#proof[
  With $theta = 2 pi alpha epsilon_1$ we have $Delta_1 = 1 - e^(i theta)$ and $Delta_0 = e^(i theta) - 1$.
  So $abs(Delta_1 e^(2 pi i h u) + Delta_0)^2 = abs(1 - e^(i theta))^2 abs(e^(2 pi i h u) - 1)^2
  = 8 (1 - cos theta) sin^2(pi h u)$. By @eq:jump,
  $1 - S = 2 (1 - cos theta) pi^(-2) sum_(abs(h) >= 2) sin^2(pi h u) \/ h^2$. The classical identity
  $sum_(h >= 1) sin^2(pi h u) \/ h^2 = pi^2 u (1 - u) \/ 2$ for $u in [0, 1]$ gives
  $sum_(abs(h) >= 2) sin^2(pi h u)\/h^2 = pi^2 u(1-u) - 2 sin^2(pi u)$. Since $cos theta = cos 2 pi alpha$,
  @eq:onefull follows.

  For $D >= 0$, use the symmetry $u <-> 1 - u$ and take $u <= 1\/2$. If $u <= 1\/3$, then
  $sin(pi u) <= pi u <= pi sqrt(u(1-u)\/2)$. If $1\/3 <= u <= 1\/2$, then
  $sin(pi u) <= 1 <= pi\/3 <= pi sqrt(u(1-u)\/2)$. Hence the worst $alpha$ is $1\/2$, which gives $g$.

  Jordan's inequality $sin(pi u) >= 2u$ on $[0, 1\/2]$ gives
  $g(u) >= 1 - 4u + (4 + 32\/pi^2) u^2$. This quadratic is minimized at $u = 2\/(4 + 32\/pi^2)$,
  where it equals $1 - 4\/(4 + 32\/pi^2) = 8\/(8 + pi^2)$.
]

#remark[
  *Relation to the poster.* The poster says the probability is "greater than $0.447$" times the
  no-glitch probability. The constant $8\/(8 + pi^2) = 0.44769$ agrees in every printed digit and
  comes from the natural argument above: reduce to $alpha = 1\/2$, then apply Jordan's inequality.
  We therefore believe this was the 2019 argument, but the original proof is not in the
  repository. Jordan's inequality is the only lossy step. The sharp continuum value is $0.6400$.
  An exhaustive search over the glitch block in Model I gives $0.6403$ for $(N, r) = (3000, 50)$ and
  $0.6400$ for $(20000, 20)$. The worst case $alpha = 1\/2$ occurs only for even $r$ and $c = r\/2$.
  Odd $r$ does slightly better.
  The theorem is stated in the continuum limit. In Model I, @lem:continuum gives
  $S^* (1) >= 0.4476 - O((1 + 1\/r)\/M)$. Since the bound holds for each $c$ separately, the
  probability of landing within $1$ of a good outcome $c M$ with $gcd(c, r) = 1$ is at least
  $0.447 dot phi(r)\/r$.
]

= Numerical evidence <sec:numerics>

== The simulator

#figure(
  table(
    columns: (auto, 1fr),
    stroke: none,
    align: (left, left),
    table.hline(),
    table.header[*Command*][*What it does*],
    table.hline(stroke: 0.5pt),
    [`glitchy.py plot`], [Plots, for one placement (given or random), the probability near every good outcome relative to ideal: the neighbourhood sum, the exact good outcome, and a heatmap over the offsets $abs(h) <= k$.],
    [`glitchy.py worst`], [Searches every placement of $k$ glitches, or a random sample, for any $k$ and sign pattern. Lists the placements with the smallest $min_c S_c (k)$.],
    [`glitchy.py check`], [Reproduces the numerical results and checks the new total-mass bounds in a few seconds.],
    [`Experiment_plots/Log_Plot{2,3,4,5}.txt`], [2019 logs: placement (`ListK`), shifts, and the 20 worst placements. The plot images were too large to commit.],
    table.hline(),
  ),
  caption: [Code in the repository. `glitchy.py` (2026) replaces the 2019 scripts `Plot.py`,
    `Find_Worst_Case.py` and `Test.py`, which remain in git history (commit `a9e171b`). It
    computes the same distribution: on test configurations, including $r divides.not N$, it agrees
    with the 2019 formula to $10^(-16)$. The `Test.py` experiment (adding one extra glitch) had
    no logged output and was not ported. Good outcomes are the integers nearest $c N\/r$, the same
    set as the 2019 rule $(y r mod N) <= floor(r\/2)$ or $>= N - floor(r\/2)$. They equal $c M$
    when $r divides N$.],
) <tab:code>

== Worst cases

#figure(
  table(
    columns: 10,
    stroke: none,
    align: (left, right, right, right, right, left, left, right, right, right),
    table.hline(),
    table.header[*Log*][$N$][$r$][$M$][$k$][*signs*][*placements*][$S^* (k)$][*at* $c\/r$][*continuum*],
    table.hline(stroke: 0.5pt),
    [—], [3000], [50], [60], [1], [ins.], [all 59], [0.6403], [0.50], [0.6400],
    [Plot2], [3000], [50], [60], [2], [ins.], [all 1711], [0.7723], [0.70], [0.7720],
    [Plot3], [3000], [100], [30], [3], [del.], [all 3654], [0.6136], [0.50], [0.6119],
    [Plot4], [1500], [50], [30], [4], [ins.], [all 23751#super[†]], [0.6888], [0.40], [0.6858],
    [Plot5], [3000], [50], [60], [5], [del.], [5063 (grid)], [not logged], [—], [—],
    table.hline(),
  ),
  caption: [Worst neighbourhood ratio $S^*(k) = r dot min_c sum_(abs(h) <= k) P_k (c M + h)$ in
    Model I, over all placements with the given signs, and the good outcome $c\/r$ where it is
    attained. The worst placements (block indices $m_1, ..., m_k$) are $(6, 54)$ for Plot2,
    $(3, 6, 27)$ for Plot3 and $(2, 4, 6, 8)$ for Plot4, whose mirror image $(22, 24, 26, 28)$
    ties. The last column evaluates the same configuration in the continuum limit.
    #super[†]The 2019 log searched only 1820 grid placements and found $0.7229$ at
    $(3, 24, 25, 28)$. The full search, the $k = 1$ row and the $c\/r$ column were computed for this
    note. Plot3 and Plot4 have $r > sqrt(N)$, which does not matter when $r divides N$ but violates
    the standing assumption. For Plot2 and Plot3 the full search finds the same minima and
    minimizers as the 2019 logs.],
) <tab:worst>

All worst cases stay well above $1\/2$. The worst good outcome is $c = r\/2$ ($alpha = 1\/2$) for
$k = 1, 3$, but $alpha = 0.7$ (equivalently $0.3$, since $S$ is even in $alpha$) for $k = 2$ and
$alpha = 0.4$ for $k = 4$. So $alpha = 1\/2$ is not always the worst case. The $k = 4$ minimizer
$(2, 4, 6, 8)$ packs the glitches at one end, like the family in @sec:packed. The values are not
monotone in $k$: $0.640, 0.772, 0.614, 0.689$ for $k = 1, 2, 3, 4$, though at different $N$ and $r$.
Any induction on $k$ has to account for the radius growing with $k$.

== Sign patterns <sec:signs>

At $alpha = 1\/2$ the phase function is $phi_j = e^(i pi s_j) = (-1)^(j-1)$ whatever the signs,
because a shift of $plus.minus 1$ changes the phase by $e^(plus.minus i pi) = -1$. So in the continuum limit
$S_(1\/2)$ does not depend on the sign pattern; signs matter only for $alpha != 1\/2$. Searching
all placements with $N = 3000$, $r = 50$ #pill([numerical evidence only], rgb("#5b3f99")) gives:

#align(center, table(
  columns: 7,
  stroke: none,
  align: (left, right, right, right, right, right, right),
  table.hline(),
  [signs], [$+ +$], [$+ -$], [$+ + +$], [$+ + -$], [$+ - +$], [$+ - -$],
  table.hline(stroke: 0.5pt),
  [$S^* (k)$], [0.7723], [0.7830], [0.6083], [0.6083], [0.6083], [0.6083],
  [at $c\/r$], [0.70], [0.50], [0.50], [0.50], [0.50], [0.50],
  table.hline(),
))

Flipping every sign gives the same values, because $S$ is even in $alpha$. In these runs mixed
signs are never worse than equal signs. For $k = 3$ the worst case sits at $alpha = 1\/2$, where
signs cannot matter.

== Figures from the poster

#figure(
  grid(
    columns: (62%, 30%),
    gutter: 10pt,
    align: horizon,
    image("poster_figures/two_glitch_probability_distribution.png", width: 100%),
    image("poster_figures/probability_distribution_legend.png", width: 100%),
  ),
  caption: [Poster figure 3: $N = 3000$, $r = 20$, $k = 2$. Each column is one good outcome $c M$.
    Blue is $P_k (c M)$. Green and yellow are $P_k (c M minus.plus 1)$, and cyan and magenta are
    $P_k (c M minus.plus 2)$. Black is the sum over $abs(h) <= 2$. The dotted line is the ideal $1\/r = 0.05$.
    Where the good outcome loses probability, its neighbours gain it.],
) <fig:p3>

#figure(
  image("poster_figures/ten_glitch_neighborhood_probability.png", width: 72%),
  caption: [Poster figure 4: $N = 11000$, $r = 100$, $k = 10$, one random placement. The black
    neighbourhood sum ($abs(h) <= 10$) stays above about $0.0077$, which is $0.77$ of the ideal
    $0.01$. The minima are near, but not at, $c = r\/2$. The placement was not logged.],
) <fig:p4>

#figure(
  image("poster_figures/nondivisible_period_probability_distributions.png", width: 100%),
  caption: [Poster figure 5: $N = 11000$, $r = 107$ ($r divides.not N$) with no glitch (left), $2$
    glitches (middle) and $3$ glitches (right). Without glitches the nearest-integer good outcomes
    already carry between roughly $4\/pi^2$ and $1$ times $1\/r$. No analysis exists for this case.],
) <fig:p5>

= Further observations made while compiling this note <sec:new>

== A packed family approaching 1/2 <sec:packed>

Place all $k$ glitches, with equal signs, equally spaced in $[0, delta)$, and take $alpha = 1\/2$
(even $r$, $c = r\/2$). Minimizing $S^"cont"_(1\/2)(k)$ over a grid of $delta$ values
#pill([numerical evidence only], rgb("#5b3f99")) gives:

#align(center, table(
  columns: 9,
  stroke: none,
  align: right,
  table.hline(),
  [$k$], [2], [4], [6], [10], [20], [50], [100], [200],
  table.hline(stroke: 0.5pt),
  [$min_delta S$], [0.782], [0.696], [0.653], [0.611], [0.571], [0.542], [0.531], [0.521],
  [$delta$], [0.245], [0.330], [0.375], [0.415], [0.455], [0.480], [0.490], [0.490],
  table.hline(),
))

Heuristically, on $[0, delta)$ the phase alternates $plus.minus 1$ every $delta \/ k$. This is a
square wave of frequency $k \/ (2 delta)$, which exceeds $k$ when $delta < 1\/2$. Its energy, about
$delta$, therefore lies outside the window $abs(h) <= k$, while the constant part on $[delta, 1)$
keeps about $1 - delta$. Letting $delta arrow.t 1\/2$ as $k -> infinity$ suggests
$inf S^"cont" (k) -> 1\/2$. *Consequence:* if @conj:half is true, its constant is asymptotically
sharp, and a proof must be tight on this family. Neither the limit nor the existence of
configurations below $1\/2$ has been settled.

== The radius must grow linearly in $k$

#proposition(name: [Radius lower bound], tags: (new, proved))[
  Let $k$ be odd, the glitches equally spaced ($u_j = j\/(k+1)$), and the signs equal.
  + If $alpha = 1\/2$, then $a_h = 0$ unless $h$ is an odd multiple of $(k+1)\/2$. Hence
    $S_(1\/2)^"cont" (L) = 0$ for $L < (k+1)\/2$, and $S_(1\/2)^"cont" (L) = 8\/pi^2$ for
    $(k+1)\/2 <= L < 3(k+1)\/2$.
  + If $r$ is odd and $c = (r-1)\/2$, which is informative because $gcd(c, r) = 1$, then
    $S_(c\/r)^"cont" (L) <= (2L + 1) (pi k \/ r)^2$ for $L < (k+1)\/2$.
] <prop:radius>
#proof[
  (1) The $k + 1$ pieces alternate $plus.minus 1$, and since $k + 1$ is even this holds across the
  wrap-around too. So $phi(u + tau) = -phi(u)$ with $tau = 1\/(k+1)$, and therefore
  $a_h = -e^(2 pi i h tau) a_h$. Hence $a_h = 0$ unless $e^(2 pi i h tau) = -1$. The energy at the
  fundamental frequency of a $plus.minus 1$ square wave is $8\/pi^2$.
  (2) $abs(a_h (alpha) - a_h (1\/2)) <= norm(phi_alpha - phi_(1\/2))_1 <= max_j abs(e^(2 pi i (alpha - 1\/2) s_j) - 1)
  <= 2 pi abs(alpha - 1\/2) k = pi k \/ r$. Then apply (1).
]

So any bound that is uniform in $c$, or that holds near the informative $c approx r\/2$, needs
radius $L >= (k+1)\/2$ (when $r >> k^(3\/2)$). The same configuration also suggests a trade-off.
For a general $alpha$ its energy concentrates near $h equiv -alpha (k+1) thick (mod k+1)$, so a radius $L$
captures roughly a fraction $2L\/(k+1)$ of the values of $c$. The total cost, (samples) $times$
(candidates per sample), is then about $k$ for every choice of $L$. This is heuristic.

== From the literal model to Model I

#proposition(name: [Transfer], tags: (new, sketch))[
  Let $f$ be a literal $k$-glitch function with $k < r$ (@sec:model, remark on idealizations), and
  round each glitch to the next period boundary to obtain a Model I configuration. Then for all
  $c$ and $L$,
  $ abs(S_c^"literal" (L) - S_c^"Model I" (L)) = O(L k (1 + k\/M) \/ M) . $
  Together with @lem:continuum, a continuum bound $S^"cont" (k) >= s_0$ gives
  $S^"literal" (k) >= s_0 - O(k^2 (1 + k\/r)\/M)$. Since $M > sqrt(N)$, this is meaningful for
  $k = o(N^(1\/4))$.
] <prop:transfer>
#proof(title: [Sketch])[
  For each class $a$ the two position sets differ in $delta_a$ points (moved, added or removed),
  with $delta_a = O(k)$ and $sum_a delta_a = O(k r)$. A glitch inside a period moves at most $r - 1$
  points, one per class. Inserted symbols and truncation contribute $O(k)$ more. New symbols form
  singleton classes, which add at most $k\/N^2$ to each $P(y)$. Since each class amplitude
  changes by at most $2 delta_a$ and has size at most $M + delta_a$,
  $abs(P(y) - P'(y)) <= N^(-2) sum_a 2 delta_a (2M + 2 delta_a) = O(k (M + k) r \/ N^2)$.
  Multiplying by $r(2L+1)$ gives the claim. The bookkeeping of $delta_a$, in particular several
  glitches in one period and mixed signs, still has to be written out.
]

= From neighbourhood mass to period recovery <sec:recovery>

The poster's hypothesis is that the period can be learned by examining every point within a
certain range of the measured outcome. Made explicit, the procedure is
#pill([proposal, not proved], rgb("#b35900")):

+ Run the period-finding circuit on $f_k$ and measure $y$.
+ For each $y'$ with $abs(y' - y) <= L$, expand $y'\/N$ in continued fractions and collect the
  convergent denominators $q <= sqrt(N)$ as candidates.
+ For each candidate, sample random $x$ and test $f_k (x) = f_k (x + q)$. Output the smallest $q$
  whose pass rate exceeds a threshold $1 - tau$.

*Why it should work* ($r divides N$, Model I):
- If $y$ lands within $L$ of $c M$ with $gcd(c, r) = 1$, then $y' = c M$ is in the window and
  $y'\/N = c\/r$ exactly, so $q = r$ is a candidate.
- For $q = r$, the test fails only when a glitch lies in $(x, x + r]$, which happens for an
  $O(k\/M)$ fraction of $x$.
- For $q equiv.not 0 mod r$, passing requires a glitch in $(x, x+q]$, which happens for at most a
  $k q \/ N <= k\/sqrt(N)$ fraction of $x$. Multiples of $r$ also pass, so the smallest passing
  candidate is $r$.
- *Cost.* Each sample needs $2L + 1$ candidates at $"poly"(log N)$ classical cost each. If
  $S_c (L) >= s_0$ for every informative $c$, the expected number of samples is
  $O(r \/ (phi(r) s_0)) = O(log log r \/ s_0)$. With $L = k$ and @conj:half this is the poster's
  "cost scales like the number of glitches".

*What is missing:* a lower bound on the mass near _informative_ $c$, which a bound uniform in $c$
would supply; a careful analysis of the classical test, including inserted symbols; and the
$r divides.not N$ case, where candidates must satisfy $abs(y'\/N - c\/r) <= 1\/(2 r^2)$.
The general-$k$ averaged guarantees in @thm:total and @prop:model-total do not supply the
informative-outcome bound.

= Open problems <sec:open>

+ *Main conjecture.* Prove @conj:half: $S^"cont" (k) >= 1\/2$. Failing that, prove any constant
  $c_0 > 0$ at radius $O(k)$, or find a configuration below $1\/2$. Known so far: it holds for
  $k = 1$ (@thm:one, sharp value $0.640$); exhaustive Model I search confirms it for $k = 2, 3, 4$
  (@tab:worst) and for every sign pattern with $k <= 3$ at $N = 3000$, $r = 50$ (@sec:signs); it holds at radius $Theta(k^2)$ (@cor:quadratic); and the constant cannot be
  improved asymptotically (@sec:packed). The step-function form @eq:jump is the natural setting.
+ *Extremal configurations.* Which $alpha$ is worst? It is $1\/2$ for $k = 1, 3$, but $0.3$ for
  $k = 2$ and $0.4$ for $k = 4$ (@tab:worst). Do worst cases cluster their glitches, as the packed
  family and the minimizers $(6, 54)$, $(3, 6, 27)$ and $(2, 4, 6, 8)$ suggest? Numerical optimization over
  $(u_1, ..., u_k)$ in the continuum model for $k <= 10$ is cheap.
+ *Mixed insertions and deletions.* Only the small runs of @sec:signs exist, and there mixed
  signs are never worse. Prove that equal signs are always the worst case, or find a
  counterexample. With mixed signs the phases $phi_j = e^(2 pi i alpha s_j)$ follow a
  $plus.minus 1$ walk instead of a geometric progression.
+ *Rigorous transfer.* Complete @prop:transfer with explicit constants, including several
  glitches in one period, and determine the admissible range of $k$ relative to $N$ and $r$.
+ *Period recovery.* Turn @sec:recovery into a theorem. Determine whether averaged bounds such
  as @thm:total can be strengthened to control informative outcomes, and whether the
  radius–samples trade-off of @prop:radius is real.
+ *$r divides.not N$, per outcome.* @thm:nondividing proves constant total mass on the
  radius-$k$ union of rounded peaks. Extend the per-outcome formula @lem:exact to non-integer
  peak positions $c N\/r$, as in @fig:p5, and control the informative outcomes.
+ *Random glitches.* An average-case model, with uniformly random positions, should be much
  easier and may give better constants. It is a different theorem from the worst-case one.
+ *Literature.* As of October 2026 the author found no public paper, preprint or manuscript on
  glitchy period finding beyond the 2019 poster. In particular, nothing is available for the
  QIP 2026 follow-up mentioned in the research prompt. Do not rely on that follow-up as a source.
  It remains worth searching for related work under other names, such as robustness of period
  finding or of the hidden subgroup problem to oracle errors, and Fourier sampling with corrupted
  oracles.
+ *Other algorithms.* The poster's future work: robustness of other Fourier-sampling or
  hidden-subgroup algorithms to glitchy oracles.

#counter(heading).update(0)
#set heading(numbering: "A.1", supplement: [Appendix])
= Reproducibility <app:repro>

The numerical simulator is `glitchy.py`, at the repository root. It needs `numpy`, and
`matplotlib` for plots. Lean proofs and diagram sources are separate, as described below. With #link("https://docs.astral.sh/uv/")[uv] the dependencies install automatically:

```
uv run glitchy.py check                                  # every number in this note
uv run glitchy.py plot  --N 11000 --r 100 --k 10         # random placement
uv run glitchy.py plot  --N 3000 --r 50 --blocks 6 54    # given placement
uv run glitchy.py worst --N 1500 --r 50 --k 4            # all placements
uv run glitchy.py worst --N 3000 --r 50 --k 3 --signs 1 -1 1
```

`check` runs a deterministic numerical test suite. It prints (1) the Model I worst cases of @tab:worst, recomputed
from the logged placements and by full search; (2) the sign-pattern table of @sec:signs; (3) the
one-glitch constants of @thm:one, with the closed form checked against direct summation; (4) the
packed-family table of @sec:packed; (5) the radius example of @prop:radius; and (6) checks of the
residue-count identity @eq:residue-mass and both bounds in @thm:total against direct Fourier
measurement probabilities. Literal tests enumerate small two- and three-edit cases and use a
fixed seed for larger $k$, including $k >= r$, $k > N$, mixed signs, adjacent edits, repeated
or new inserted values, and truncation. They also verify the worked example in
@fig:residue-example and equality cases of @prop:model-total. (7) The nondivisible-period tests
check @eq:matching-moment, the edit-pair bound @eq:matching-bound, rounding separation, the
Fejér-kernel expectation and both cases of @thm:nondividing against the exact literal-prefix
Fourier distribution. They cover half-integer peaks, mixed signs, old/new inserted symbols,
truncation and overlapping windows. The theorem's proof, not these finite checks, supplies the
uniform guarantee.

The poster figures were cropped from #raw("Summer_Symposium_Poster.pdf"). The explanatory
figures in @sec:total have editable sources in #raw("report_figures/"). Rebuild them with:

```sh
mmdc -i report_figures/general_k_proof.md -o report_figures/general_k_proof.svg
python3 report_figures/build_general_k_example.py
typst compile glitchy_period_finding_report.typ
```

The Mermaid command emits #raw("general_k_proof-1.svg") from the single fenced diagram. The
committed SVGs are sufficient to compile the report; Mermaid is needed only to regenerate its
proof map.

*Lean verification.* With #link("https://github.com/leanprover/elan")[elan] installed:

```sh
cd formal
lake exe cache get       # first use: download precompiled mathlib dependencies
lake build
lake env lean Audit.lean
```

The toolchain and dependencies are pinned in #raw("formal/lean-toolchain"),
#raw("formal/lakefile.toml"), and #raw("formal/lake-manifest.json"). The axiom audit covers the
literal-edit theorem, both sharper bounds, their neighbourhood and finite-size corollaries,
the Model I bound and sharpness, the worked example, and the complete nondivisible-period
matching/Fejér argument in @thm:nondividing. If macOS rejects the cache executable,
#raw("formal/README.md") documents an interpreter-based cache-fetch command. The verification
itself uses the same Lean kernel in either case.
