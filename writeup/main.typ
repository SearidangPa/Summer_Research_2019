// Quantum Glitchy Period Finding: status report and open problems.
// Compile with:  typst compile writeup/main.typ
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
  would be asymptotically sharp. A list of open problems closes the note.
]

#v(4pt)
#block(inset: (x: 18pt))[
  #set text(size: 9pt)
  *How to read the status labels.*
  #from2019 the statement or formula appears in the 2019 poster, code or logs.
  #claim2019 stated on the 2019 poster; no proof is in the repository.
  #new added while compiling this note.
  #proved a complete argument is given here.
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
- *Not done.* No multi-glitch theorem at radius $O(k)$. No period-recovery theorem
  (@sec:recovery). No analysis of $r divides.not N$. No public follow-up paper or manuscript exists as of October 2026
  (@sec:open, item 8).

= Standard period finding <sec:background>

Let $f : ZZ_N -> W$ have period $r$, with $r < sqrt(N)$. For most of this note we also assume
#footnote[The poster's figure 5 treats $r divides.not N$; nothing is proved there.]
$ r divides N, quad M := N\/r, quad f(x) = w_(x mod r) "with" w_0, ..., w_(r-1) "distinct". $ <eq:standing>
The algorithm prepares $N^(-1\/2) sum_x lr(|x chevron.r) lr(|f(x) chevron.r)$, applies the quantum Fourier
transform to the first register and measures. Outcome $y$ appears with probability
$ P_f (y) = 1/N^2 sum_(w in W) abs(sum_(x : f(x) = w) e^(2 pi i x y \/ N))^2 . $ <eq:pf>
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
    image("figures/poster_fig1_insertion.png", width: 100%),
    image("figures/poster_fig2_deletion.png", width: 100%),
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
equals $(1\/r) sum_c S_c (k) >= S^* (k)$ when $k < M\/2$. Every statement below is per $c$ and
therefore also bounds $S_k$.

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
weaker than @conj:half.

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
    [`glitchy.py check`], [Reproduces every number in this note in about two seconds.],
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
    image("figures/poster_fig3.png", width: 100%),
    image("figures/poster_legend.png", width: 100%),
  ),
  caption: [Poster figure 3: $N = 3000$, $r = 20$, $k = 2$. Each column is one good outcome $c M$.
    Blue is $P_k (c M)$. Green and yellow are $P_k (c M minus.plus 1)$, and cyan and magenta are
    $P_k (c M minus.plus 2)$. Black is the sum over $abs(h) <= 2$. The dotted line is the ideal $1\/r = 0.05$.
    Where the good outcome loses probability, its neighbours gain it.],
) <fig:p3>

#figure(
  image("figures/poster_fig4.png", width: 72%),
  caption: [Poster figure 4: $N = 11000$, $r = 100$, $k = 10$, one random placement. The black
    neighbourhood sum ($abs(h) <= 10$) stays above about $0.0077$, which is $0.77$ of the ideal
    $0.01$. The minima are near, but not at, $c = r\/2$. The placement was not logged.],
) <fig:p4>

#figure(
  image("figures/poster_fig5.png", width: 100%),
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
+ *Period recovery.* Turn @sec:recovery into a theorem. Determine whether bounds averaged over
  $c$ suffice and whether the radius–samples trade-off of @prop:radius is real.
+ *$r divides.not N$.* Extend @lem:exact to peaks at non-integer positions $c N\/r$, as in @fig:p5.
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

All code is in one file, `glitchy.py`, at the repository root. It needs `numpy`, and `matplotlib`
for plots. With #link("https://docs.astral.sh/uv/")[uv] the dependencies install automatically:

```
uv run glitchy.py check                                  # every number in this note
uv run glitchy.py plot  --N 11000 --r 100 --k 10         # random placement
uv run glitchy.py plot  --N 3000 --r 50 --blocks 6 54    # given placement
uv run glitchy.py worst --N 1500 --r 50 --k 4            # all placements
uv run glitchy.py worst --N 3000 --r 50 --k 3 --signs 1 -1 1
```

`check` runs in about two seconds. It prints (1) the Model I worst cases of @tab:worst, recomputed
from the logged placements and by full search; (2) the sign-pattern table of @sec:signs; (3) the
one-glitch constants of @thm:one, with the closed form checked against direct summation; (4) the
packed-family table of @sec:packed; and (5) the radius example of @prop:radius. The poster
figures were cropped from `Summer_Symposium_Poster.pdf`.
