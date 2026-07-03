# ML Parameter Estimation in Clann: β, η, and the Normalising Constant

## Overview

Clann implements the maximum-likelihood supertree criterion of Steel & Rodrigo
(2008). This note explains the statistical model, derives the closed-form
estimator for the rate parameter β, describes the experimental tree-size
scaling exponent η (formerly called α in earlier versions — renamed to avoid
confusion with the normalising constant α used in Steel & Rodrigo 2008 and
Akanni *et al.* 2014), and discusses the normalising constant Z_T whose
omission in the original paper was corrected by Bryant & Steel (2008).

---

## 1. The Steel & Rodrigo (2008) model

The model treats each source tree T_i as an independent observation drawn from
an exponential distribution parameterised by the RF distance to the supertree S:

    P(T_i | S, β) = (1 / Z_T) · exp(−β · d(T_i, S|X_i))

where d(T_i, S|X_i) is the Robinson-Foulds distance between T_i and the
restriction of S to T_i's taxa, β > 0 is a global rate parameter, and Z_T is
the normalising constant (see Section 4 below).

Large β implies high confidence in the source trees; small β implies high
tolerance for disagreement.

### 1.1 The normalising constant

Z_T is the sum of exp(−β·d) over *all* fully-resolved trees on the same taxa
as T_i:

    Z_T = Σ_{T'} exp(−β · d(T', S|X_i))

The original Steel & Rodrigo (2008) paper omitted Z_T, claiming it was
approximately constant across candidate supertrees. Bryant & Steel (2008)
showed this is only approximately true — Z_T varies with the shape of S|X_i
(in particular its cherry count) — and that the omission can affect topology
rankings when β is near 1.5 and source trees have many taxa. See Section 4
for details.

When Z_T is ignored, the joint log-likelihood over all source trees, allowing
per-tree weights w_i, simplifies to:

    log L(S, β) = W·log(β) + Σ w_i·log(w_i) − β · Σ_i w_i · d(T_i, S|X_i)
                = W·log(β) + H_w − β · WD

where W = Σ_i w_i, WD = Σ_i w_i · d(T_i, S|X_i) is the weighted RF sum,
and H_w = Σ w_i·log(w_i) is a constant determined by the weights. (The
Σ w_i·log(w_i) term arises because per-tree β_i = β·w_i when weights are used
as β multipliers; it is zero when all weights equal 1.)

### 1.2 Closed-form MLE for β

Setting ∂ log L / ∂β = 0:

    W/β − WD = 0
    β̂ = W / WD

This is the unique global maximum (the log-likelihood is strictly concave in β
for WD > 0). No numerical optimisation is required. The maximum log-likelihood
is:

    log L(β̂) = W·log(W/WD) − W + H_w

### 1.3 Using mlscores

The `mlscores` command computes β̂ and log L(β̂) for the current supertree and
source trees, then updates `ml_beta` so that all subsequent `hs` and `boot`
runs use the estimated value. Optionally, a log-likelihood profile over a
log-uniform grid of β values can be written to a file (`outfile=` option) for
visual inspection of the likelihood surface.

    mlscores
    mlscores outfile=beta_profile.txt scan=200

---

## 2. The tree-size problem

The Steel model gives every split disagreement equal cost regardless of which
source tree it comes from. A source tree with k_i = 50 internal splits can
contribute up to 50 + |S| disagreements; one with k_i = 5 can contribute at
most 5 + |S|. In a heterogeneous dataset — typical in phylogenomics, where gene
trees range from 4-taxon quartets to near-complete trees — large trees
disproportionately dominate WD, and therefore dominate both β̂ and the search
objective.

This is not necessarily wrong: a tree with more bipartitions carries more
phylogenetic information per se. However, it conflates two distinct quantities —
informational content (which may justify higher weight) and mere size (which
should not). If large trees are systematically less reliable (more incomplete
lineage sorting, more estimation error, lower bootstrap support), the raw RF
model will be biased in their favour.

The `mleta` parameter provides a principled way to investigate and correct
this bias.

---

## 3. The tree-size scaling exponent η

### 3.1 Derivation from the probability model

Rather than bolting a scaling factor onto the distance after the fact, η is
derived by modifying the generative model itself. Suppose the rate of the
exponential for tree i is β / k_i^η rather than β:

    P(d_i | S, β, η) = (β / k_i^η) · exp(−(β / k_i^η) · d_i)

This is still a valid exponential distribution for each tree; the rate is now
tree-specific and depends on size. The joint log-likelihood is:

    log L(β, η) = n·log(β) − η · Σ_i log(k_i) − β · Σ_i w_i · d_i / k_i^η
                = n·log(β) − η · P − β · WD(η)

where P = Σ_i log(k_i) is a fixed penalty constant (determined by the data)
and WD(η) = Σ_i w_i · d_i · k_i^{−η}.

### 3.2 Why η does not degenerate

The term −η · P is a natural penalty that grows linearly with η. It arises
directly from the normalisation constant of the per-tree exponential: as the
rate β/k_i^η is made smaller (larger k_i^η), the distribution spreads out and
the probability of observing any particular d_i decreases. This penalty
counteracts the gain from reducing WD(η), creating an interior optimum
for η.

This is the critical difference from naive weighting schemes. If one simply
multiplies d_i by an ad hoc factor that decreases with k_i, and then
re-estimates β, the likelihood will always increase as the factor approaches
zero because WD decreases without any corresponding penalty term. The η
parameter avoids this because it is derived from a fully specified probability
model.

### 3.3 Profiling β at each η

For any fixed η, the MLE for β is still closed-form:

    β̂(η) = W / WD(η)

Substituting back:

    log L(η) = W·log(W / WD(η)) − W − η · P

This reduces the joint optimisation over (β, η) to a 1-D search over η alone.
At each grid point η_k, WD(η_k) is computed from the raw RF distances and
β̂(η_k) is obtained immediately. The η with the highest log L(η) is the MLE.

### 3.4 Interpretation of η*

| η value | Interpretation |
|---------|----------------|
| 0 | Original Steel & Rodrigo (2008) model — raw RF, no size correction |
| (0, 1) | Partial correction — large trees still contribute more, but less than raw |
| 1 | Full normalisation — each tree contributes equally per unit weight, regardless of size |
| > 1 | Active down-weighting of large trees — appropriate if large trees are less reliable |
| → ∞ | All weight collapses onto the smallest trees (resisted by the −η·P penalty) |

The MLE condition ∂ log L / ∂η = 0 can be written as:

    n · 〈log k〉_{w(η)} = P = Σ_i log(k_i)

where 〈log k〉_{w(η)} is the weighted average of log(k_i) with weights
proportional to w_i · d_i · k_i^{−η}. This is a moment-matching condition: the
model fits when the weighted average log-size (weighted by disagreement) equals
the unweighted sum of log-sizes. Whether this has an interior solution depends
on the (d_i, k_i) distribution in the dataset.

If large trees tend to disagree more in absolute terms but less in relative
terms (d_i / k_i smaller for larger k_i), there is a meaningful η* > 0. If
relative disagreement is uniform across sizes, the penalty and gain cancel and
η* will be near 0.

### 3.5 Confidence intervals and identifiability

The curvature of log L(η) near η* determines the precision of the estimate. A
flat profile (low curvature) indicates that the data are not informative about
η — tree-size effects are absent or confounded. A sharp peak indicates a
well-identified size effect. The profile written by `mlscores outfile=` provides
the full likelihood surface for visual assessment.

---

## 4. The normalising constant Z_T (Bryant & Steel 2008)

### 4.1 Why Z_T matters

The full model probability for source tree T_i is:

    P(T_i | S) = exp(−β · d_i) / Z_T

where Z_T = Σ_{T'} exp(−β · d(T', S|X_i)) sums over all fully-resolved trees
on T_i's taxa. The standard ML search maximises Σ w_i · log P(T_i|S), which
means the term −Σ w_i · log(Z_{T|X_i}) should appear in the lnL. Omitting it
introduces two types of error:

1. **Absolute lnL bias**: every source tree's contribution is wrong by
   log(Z_{T|X_i}), typically ~log(2–3) nats per tree at β ≈ 1.9, n ≈ 33.
2. **Topology ranking bias**: different candidate supertrees have different
   Z_{T|X_i} values because they produce differently-shaped restrictions T|X_i.
   If the variation |log Z_{T1} − log Z_{T2}| ≥ 2β for any pair of candidates,
   the topology ranking can change.

### 4.2 Magnitude of the effect

Bryant & Steel (2008) showed that Z_T depends only weakly on tree shape in most
regimes. The key results are:

- **Small β (β < 0.03):** Z_T ≈ b(n) = (2n−5)!! (total number of binary
  trees). Astronomically large; varies slightly with cherry count.
- **Large β (β > 3):** Z_T → 1. The correction becomes negligible.
- **Middle range β ∈ [0.03, 3]:** Z_T changes rapidly across many orders of
  magnitude; topology-dependent variation peaks near β ≈ 1.5.

The condition for topology rankings to change is |log Z_{T1} − log Z_{T2}| ≥ 2β,
where the difference is driven by the cherry-count difference Δc_T between the
two candidate topologies' restrictions. Bryant & Steel showed:

- **n ≤ 20 source trees:** No value of β causes a ranking change (safe to ignore Z_T).
- **n = 50, β ∈ [1.25, 1.86]:** Rankings can potentially change.
- **n = 33, β ≈ 1.9:** Variation ≈ 0.045 ≪ 2β = 3.8; rankings unaffected.

### 4.3 The large-β approximation

Clann implements a truncated large-β expansion for log Z_T:

    Z_T ≈ 1 + b₂ · ε + b₄(c_T) · ε²

where:
- ε = e^{−2β}
- b₂ = 2(n−3): number of trees at RF distance 2 (same for all binary n-leaf trees)
- b₄(c_T) = 4·C(n−3, 2) + 6·(n−6+c_T): depends on cherry count c_T of T|X_i
- c_T: number of cherry pairs (splits with exactly 2 leaves on one side)

Therefore:

    log Z_T ≈ log(1 + 2(n−3)·e^{−2β} + [4·C(n−3,2) + 6·(n−6+c_T)]·e^{−4β})

This is accurate when β > ~1.5 (ε < 0.05). For β < 1.5, the series converges
slowly and the exact O(n⁵) algorithm of Bryant & Steel (2008) would be needed
for accuracy, though practical topology impact remains small for n < 30.

The cherry count c_T is computed by parsing the pruned Newick string T|X_i
directly: it counts internal nodes whose subtree contains exactly 2 taxa.

### 4.4 The `normcorrect` option

The correction is applied in `usertrees` when the `normcorrect` flag is set:

    usertrees candidates.ph normcorrect
    usertrees candidates.ph tests=yes normcorrect

For each candidate supertree T and each source tree T_i, Clann:
1. Builds the restriction T|X_i (already done for RF scoring)
2. Counts cherries c_T in T|X_i from the pruned Newick
3. Computes log Z_{T_i} ≈ log(1 + b₂ε + b₄ε²)
4. Subtracts w_i · log Z_{T_i} from source tree i's lnL contribution

The corrected total lnL is:

    lnL_corrected = W·log(β) + Σ w_i·log(w_i) − β·WD − η·P − Σ w_i·log(Z_{T_i})

**Important caveats:**
- The correction is only applied during `usertrees`, not during `hs`/`boot`
  (the Z_T correction is too expensive to apply at every heuristic search step,
  and has negligible effect on topology rankings at typical β values).
- KH and SH test statistics (which use differences in lnL between candidate
  trees) are affected only by the topology-dependent part of Z_T — the
  cherry-count term. For typical datasets the change is very small.
- `normcorrect` is only active with `criterion=ml` and `mlscale=lnl`.

### 4.5 Naming note: η vs α

Both Steel & Rodrigo (2008) and L.U.St (Akanni *et al.* 2014) use α to denote
the normalising constant a_i = 1/Z_T (or equivalently, α ≈ 1 when Z_T ≈ 1 for
large β). Clann formerly used `ml_alpha` for the tree-size scaling exponent
described in Section 3 — a completely different quantity. To avoid this naming
conflict, Clann's tree-size exponent was renamed η (`ml_eta`, option `mleta`)
in Clann v5.

---

## 5. Using η estimation in practice

### 5.1 Estimating η

    set criterion ml
    exe gene_trees.ph
    hs
    mlscores eta=auto escan=100 etamax=3.0 outfile=profile.txt

This runs the standard β MLE first, then performs a 100-point grid search over
η ∈ [0, 3]. Both `ml_eta` and `ml_beta` are updated to the joint MLE values.
The output file contains the β profile followed by the η profile.

### 5.2 Setting η manually

If you have a prior expectation or want to fix η for comparison:

    set mleta 1.0    # full normalisation
    hs               # search uses eta=1 throughout
    mlscores         # estimate beta conditional on eta=1

### 5.3 Interpreting the result

- If η* ≈ 0: tree-size heterogeneity is not biasing your result; the original
  Steel model is appropriate.
- If η* ≈ 1: large trees dominate purely because of their size; full
  normalisation is warranted.
- If η* is intermediate: a partial size correction is supported by the data.
- If η* > 1: the data favour actively down-weighting large trees, suggesting
  that size correlates with unreliability in your dataset.
- If the optimal supertree topology changes substantially between η = 0 and η*,
  tree-size heterogeneity is a confounding factor and should be reported.

---

## 6. Relationship to tree weights

The `mleta` mechanism is conceptually distinct from per-tree weights (w_i).
Weights encode prior confidence in individual trees (e.g. from bootstrap
support via `autoweight=bs`). η encodes a systematic relationship between tree
size and reliability that applies uniformly across all trees. The two can be
used together: weights handle tree-specific variation; η handles the systematic
size trend.

---

## 7. Summary of commands

| Command / option | Effect |
|---|---|
| `mlscores` | Estimate β̂ = W/WD for current supertree; update `ml_beta` |
| `mlscores outfile=<f>` | Also write β log-likelihood profile to file |
| `mlscores eta=auto` | Grid search over η ∈ [0, etamax]; update `ml_eta` and `ml_beta` |
| `mlscores escan=<n>` | Number of grid points for η search (default 50) |
| `mlscores etamax=<f>` | Upper bound of η grid (default 3.0) |
| `mlscores fixbeta` | Hold β fixed at current `ml_beta`; do not recompute |
| `mlscores sourcescores=<f>` | Write per-source-tree lnL contributions to file |
| `set mleta <f>` | Fix η to a specific value for subsequent hs/boot runs |
| `set mlbeta <f>` | Fix β to a specific value manually |
| `usertrees <f> normcorrect` | Apply Bryant & Steel (2008) Z_T correction to lnL in usertrees |

---

## 8. References

Bryant, D. and Steel, M., 2009. Computing the distribution of a tree metric. IEEE/ACM transactions on computational biology and bioinformatics, 6(3), pp.420-426.

Steel, M. and Rodrigo, A., 2008. Maximum likelihood supertrees. Systematic biology, 57(2), pp.243-250.

Akanni, W.A., Creevey, C.J., Wilkinson, M. and Pisani, D., 2014. LU St: a tool for approximated maximum likelihood supertree reconstruction. BMC bioinformatics, 15(1), p.183.
