# ML Parameter Estimation in Clann: β and α

## Overview

Clann implements the maximum-likelihood supertree criterion of Steel & Rodrigo
(2008). This note explains the statistical model, derives the closed-form
estimator for the rate parameter β, and describes the experimental tree-size
scaling exponent α — including its theoretical justification and how it differs
from naive weighting schemes that degenerate.

---

## 1. The Steel & Rodrigo (2008) model

The model treats each source tree T_i as an independent observation drawn from
an exponential distribution parameterised by the RF distance to the supertree S:

    P(T_i | S, β) = β · exp(−β · d(T_i, S))

where d(T_i, S) is the Robinson-Foulds distance (number of bipartitions in the
symmetric difference of T_i and S) and β > 0 is a global rate parameter
controlling how steeply the likelihood decays with disagreement. Large β implies
high confidence in the source trees; small β implies high tolerance for
disagreement.

The joint log-likelihood over all source trees, allowing per-tree weights w_i,
is:

    log L(S, β) = W·log(β) − β · Σ_i w_i · d(T_i, S)
                = W·log(β) − β · WD

where W = Σ_i w_i and WD = Σ_i w_i · d(T_i, S) is the weighted RF sum.

### 1.1 Closed-form MLE for β

Setting ∂ log L / ∂β = 0:

    W/β − WD = 0
    β̂ = W / WD

This is the unique global maximum (the log-likelihood is strictly concave in β
for WD > 0). No numerical optimisation is required. The maximum log-likelihood
is:

    log L(β̂) = W·log(W/WD) − W

### 1.2 Using mlscores

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

The `mlalpha` parameter provides a principled way to investigate and correct
this bias.

---

## 3. The tree-size scaling exponent α

### 3.1 Derivation from the probability model

Rather than bolting a scaling factor onto the distance after the fact, α is
derived by modifying the generative model itself. Suppose the rate of the
exponential for tree i is β / k_i^α rather than β:

    P(d_i | S, β, α) = (β / k_i^α) · exp(−(β / k_i^α) · d_i)

This is still a valid exponential distribution for each tree; the rate is now
tree-specific and depends on size. The joint log-likelihood is:

    log L(β, α) = n·log(β) − α · Σ_i log(k_i) − β · Σ_i w_i · d_i / k_i^α
                = n·log(β) − α · P − β · WD(α)

where P = Σ_i log(k_i) is a fixed penalty constant (determined by the data)
and WD(α) = Σ_i w_i · d_i · k_i^{−α}.

### 3.2 Why α does not degenerate

The term −α · P is a natural penalty that grows linearly with α. It arises
directly from the normalisation constant of the per-tree exponential: as the
rate β/k_i^α is made smaller (larger k_i^α), the distribution spreads out and
the probability of observing any particular d_i decreases. This penalty
counteracts the gain from reducing WD(α), creating a genuine interior optimum
for α.

This is the critical difference from naive weighting schemes. If one simply
multiplies d_i by an ad hoc factor that decreases with k_i, and then
re-estimates β, the likelihood will always increase as the factor approaches
zero because WD decreases without any corresponding penalty term. The α
parameter avoids this because it is derived from a fully specified probability
model.

### 3.3 Profiling β at each α

For any fixed α, the MLE for β is still closed-form:

    β̂(α) = W / WD(α)

Substituting back:

    log L(α) = W·log(W / WD(α)) − W − α · P

This reduces the joint optimisation over (β, α) to a 1-D search over α alone.
At each grid point α_k, WD(α_k) is computed from the raw RF distances and
β̂(α_k) is obtained immediately. The α with the highest log L(α) is the MLE.

### 3.4 Interpretation of α*

| α value | Interpretation |
|---------|----------------|
| 0 | Original Steel & Rodrigo (2008) model — raw RF, no size correction |
| (0, 1) | Partial correction — large trees still contribute more, but less than raw |
| 1 | Full normalisation — each tree contributes equally per unit weight, regardless of size |
| > 1 | Active down-weighting of large trees — appropriate if large trees are less reliable |
| → ∞ | All weight collapses onto the smallest trees (resisted by the −α·P penalty) |

The MLE condition ∂ log L / ∂α = 0 can be written as:

    n · 〈log k〉_{w(α)} = P = Σ_i log(k_i)

where 〈log k〉_{w(α)} is the weighted average of log(k_i) with weights
proportional to w_i · d_i · k_i^{−α}. This is a moment-matching condition: the
model fits when the weighted average log-size (weighted by disagreement) equals
the unweighted sum of log-sizes. Whether this has an interior solution depends
on the (d_i, k_i) distribution in the dataset.

If large trees tend to disagree more in absolute terms but less in relative
terms (d_i / k_i smaller for larger k_i), there is a meaningful α* > 0. If
relative disagreement is uniform across sizes, the penalty and gain cancel and
α* will be near 0.

### 3.5 Confidence intervals and identifiability

The curvature of log L(α) near α* determines the precision of the estimate. A
flat profile (low curvature) indicates that the data are not informative about
α — tree-size effects are absent or confounded. A sharp peak indicates a
well-identified size effect. The profile written by `mlscores outfile=` provides
the full likelihood surface for visual assessment.

---

## 4. Using α estimation in practice

### 4.1 Estimating α

    set criterion ml
    exe gene_trees.ph
    hs
    mlscores alpha=auto ascan=100 alphamax=3.0 outfile=profile.txt

This runs the standard β MLE first, then performs a 100-point grid search over
α ∈ [0, 3]. Both `ml_alpha` and `ml_beta` are updated to the joint MLE values.
The output file contains the β profile followed by the α profile.

### 4.2 Setting α manually

If you have a prior expectation or want to fix α for comparison:

    set mlalpha 1.0    # full normalisation
    hs                 # search uses alpha=1 throughout
    mlscores           # estimate beta conditional on alpha=1

### 4.3 Interpreting the result

- If α* ≈ 0: tree-size heterogeneity is not biasing your result; the original
  Steel model is appropriate.
- If α* ≈ 1: large trees dominate purely because of their size; full
  normalisation is warranted.
- If α* is intermediate: a partial size correction is supported by the data.
- If α* > 1: the data favour actively down-weighting large trees, suggesting
  that size correlates with unreliability in your dataset.
- If the optimal supertree topology changes substantially between α = 0 and α*,
  tree-size heterogeneity is a confounding factor and should be reported.

---

## 5. Relationship to tree weights

The `mlalpha` mechanism is conceptually distinct from per-tree weights (w_i).
Weights encode prior confidence in individual trees (e.g. from bootstrap
support via `autoweight=bs`). α encodes a systematic relationship between tree
size and reliability that applies uniformly across all trees. The two can be
used together: weights handle tree-specific variation; α handles the systematic
size trend.

---

## 6. Summary of commands

| Command / option | Effect |
|---|---|
| `mlscores` | Estimate β̂ = W/WD for current supertree; update `ml_beta` |
| `mlscores outfile=<f>` | Also write β log-likelihood profile to file |
| `mlscores alpha=auto` | Grid search over α ∈ [0, alphamax]; update `ml_alpha` and `ml_beta` |
| `mlscores ascan=<n>` | Number of grid points for α search (default 50) |
| `mlscores alphamax=<f>` | Upper bound of α grid (default 3.0) |
| `set mlalpha <f>` | Fix α to a specific value for subsequent hs/boot runs |
| `set mlbeta <f>` | Fix β to a specific value manually |

---

## 7. References

Steel, M. and Rodrigo, A., 2008. Maximum likelihood supertrees. *Systematic biology*, *57*(2), pp.243-250.

Akanni, W.A., Creevey, C.J., Wilkinson, M. and Pisani, D., 2014. LU St: a tool for approximated maximum likelihood supertree reconstruction. *BMC bioinformatics*
