# Clann v5.0 Tutorial

Copyright (C) 2003-2026 Chris Creevey <chris.creevey@gmail.com>
Distributed under the GNU General Public License v2 or later.
See COPYING.txt for full licence terms.

---

## Introduction

This tutorial walks through the major new features introduced in Clann v5.0,
using small, realistic mammalian gene-tree datasets included in the `examples/`
directory. You will need a compiled Clann binary on your PATH before starting.

**Tutorial datasets**

| File | Contents |
|------|----------|
| `examples/tutorial_single.ph` | 28 single-copy gene trees, 4–9 taxa (primates + rodents + carnivores) |
| `examples/tutorial_multicopy.ph` | 8 trees: 3 single-copy, 3 with prunable same-species clades, 2 with deep paralogs |
| `examples/tutorial_candidates.ph` | 4 candidate topologies for topology testing |

**Taxon set:** Human, Chimp, Gorilla, Orangutan, Macaque, Mouse, Rat, Dog, Cat.
The broadly accepted tree for this group is:
```
(((((Human,Chimp),Gorilla),Orangutan),Macaque),((Mouse,Rat),(Dog,Cat)));
```
Most gene trees in the tutorial datasets are consistent with this topology,
but with deliberate conflict to make the analyses non-trivial.

---

## Part 1: Your first supertree

Run a basic heuristic search using the default DFIT criterion:

```bash
clann hs examples/tutorial_single.ph
```

Clann will load the 28 gene trees, run 5 replicate searches (default `nreps=5`),
and print a result like:

```
Best score found: 7.3421
Best tree: (((((Human,Chimp),Gorilla),Orangutan),Macaque),((Mouse,Rat),(Dog,Cat)));
```

The **DFIT score** is the sum of normalised RF distances between the supertree
and each gene tree — lower is better. A score of 0 means perfect agreement with
all gene trees.

**Interactive mode** achieves the same result:
```
clann> exe examples/tutorial_single.ph
clann> hs nreps=5
```

---

## Part 2: Parallel search with OpenMP

Clann v5.0 can run independent search replicates in parallel using OpenMP.
By default it uses all available CPU cores. You can control this with `nthreads`:

```bash
# Use 4 threads (4 replicates in parallel)
clann hs examples/tutorial_single.ph nthreads=4 nreps=20

# Single-threaded for comparison / reproducibility
clann hs examples/tutorial_single.ph nthreads=1 nreps=20 seed=42
```

With a fixed `seed`, single-threaded runs are fully reproducible. Parallel runs
use per-thread seeds derived from the master seed and will differ between
thread counts, but the best score found should be the same or better than the
serial run given enough replicates.

**Parallel bootstrap**

The `boot` command also supports `nthreads`: independent bootstrap replicates
run in parallel, each doing a single-threaded heuristic search internally.

```bash
clann boot examples/tutorial_single.ph nreps=100 nthreads=8
```

**When does parallelism help?**
Parallelism scales linearly up to `nreps` threads. For small datasets like this
tutorial, I/O and tree-building dominate and the speedup is modest. On large
datasets (100+ gene trees, 30+ taxa) with `nreps=100`, using all available cores
can reduce wall-clock time from hours to minutes.

---

## Part 3: Comparing scoring criteria

Clann v5.0 adds two new scoring criteria alongside the existing DFIT, SFIT,
QFIT, and AVCON options.

### Robinson-Foulds (RF) criterion

The RF criterion counts shared bipartitions (splits) between the supertree
and each gene tree, normalised to [0, 1]. It is a fast, intuitive measure of
topological agreement.

```bash
clann hs examples/tutorial_single.ph criterion=rf nreps=10
```

### Maximum-likelihood (ML) criterion

The ML criterion uses an exponential likelihood model (Steel & Rodrigo 2008):
P(G_i | T) ∝ exp(−β · d_i), where d_i is the raw RF distance between the
supertree T and gene tree G_i. The score is reported as lnL = −β · Σd_i
(negative; higher is better, consistent with standard ML convention).

```bash
clann hs examples/tutorial_single.ph criterion=ml nreps=10
```

**Tuning the β parameter (`mlbeta`)**

β controls how steeply the likelihood penalises disagreement. Larger β
emphasises topological accuracy over consistency across trees:

```bash
# Gentler penalty
clann hs examples/tutorial_single.ph criterion=ml mlbeta=0.5 nreps=10

# Stronger penalty
clann hs examples/tutorial_single.ph criterion=ml mlbeta=2.0 nreps=10
```

**Score display (`mlscale`)**

By default scores are shown as lnL. Two alternative conventions are available:

```bash
# Paper convention (Steel & Rodrigo 2008): β·Σd_i, positive, lower is better
clann hs examples/tutorial_single.ph criterion=ml mlscale=paper

# Akanni et al. 2014 L.U.st convention
clann hs examples/tutorial_single.ph criterion=ml mlscale=lust
```

**Comparing criteria on the same data**

Run all three criteria and compare the resulting trees:

```bash
for crit in dfit rf ml; do
    echo "=== $crit ==="
    clann hs examples/tutorial_single.ph criterion=$crit nreps=10
done
```

For this dataset the inferred tree is typically the same across criteria
because the signal is strong. On real data with noisy or conflicting gene trees,
the choice of criterion can matter.

---

## Part 4: Testing alternative topologies (ML topology tests)

Given competing hypotheses for a phylogeny, `usertrees` with `tests=yes` performs
three standard statistical tests to determine whether the best-scoring topology
is significantly better than the alternatives.

The `examples/tutorial_candidates.ph` file contains four candidate topologies:

| Topology | Description |
|----------|-------------|
| T1 | Standard topology — matches expected mammalian tree |
| T2 | Human basal within great apes (Human sister to Chimp+Gorilla) |
| T3 | Gorilla sister to Orangutan (rather than to Human+Chimp) |
| T4 | Completely different placement of rodents and Orangutan |

```bash
clann usertrees \
    --source=examples/tutorial_single.ph \
    --topologies=examples/tutorial_candidates.ph \
    criterion=ml tests=yes nboot=1000
```

The positional form also works and is equivalent:
```bash
clann usertrees examples/tutorial_single.ph examples/tutorial_candidates.ph \
    criterion=ml tests=yes nboot=1000
```

**Output** (abridged):

```
ML topology tests  (T1 = tree 1,  lnL = -23.0000,  beta = 1.00)
=============================================================================
  Tree    lnL          WinSites p    KH z / p           SH p
-----------------------------------------------------------------------------
  T2      lnL=-26.000  p=0.1445 ns   z=1.18 p=0.119 ns   p=0.2200 ns
  T3      lnL=-28.000  p=0.0412 *    z=2.01 p=0.022 *     p=0.0480 *
  T4      lnL=-41.000  p=0.0020 **   z=4.53 p=0.000 ***   p=0.0010 ***
-----------------------------------------------------------------------------
  * p<alpha=0.05;  ** p<0.01;  *** p<0.001;  ns = not significant
```

**Interpreting the three tests:**

- **Winning Sites (WS):** Counts how many gene trees prefer T1 over the
  candidate (n+) versus the candidate over T1 (n−). Reports an exact
  two-sided binomial p-value. Simple and distribution-free.

- **Kishino–Hasegawa (KH):** A parametric z-test on the per-gene-tree
  log-likelihood differences. More powerful than WS when differences are
  consistent in magnitude, but assumes normality of the score distribution.

- **Shimodaira–Hasegawa (SH):** A bootstrap version of KH that corrects for
  the bias of having selected T1 as the best tree. More conservative than KH;
  the recommended test when comparing more than two topologies.

A detailed per-gene-tree breakdown is written to `mltest_results.txt`. Use
`testsfile=` to redirect it:

```bash
clann usertrees \
    --source=examples/tutorial_single.ph \
    --topologies=examples/tutorial_candidates.ph \
    criterion=ml tests=yes nboot=5000 testsfile=my_results.txt
```

---

## Part 5: Multicopy gene trees and autoprunemono

Real phylogenomic datasets often contain gene families with multiple copies per
species. By default, Clann excludes multicopy trees from supertree searches
(`hs`, `nj`, `alltrees`), but includes them in `reconstruct`. Many multicopy
trees contain only short, lineage-specific duplications: a monophyletic clade of
copies from the same species. These can be safely collapsed to a single
representative, often yielding a tree that is effectively single-copy.

`autoprunemono=yes` performs this pruning automatically at load time.

### Seeing the problem

Load the multicopy dataset without autoprunemono:

```bash
clann exe examples/tutorial_multicopy.ph
clann> hs nreps=5
```

Clann will report that only 3 of the 8 trees are being used (the 5 multicopy
trees are excluded by the single-copy filter).

### Enabling autoprunemono

```bash
clann exe examples/tutorial_multicopy.ph autoprunemono=yes
```

Output:
```
Autoprunemono: pruned monophyletic same-species clades in 3 multicopy trees.
  Promoted to single-copy pool:   3 trees
  Still multicopy after pruning:  2 trees (retained for reconstruct)
  Original (unpruned) trees stored for reconstruct.
```

Now run the search — all 6 prunable/single-copy trees contribute:

```bash
clann> hs nreps=5
```

### Gene-tree reconciliation

`reconstruct` maps each gene tree onto a species tree using a
duplication-and-loss model. It automatically uses the **original unpruned**
trees (so no information is lost), regardless of whether `autoprunemono` was
used:

```bash
# In the same session after hs:
clann> reconstruct speciestree=memory
```

Or point to an external species tree:

```bash
clann> reconstruct speciestree=examples/tutorial_single.ph dups=1.0 losses=1.0
```

Output files include:
- `*.recon` — reconciliation mapping per gene tree
- `*.onetoone` — putative one-to-one orthologs
- `*.recon.dist` — distribution of duplication/loss events

**Note on the delimiter:** The multicopy tutorial trees use `.` (dot) as
the species delimiter, which is the Clann default. `Human.alpha` and `Human.beta`
are treated as two copies of species `Human`. If your data uses a different
delimiter, set it with `set delimiter_char=X` before loading trees.

---

## Part 6: Exploring tree space with visitedtrees

This part demonstrates the `visitedtrees=` option, which records every unique
topology visited during the heuristic search — across all replicates and threads
— to a tab-separated file. The resulting data can be used for post-hoc analysis
of the score landscape in tree space, helping to distinguish local optima from
the global optimum and understand search behaviour.

### 6.1 Running the search with landscape recording

```bash
clann hs examples/tutorial_single.ph \
    criterion=dfit nreps=20 nthreads=4 \
    visitedtrees=landscape.tsv
```

Or in interactive mode:
```
exe examples/tutorial_single.ph
hs nreps=20 nthreads=4 visitedtrees=landscape.tsv
```

At the end of the search, Clann reports something like:
```
Number of unique topologies scored: 537
Visited topology landscape written to: landscape.tsv
  Unique topologies recorded: 537
```

### 6.2 Understanding the output file

The file `landscape.tsv` has three tab-separated columns:

```
newick	score	visit_count
((Gorilla,Orangutan),((Chimp,Human),Macaque),((Mouse,Rat),(Dog,Cat)));	21.706350	163
(((Cat,Mouse),((Dog,Gorilla),(Human,Chimp))),Macaque,(Rat,Orangutan));	34.956352	2
...
```

- **`newick`** — named-taxon unrooted Newick string; each row is a distinct topology
- **`score`** — criterion score on first encounter (lower = better for dfit/RF; more negative = better for ML)
- **`visit_count`** — total times this topology was proposed, including revisits within and across replicates

A quick sanity check — total visits should far exceed the number of unique topologies,
confirming that the search revisited regions of tree space:

```bash
awk -F'\t' 'NR>1{u++; v+=$3} END{print "unique:", u, "total visits:", v}' landscape.tsv
# unique: 537   total visits: 8421
```

The best-scoring topology (the supertree) should appear in the file with a high
`visit_count`, indicating that the search repeatedly converged on that region:

```bash
sort -t$'\t' -k2 -n landscape.tsv | head -3
```

### 6.3 Computing pairwise SPR distances

To embed the visited topologies in a metric space, compute pairwise SPR distances
between them. The SPR distance is particularly meaningful here because it mirrors
the search moves Clann makes. Any tree-distance tool can be used; we illustrate
with `rspr` (Whidden & Matsen 2014, freely available):

```bash
# Extract just the Newick strings (skip header line)
tail -n +2 landscape.tsv | cut -f1 > visited_trees.ph

# Compute pairwise SPR distances
# rspr outputs an N×N matrix
rspr visited_trees.ph > spr_matrix.txt
```

For datasets with hundreds of topologies, `rspr` runs in seconds. For thousands
of topologies (large taxa sets, many reps), consider using the approximate mode
or a quartet/RF distance as a faster proxy.

### 6.4 Visualising the landscape

Install the required Python packages if needed:

```bash
pip install pandas numpy scikit-learn matplotlib scipy networkx
```

#### Scatter plot (recommended first step)

Save as `landscape_plot.py` and run with `python3 landscape_plot.py`:

```python
import pandas as pd
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

# Load Clann landscape output
df = pd.read_csv("landscape.tsv", sep="\t")

# Load pairwise SPR distance matrix (N×N, no header)
D = np.loadtxt("spr_matrix.txt")

# Non-metric MDS: embed N topologies in 2D using SPR distances
mds = MDS(n_components=2, metric=False, dissimilarity="precomputed",
          random_state=42, n_init=4, max_iter=500)
coords = mds.fit_transform(D)

# Scatter plot
fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(
    coords[:, 0], coords[:, 1],
    s=np.sqrt(df["visit_count"]) * 6,   # size ∝ visit frequency
    c=df["score"],
    cmap="viridis_r",                    # low score (better) = bright yellow
    alpha=0.75, linewidths=0.3, edgecolors="grey"
)
plt.colorbar(sc, ax=ax, label="Score (lower = better for dfit)")
ax.set_xlabel("MDS axis 1")
ax.set_ylabel("MDS axis 2")
ax.set_title(f"Tree-space landscape  ({len(df)} topologies, MDS on SPR distances)")
plt.tight_layout()
plt.savefig("landscape_scatter.png", dpi=150)
plt.show()
```

**How to read this plot:**
- Each point is a unique topology the search visited.
- **Point size** reflects how many times it was visited: large points = attractors.
- **Point colour** reflects score: bright yellow = best (lowest) score, dark purple = poor score.
- **Clusters** of large, bright points are basins of attraction — regions the search
  repeatedly converges on. Ideally there is one dominant basin (the global optimum).
- **Isolated dark points** are topologies visited only a few times (dead-ends or brief
  local optima the search escaped quickly).

#### Score histogram (quick diagnostic)

```python
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(df["score"], bins=40, weights=df["visit_count"],
        color="steelblue", edgecolor="white", linewidth=0.5)
ax.set_xlabel("Score (lower = better)")
ax.set_ylabel("Total visits (visit_count weighted)")
ax.set_title("Distribution of search effort across score values")
plt.tight_layout()
plt.savefig("score_histogram.png", dpi=150)
plt.show()
```

A narrow, unimodal histogram — most visits concentrated near the best score —
indicates the search converged well. A broad or bimodal histogram suggests the
search is split between two competing regions; try increasing `nreps`.

#### Neighbour graph

Edges connect topologies that are exactly one SPR move apart. Cliques in the
graph reveal tight neighbourhoods that the search traversed freely:

```python
import networkx as nx

G = nx.Graph()
n = len(df)
for i in range(n):
    G.add_node(i, score=float(df["score"].iloc[i]),
                  visits=int(df["visit_count"].iloc[i]))
for i in range(n):
    for j in range(i + 1, n):
        if D[i, j] == 1:           # directly adjacent in SPR space
            G.add_edge(i, j)

sizes  = [G.nodes[v]["visits"] ** 0.5 * 15 for v in G]
colors = [G.nodes[v]["score"] for v in G]
pos    = nx.spring_layout(G, seed=42, k=0.5)
fig, ax = plt.subplots(figsize=(9, 7))
nx.draw_networkx(G, pos, ax=ax, node_size=sizes, node_color=colors,
                 cmap="viridis_r", with_labels=False,
                 edge_color="lightgrey", width=0.5)
plt.colorbar(plt.cm.ScalarMappable(cmap="viridis_r"), ax=ax, label="Score")
ax.set_title("SPR neighbour graph of visited topologies")
ax.axis("off")
plt.tight_layout()
plt.savefig("neighbour_graph.png", dpi=150)
plt.show()
```

#### 3D score surface

A surface plot adds height (score) to the MDS layout, giving an intuitive
"fitness landscape" view:

```python
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

xi = np.linspace(coords[:, 0].min(), coords[:, 0].max(), 200)
yi = np.linspace(coords[:, 1].min(), coords[:, 1].max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata(coords, df["score"].values, (xi, yi), method="linear")

fig = plt.figure(figsize=(10, 7))
ax  = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(xi, yi, zi, cmap="viridis_r", alpha=0.85)
fig.colorbar(surf, ax=ax, shrink=0.5, label="Score")
ax.set_xlabel("MDS 1"); ax.set_ylabel("MDS 2"); ax.set_zlabel("Score")
ax.set_title("Tree-space score surface")
plt.tight_layout()
plt.savefig("landscape_3d.png", dpi=150)
plt.show()
```

Valleys in the surface correspond to low-score (good) regions; peaks are
high-score (poor) regions. Deep, narrow valleys are local optima the search
may get trapped in.

### 6.5 Identifying local optima

Sort the landscape file by score and look at the top-scoring topologies:

```bash
# Best 10 topologies by score
sort -t$'\t' -k2 -n landscape.tsv | head -11
```

To check whether the top-scoring topologies are all the same shape or genuinely
different topologies, load them back into Clann:

```bash
# Write best topologies to a file
sort -t$'\t' -k2 -n landscape.tsv | awk -F'\t' 'NR>1 && NR<=11{print $1}' > best10.ph
clann rfdists best10.ph
```

If the RF distances among the top 10 are all zero, the search found a single
well-defined optimum. Non-zero RF distances indicate that several distinct
topologies score nearly as well — worth investigating with `usertrees tests=yes`.

### 6.6 Tips for better landscape coverage

| Problem | Solution |
|---------|----------|
| Too few unique topologies recorded | Increase `nreps`; decrease `maxskips` or set `maxskips=0` |
| All visits concentrated at one score | Search converged — the landscape is simple; trust the result |
| Bimodal score histogram | Two competing optima; increase `sample` to improve starting trees |
| Very high total visit count but few unique topologies | Search is revisiting the same region — add `swap=tbr` for more exploratory moves |

---

## Part 7: Using Clann in pipelines

Clann v5.0's direct CLI interface makes it straightforward to integrate into
analysis pipelines.

### One-liners

```bash
# Heuristic search, 20 replicates, 8 threads, ML criterion
clann hs gene_trees.ph criterion=ml nreps=20 nthreads=8

# Test candidate topologies
clann usertrees gene_trees.ph candidates.ph criterion=ml tests=yes nboot=1000

# Load with autoprunemono, then search (interactive batch)
printf 'exe data.ph autoprunemono=yes\nhs nreps=20 criterion=ml\nquit\n' | clann
```

### Snakemake example

```python
rule supertree:
    input:
        trees = "results/{dataset}/gene_trees.ph"
    output:
        result = "results/{dataset}/Heuristic_result.txt"
    params:
        nreps  = 20,
        nthreads = 8,
        criterion = "ml"
    shell:
        """
        clann hs {input.trees} \
            criterion={params.criterion} \
            nreps={params.nreps} \
            nthreads={params.nthreads}
        mv Heuristic_result.txt {output.result}
        """

rule topology_tests:
    input:
        trees      = "results/{dataset}/gene_trees.ph",
        candidates = "results/{dataset}/candidates.ph"
    output:
        tests = "results/{dataset}/mltest_results.txt"
    shell:
        """
        clann usertrees \
            --source={input.trees} \
            --topologies={input.candidates} \
            criterion=ml tests=yes nboot=1000 \
            testsfile={output.tests}
        """
```

### Exit codes

Clann returns:
- `0` — success
- `1` — error (bad arguments, file not found, parse error)

---

## Reference: Output files

| File | Produced by | Contents |
|------|------------|---------|
| `Heuristic_result.txt` | `hs` | Best supertree topology and score |
| `landscape.tsv` (user-named) | `hs visitedtrees=` | All visited topologies: newick, score, visit_count (see Part 6) |
| `NJ-tree.ph` | `nj` | Neighbour-joining supertree |
| `consensus.ph` | `consensus` | Consensus tree |
| `Usertrees_result.txt` | `usertrees` | Score for each candidate topology |
| `mltest_results.txt` | `usertrees tests=yes` | Per-gene-tree δ table for WS/KH/SH tests |
| `sourcetree_scores.txt` | `usertrees` | Per-source-tree contribution to each candidate score |
| `*.recon` | `reconstruct` | Reconciliation mapping per gene tree |
| `*.recon.dist` | `reconstruct` | Duplication/loss counts per gene tree |
| `*.onetoone` | `reconstruct` | Putative one-to-one orthologs |
| `*.strict.onetoone` | `reconstruct` | Strict one-to-one orthologs |
| `*.gene.births` | `reconstruct` | Gene birth-point assignments |

---

## Further reading

- Creevey & McInerney (2005) *Bioinformatics* 21(3):390–392 — original Clann paper
- Steel & Rodrigo (2008) *Systematic Biology* 57(2):197–211 — ML supertree criterion
- Kishino & Hasegawa (1989) *J Mol Evol* 29(2):170–179 — KH topology test
- Shimodaira & Hasegawa (1999) *Mol Biol Evol* 16(8):1114–1116 — SH topology test
- Whidden & Matsen (2014) *IEEE/ACM TCBB* 12(1):9–20 — fast SPR distance computation (`rspr`)

For help and bug reports: https://github.com/ChrisCreevey/clann/issues
