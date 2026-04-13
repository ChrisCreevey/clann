# Tree-space landscape analysis with `visitedtrees`

The `visitedtrees=<filename>` option on the `hs` command records every unique
topology visited during the heuristic search — across all replicates and threads
— to a tab-separated file. This enables post-hoc analysis of the score landscape
in tree space: which regions were explored, how often were they revisited, and
where do local and global optima sit relative to each other?

---

## 1. Generating the data

```bash
# CLI
clann hs gene_trees.ph criterion=dfit nreps=20 nthreads=8 visitedtrees=landscape.tsv

# Interactive / batch
printf 'exe gene_trees.ph\nhs nreps=20 nthreads=8 visitedtrees=landscape.tsv\nquit\n' | clann
```

At the end of the search Clann reports:
```
Number of unique topologies scored: 537
Visited topology landscape written to: landscape.tsv
  Unique topologies recorded: 537
```

---

## 2. Output format

The file has a header row followed by one row per unique topology:

```
newick	score	visit_count
(((((Human,Chimp),Gorilla),Orangutan),Macaque),((Mouse,Rat),(Dog,Cat)));	7.297618	163
((Gorilla,Orangutan),((Chimp,Human),Macaque),((Mouse,Rat),(Dog,Cat)));	21.706350	42
...
```

| Column | Type | Description |
|--------|------|-------------|
| `newick` | string | Named-taxon unrooted Newick string; each row is a distinct topology |
| `score` | float | Criterion score on first visit — deterministic for a given topology and source-tree set (lower = better for dfit/RF; more negative = better for ML) |
| `visit_count` | int | Total times this topology was proposed, including revisits within and across all replicates |

The sum of all `visit_count` values equals the total number of SPR moves scored
(unique + revisits) — a measure of search effort per region of tree space.

Quick sanity check:
```bash
awk -F'\t' 'NR>1{u++; v+=$3} END{print "unique:", u, "total visits:", v}' landscape.tsv
# unique: 537   total visits: 8421
```

The best-scoring topology should appear with a high visit count, confirming that
the search repeatedly converged on that region:

```bash
sort -t$'\t' -k2 -n landscape.tsv | head -5
```

---

## 3. Computing pairwise SPR distances

To embed topologies in a metric space you first need pairwise distances. SPR
distance is the natural choice — it mirrors the moves Clann makes during search.

| Tool | Notes |
|------|-------|
| `rspr` | Whidden & Matsen (2014); exact and approximate modes; fast for hundreds of trees |
| `tqdist` | Quartet-distance; good alternative metric |
| `clann sprdists` | Clann's built-in SPR distance command for loaded source trees |

Using `rspr`:
```bash
# Extract Newick strings (skip header)
tail -n +2 landscape.tsv | cut -f1 > visited_trees.ph

# Compute pairwise SPR distances — outputs an N×N matrix
rspr -pairwise -unrooted -no-symmetric-pairwise < visited_trees.ph > spr_matrix.txt
```

For hundreds of topologies `rspr` runs in seconds. For thousands (large taxon
sets, many reps) consider using approximate mode (`rspr -split`) or RF distance
as a faster proxy.

---

## 4. Visualising the landscape

Install Python dependencies:
```bash
pip install pandas numpy scikit-learn matplotlib scipy networkx
```

### 4.1 MDS scatter plot (recommended first step)

Non-metric MDS embeds the N×N SPR distance matrix in 2D. Each point is a unique
topology; size reflects visit frequency; colour reflects score.

```python
import pandas as pd
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

df = pd.read_csv("landscape.tsv", sep="\t")
D  = np.loadtxt("spr_matrix.txt")          # N×N SPR distance matrix

mds = MDS(n_components=2, metric=False, dissimilarity="precomputed",
          random_state=42, n_init=4, max_iter=500)
coords = mds.fit_transform(D)

fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(
    coords[:, 0], coords[:, 1],
    s=np.sqrt(df["visit_count"]) * 6,   # size ∝ visit frequency
    c=df["score"],
    cmap="viridis_r",                    # low score (better) = bright yellow
    alpha=0.75, linewidths=0.3, edgecolors="grey"
)
plt.colorbar(sc, ax=ax, label="Score (lower = better)")
ax.set_xlabel("MDS axis 1")
ax.set_ylabel("MDS axis 2")
ax.set_title(f"Tree-space landscape  ({len(df)} topologies, MDS on SPR distances)")
plt.tight_layout()
plt.savefig("landscape_scatter.png", dpi=150)
plt.show()
```

**How to read this plot:**
- **Point size** reflects how often a topology was visited: large points = attractors.
- **Point colour** reflects score: bright yellow = best (lowest), dark purple = worst.
- **Clusters of large bright points** are basins of attraction — regions the search
  converged on repeatedly. Ideally one dominant basin = the global optimum.
- **Isolated dark points** are topologies visited only a few times: dead-ends or
  brief local optima the search escaped quickly.

### 4.2 Score histogram (quick diagnostic)

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("landscape.tsv", sep="\t")

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

A narrow, unimodal histogram (most visits near the best score) indicates the
search converged well. A broad or bimodal distribution suggests two competing
regions; try increasing `nreps` or `sample`.

### 4.3 SPR neighbour graph

Nodes are unique topologies; edges connect pairs exactly one SPR move apart.
Node size reflects visit count; colour reflects score. Cliques reveal tight
neighbourhoods the search traversed freely; isolated nodes are topologies
reached only via a long path.

```python
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

df = pd.read_csv("landscape.tsv", sep="\t")
D  = np.loadtxt("spr_matrix.txt")

G = nx.Graph()
n = len(df)
for i in range(n):
    G.add_node(i, score=float(df["score"].iloc[i]),
                  visits=int(df["visit_count"].iloc[i]))
for i in range(n):
    for j in range(i + 1, n):
        if D[i, j] == 1:
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

### 4.4 3D score surface

MDS axes as x/y, score interpolated onto a grid as z. Valleys are good regions;
peaks are poor. Deep, narrow valleys are local optima the search may get trapped in.

```python
import pandas as pd
import numpy as np
from sklearn.manifold import MDS
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib.pyplot as plt

df = pd.read_csv("landscape.tsv", sep="\t")
D  = np.loadtxt("spr_matrix.txt")

mds    = MDS(n_components=2, metric=False, dissimilarity="precomputed",
             random_state=42, n_init=4, max_iter=500)
coords = mds.fit_transform(D)

xi = np.linspace(coords[:, 0].min(), coords[:, 0].max(), 200)
yi = np.linspace(coords[:, 1].min(), coords[:, 1].max(), 200)
xi, yi = np.meshgrid(xi, yi)
zi = griddata(coords, df["score"].values, (xi, yi), method="linear")

fig  = plt.figure(figsize=(10, 7))
ax   = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(xi, yi, zi, cmap="viridis_r", alpha=0.85)
fig.colorbar(surf, ax=ax, shrink=0.5, label="Score")
ax.set_xlabel("MDS 1"); ax.set_ylabel("MDS 2"); ax.set_zlabel("Score")
ax.set_title("Tree-space score surface")
plt.tight_layout()
plt.savefig("landscape_3d.png", dpi=150)
plt.show()
```

---

## 5. Identifying local optima

Sort by score to find the best-scoring topologies:
```bash
sort -t$'\t' -k2 -n landscape.tsv | head -11
```

Check whether the top topologies are genuinely distinct or all the same shape by
loading them into Clann and computing RF distances:
```bash
sort -t$'\t' -k2 -n landscape.tsv \
    | awk -F'\t' 'NR>1 && NR<=11{print $1}' > best10.ph
clann rfdists best10.ph
```

- All RF distances zero → single well-defined global optimum; the search converged.
- Non-zero RF distances among top topologies → several distinct topologies score
  similarly. Investigate further with `usertrees tests=yes` to determine whether
  any can be statistically rejected.

---

## 6. Tips for better landscape coverage

| Symptom | Cause | Solution |
|---------|-------|----------|
| Very few unique topologies recorded | Search converged quickly | Increase `nreps`; set `maxskips=0` to disable early stopping |
| All visits concentrated at one score | Landscape is simple and well-converged | Trust the result |
| Bimodal score histogram | Two competing basins of attraction | Increase `sample` for better starting trees; increase `nreps` |
| High visit count but few unique topologies | Search revisiting same small region | Use `swap=tbr` for more exploratory moves |
| Thousands of topologies (slow MDS) | Large taxon set or many reps | Use RF distance instead of SPR; or subsample to the N best-scoring topologies |

---

## 7. References

- Whidden, C. and Matsen, F.A., 2018. Calculating the unrooted subtree prune-and-regraft distance. IEEE/ACM transactions on computational biology and bioinformatics, 16(3), pp.898-911. — `rspr` SPR distance tool
- Estabrook, G.F., McMorris, F.R. and Meacham, C.A., 1985. Comparison of undirected phylogenetic trees based on subtrees of four evolutionary units. Systematic Zoology, 34(2), pp.193-200. — quartet distance (tqdist)
- Maddison, W.P., 1997. Gene trees in species trees. Systematic biology, 46(3), pp.523-536. — gene-tree / species-tree reconciliation
