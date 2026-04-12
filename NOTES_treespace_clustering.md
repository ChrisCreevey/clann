# Tree-space Topology Clustering

The `clusterlandscape` option on `hs` (and the standalone `recluster` command)
groups the topologies recorded in the landscape file into clusters of similar
trees. This note explains the algorithm, its parameters, the output format, and
practical guidance for choosing a threshold.

---

## 1. Why cluster?

A landscape file from a long `hs` run can contain hundreds or thousands of
distinct topologies. Most represent minor rearrangements around the same
fundamental solution. Clustering replaces this cloud of topologies with a small
table of cluster representatives, each summarising a distinct region of tree
space. This makes it easy to answer questions such as:

- Was there one dominant attractor, or do several topologically distinct trees
  score similarly?
- How much of the search effort was spent around the best tree vs alternative
  solutions?
- Are there well-separated secondary optima that warrant further investigation?

---

## 2. Algorithm overview

The implementation is in `treecluster.c`, function `lm_cluster`.

### Step 1 — Sort entries

All topologies in the landscape map are collected into a flat array and sorted
by one of two keys (controlled by the `clusterorderby` option):

| `clusterorderby` | Sort order | Effect |
|------------------|------------|--------|
| `score` (default) | ascending (best first) | The best-scoring topology always becomes a cluster representative |
| `visits` | descending (most-visited first) | The most-visited (most probable basin) topology becomes a representative |

### Step 2 — Greedy sweep

The algorithm makes a single left-to-right pass through the sorted array:

```
clusters = []

for each topology T in sorted order:

    compute bipartitions(T)

    best_rf   = infinity
    best_clus = none

    for each existing cluster C:
        rf = normalised_RF(bipartitions(T), C.representative_bipartitions)
        if rf < best_rf:
            best_rf   = rf
            best_clus = C

    if best_clus exists AND best_rf <= threshold:
        assign T to best_clus
        update best_clus.best_score if T.score < best_clus.best_score
        best_clus.member_count += 1
        best_clus.total_visits += T.visit_count
    else:
        open a new cluster with T as representative
        record bipartitions(T) as the representative bipartitions
```

This is a **greedy single-pass algorithm**: once a cluster is opened, its
representative never changes. Topologies are assigned to the closest existing
cluster whose representative is within the threshold. If no cluster is close
enough, a new one is started. The order matters: because entries are sorted
by score (or visits) before the sweep, the representative of each cluster is
always the best-scoring (or most-visited) topology in that cluster.

### Step 3 — Sort and write clusters

After the sweep, clusters are sorted by `member_count` descending (largest
first) and written to the output TSV file.

---

## 3. Bipartition hashing

Comparing topologies is done entirely through sorted arrays of 64-bit
bipartition hashes — the same scheme used throughout Clann for RF distance
scoring. The full method is described in
[NOTES_rf_hash_method.md](NOTES_rf_hash_method.md); the key points for
clustering are:

- Each taxon `i` has a splitmix64-derived weight `taxon_hash_vals[i]`.
- A bipartition is represented as `min(XOR of left-side weights,
  XOR of right-side weights)` — a single canonical 64-bit integer.
- An unrooted binary tree on *n* leaves has exactly *n* − 3 non-trivial
  bipartitions; `collect_biparts_named` extracts all of them with a single
  stack-based pass over the Newick string.
- The two sorted arrays for any pair of trees are compared in O(*n*) time
  with a merge-walk, counting the number of shared bipartitions.

The **normalised RF distance** between two trees is:

```
raw_RF    = |biparts(T1)| + |biparts(T2)| - 2 * |biparts(T1) ∩ biparts(T2)|
max_RF    = 2 * max(|biparts(T1)|, |biparts(T2)|)
            = 2 * (n - 3)   for fully bifurcating trees on n leaves
norm_RF   = raw_RF / max_RF   ∈ [0, 1]
```

Two identical topologies have `norm_RF = 0`; two maximally different
topologies (no shared bipartitions) have `norm_RF = 1`.

---

## 4. The `clusterthreshold` parameter

`clusterthreshold` is the maximum normalised RF distance at which a topology
is assigned to an existing cluster rather than starting a new one. It controls
the granularity of the clustering:

| Value | Effect |
|-------|--------|
| 0.0 | Only topologically identical trees are grouped; every distinct topology is its own cluster |
| 0.1 | Topologies that differ by up to 10% of the maximum possible RF distance are grouped; tight clusters around nearly-identical trees |
| 0.2 (default) | Moderate grouping; works well for most datasets |
| 0.5 | Loose grouping; topologies that share more than half their bipartitions are placed together |
| 1.0 | All topologies are forced into one cluster |

**Choosing a good threshold:**

- Start with the default (0.2) and inspect the number of clusters and the
  `member_count` distribution in the output TSV.
- If you get one giant cluster with everything in it, try a smaller threshold.
- If you get hundreds of clusters with one member each, try a larger threshold.
- A well-chosen threshold typically produces a small number of clusters (5–20)
  with a clear dominant cluster (high `member_count` and `total_visits`)
  corresponding to the main attractor, plus a few smaller clusters for secondary
  optima.
- Use the standalone `recluster` command to re-run clustering at different
  thresholds without re-running the search (see Section 6).

---

## 5. Output format

The cluster TSV has a header row followed by one row per cluster, sorted by
`member_count` descending (most populated cluster first):

```
cluster_id	rep_newick	member_count	total_visits	best_score	rep_score
1	((Human,Chimp),(Gorilla,Orangutan),Macaque);	312	4821	7.1234	7.1234
2	((Human,(Chimp,Gorilla)),(Orangutan,Macaque));	84	1103	7.8901	7.8901
3	((Human,Chimp),(Gorilla,(Orangutan,Macaque)));	17	212	8.2341	8.4102
...
```

| Column | Type | Description |
|--------|------|-------------|
| `cluster_id` | int | Sequential integer label (1-based), in order of decreasing `member_count` |
| `rep_newick` | string | Newick string of the cluster representative (best-scoring or most-visited topology, depending on `clusterorderby`) |
| `member_count` | int | Number of distinct topologies in this cluster |
| `best_score` | float | Score of the best-scoring topology in the cluster (may differ from `rep_score` when `clusterorderby=visits`) |
| `rep_score` | float | Score of the representative topology |
| `total_visits` | int | Sum of `visit_count` across all member topologies — a measure of how much search effort this region of tree space attracted |

**Note on `best_score` vs `rep_score`:** when `clusterorderby=score` (the
default), the representative is always the best-scoring topology in its
cluster, so `best_score == rep_score`. When `clusterorderby=visits`, the
representative is the most-visited topology, which may not be the best-scoring
one; in that case `best_score ≤ rep_score` (lower is better).

---

## 6. Running clustering inline vs standalone

### Inline (within `hs`)

```
hs nreps=20 nthreads=8 visitedtrees=landscape.tsv \
   clusterlandscape=yes clusterthreshold=0.2 clusteroutput=clusters.tsv
```

Clustering runs immediately after the search finishes, while the landscape
map is still in memory. No disk round-trip occurs for the bipartition data.

### Standalone (`recluster`)

```
recluster landscape.tsv clusterthreshold=0.1 clusteroutput=tight_clusters.tsv
```

or via the CLI:

```bash
clann recluster landscape.tsv clusterthreshold=0.3 clusterorderby=visits
```

`recluster` reads the landscape TSV back into a fresh `LandscapeMap` and then
runs `lm_cluster` exactly as the inline path does. This lets you experiment
with different thresholds and orderings without re-running the (potentially
slow) heuristic search.

When no gene trees are loaded, `recluster` derives the taxon list automatically
from the Newick strings in the landscape file, so no `exe` command is required.

---

## 7. Computational complexity

The greedy algorithm is **O(N × K)** where *N* is the number of unique
topologies in the landscape and *K* is the number of clusters opened so far.
Each iteration compares the current topology against all existing cluster
representatives.

In practice *K* ≪ *N* for reasonable thresholds, so the algorithm is nearly
linear in *N*. For pathological cases (very small threshold on a large, diverse
landscape), *K* approaches *N* and the runtime approaches O(N²). Each
individual comparison is O(n) in the number of taxa (merge-walk on sorted
bipartition arrays), so the total worst-case runtime is O(N² × n).

For typical `hs` outputs (a few hundred to a few thousand unique topologies,
threshold 0.1–0.3) the clustering completes in milliseconds.

---

## 8. Relationship to other clustering methods

The algorithm is a form of **leader algorithm** or **sequential clustering**:
the first topology in the sorted order becomes the leader of the first cluster;
each subsequent topology either joins the nearest existing cluster (if within
threshold) or starts a new one. It is deterministic and order-dependent.

This differs from hierarchical clustering (which builds a full linkage tree)
and k-means (which requires specifying *k* in advance). The main advantages
for tree-space clustering are:

- **No need to specify k** — the number of clusters emerges naturally from
  the threshold and the data.
- **Consistent representatives** — the representative of each cluster is
  always the first (best-scoring) topology assigned to it, not a centroid.
  Centroids are not meaningful for tree topologies.
- **Speed** — a single pass is sufficient; no iterations are required.

The main disadvantage is that the result depends on the sort order and the
threshold, and there is no guarantee of optimality (a globally optimal
partition would require comparing all pairs of topologies). For the purpose
of summarising the search landscape this is an acceptable trade-off.

---

## 9. References

- Robinson & Foulds (1981) *Math Biosci* 53(1-2):131–147 — RF distance
- Estabrook *et al.* (1985) *Math Biosci* 77:207–216 — quartet distance
- Whidden & Matsen (2014) *IEEE/ACM TCBB* 12(1):9–20 — SPR distance tool (`rspr`)
- The bipartition hashing method is described in detail in
  [NOTES_rf_hash_method.md](NOTES_rf_hash_method.md)
- The landscape recording workflow is described in
  [NOTES_treespace_landscape.md](NOTES_treespace_landscape.md)
