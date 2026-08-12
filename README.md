# Cophylogenetic Reconciliation Algorithm

A computational biology project implemented for the **Evolution through Programming** course. This repository contains an implementation of a cophylogenetic reconciliation algorithm using the **CDTL model** (Cospeciation, Duplication, Transfer/Host-Switching, Loss) to map parasite phylogenetic trees onto host trees.

## Overview
The algorithm applies **Dynamic Programming** to compute the most parsimonious evolutionary history based on a configurable cost matrix for four key evolutionary events:
- **Cospeciation ($C_c$):** Simultaneous speciation of host and parasite lineages.
- **Duplication ($C_d$):** Parasite speciation on a single host branch.
- **Host-Switching ($C_t$):** Parasite transfer to a different host lineage.
- **Loss ($C_l$):** Extinction or absence of a parasite on a host branch.

The project analyzes how varying the penalty weight for **Host-Switching ($C_t$)** alters inferred evolutionary dynamics on both predefined tree topologies and empirical datasets

## Data
Data is taken from [here](https://www.cs.hmc.edu/~hadas/jane/SampleTrees.html)

## Visualization
`visualization.ipynb` draws a **tanglegram** for any dataset — the host tree and the parasite tree side by side, with a line connecting each parasite to the host it lives on. Crossing lines suggest host-switching (incongruence); near-parallel lines suggest cospeciation.

```python
draw_tanglegram('gopher_louse')
```

Available datasets: `gopher_louse`, `ficus_ceratosolen`, `seabird_louse`, `finches_african_brood_parasites`.

## Reconciliation algorithm
`reconciliation.ipynb` implements the DTL reconciliation with dynamic programming. For a parasite node `p` placed on host node `h`, it fills `cost[p][h]` bottom-up (leaves → root) and returns the cheapest total cost. It uses two tables: `cost` (`p` exactly at `h`) and `best_below` (`p` at `h` or below it, paying one **loss** per step down). A **host-switch** may target any host that is neither an ancestor nor a descendant of the current one.

```python
reconcile(host_tree, parasite_tree, mapping, C)   # -> minimum reconciliation cost
```

The notebook also **sweeps the host-switch penalty** (the optimal number of switches drops as it gets more expensive), **compares all datasets** using cost per parasite (a size-independent measure — finches/brood-parasites are the least faithful, ficus and gopher the most), and uses **traceback** to recover the individual events and colour each ancestor of the parasite tree by its event.
