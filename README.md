This repository contains the code for the paper "Cycles and cuts in supersingular $L$-isogeny graphs".
Each Jupyter notebook filed (.ipynb) is a SageMath 10.4 file.
The file contents (in alphabetical order) are as follows:

* `cycle_counts.ipynb`: code for counting the number of $\{\ell_1^{e_1},\dots,\ell_r^{e_r}\}$-isogeny cycles in the $L$-isogeny graph.
  1. Brandt Matrix Trace Formula (Theorem 24);
  2. The Brandt approximation (Equation (8));
  3. Ideal-theoretic count (only applies to $\{\ell_1,\dots,\ell_r\}$-isogeny cycles) (Theorem 33);
  4. Brute-force algorithm.
* `dist.ipynb`: code for computing the distance distributions in L-isogeny graphs.
* `example_nonprincipalcycle.ipynb`: code to accompany Example 23 in the paper, computing the canonical decomposition of a non-principal isogeny cycle.
* `fiedler.ipynb`: code to compute the edge expansion of the cut given by Fiedler's algorithm.
* `fiedler_comparison.ipynb`: code to compute the edge expansion for the following cuts: neighbour ordering, greedy neighbour ordering, and Fiedler (See Figure 7).
* `fiedler_viz.ipynb`: visualisations of graph cuts using reorderings of the vertices in the adjacency matrix of the graph (See Figure 8).
* `utils.sage`: This sage file contains helper functions for the computations in `cycle_counts.ipynb`, `dist.ipynb`, `fiedler.ipynb`, `fiedler_comparison.ipynb`, and `fiedler_viz.ipynb`.