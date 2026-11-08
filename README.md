# Optimized 2-Approximation of Treewidth — C++ implementation

Paper: M. Belbasi, M. Fürer, M. Kumar, *Optimized 2-Approximation of Treewidth*, arXiv:2411.16918v3
(<https://arxiv.org/abs/2411.16918>; the PDF is not included in this repository because it is copyrighted by the authors.)

Given a graph G, the program either constructs a tree decomposition of width ≤ 2k+1 or concludes that tw(G) > k (Theorem 1).
In the default (auto) mode it searches for k and outputs a tree decomposition of **width ≤ 2·tw(G)+1** together with a lower bound on tw.

All of the code in this repository was implemented by Claude Opus 5.5.

## Build and run

```
make                       # bin/tw2approx, bin/verify_td, bin/gen
./run_limited.sh ./bin/tw2approx graph.gr -o out.td     # run with a 32GB / 10 min limit
./bin/verify_td graph.gr out.td                          # verify the tree decomposition
tests/battery.sh                                         # randomized tests (with --check)
```

`run_limited.sh` runs a command under `ulimit -v 32GB` and `timeout 600`.

### Options of tw2approx

| Option | Meaning |
|---|---|
| `-k K` | Fix k. Outputs a tree decomposition of width ≤ 2K+1, or `tw > K` (exit code 1). Requires an initial width ≤ 4K+3 |
| (no `-k`) | Auto mode. Applies Theorem 1 repeatedly (see below) |
| `--init mindeg\|random\|mixed` | Elimination ordering used to build the initial tree decomposition: min-degree (default) / random / min-degree that eliminates a random vertex with probability `--mix Q` (for testing; deliberately produces a poor initial decomposition) |
| `--no-lb` | Do not use the lower bound tw ≥ degeneracy (for testing the failure path) |
| `--mem-gb G` | Memory limit for the DP tables (default 24GB). If the limit would be exceeded, the program exits with an error before starting the computation |
| `--check` | After every Move and Split, check the invariants of the grouped TD and compare **all DP tables with tables recomputed from scratch** (slow; for testing) |
| `-v` | Show the rounds and the values of k tried |

The input is in PACE `.gr` format and the output is in PACE `.td` format (only main nodes are output). Statistics are printed to stderr.

`bin/gen` is a test graph generator: `gen ktree N K P SEED` (a partial K-tree, tw ≤ K), `gen grid R C`, `gen gnm N M SEED`.

## Correspondence between the paper and the implementation

| Paper | Implementation (`src/tw2approx.cpp`) |
|---|---|
| Def 6 grouped TD / Lemma 7 | `GTD::build`. A tree decomposition built from an elimination ordering already has unique homes (the transformation of Lemma 5 is not needed). Intersection children are grouped by bag |
| Def 14 DP table `A_x` | `Node::E`. L^{|B|} entries; the label of the vertex at position p in the bag is the p-th digit of a base-L number (0 = X, 1..3 = W_i) |
| 3.1 Initialization / Addition / Subtraction / Update | `legal_init`, `add_im` (main←inter), `add_mi` (inter←main), with sign ±1 |
| 3.2 DFS instance | `dfs_find` |
| Alg 2 Move / Alg 3 Rotate | `move_to_root`, `rotate` (forget chains: `build_chain`) |
| Alg 4 Merge | The merge inside `rotate` and `attach_inter_*` (merges sibling inter nodes that have the same bag) |
| Def 10/11, Alg 5/6 Split | `split`, `process_pending`, `finish_split` |
| Remark 24 (3-way and 4-way partitions) | In a round with \|B_r\| = w+1 ≥ 3k+4, L=3 (C1,C2,S); otherwise L=4 (C1,C2,C3,S). All tables are recomputed when switching |

## Gaps in the paper filled and corrections made

1. **DP tables store the "excess"**: `E[P] = s(A[P]) − |X_P|`. With this representation, Addition becomes
   `E_m[P] += E_z[P|B_z]` for main←inter and `E_z[P] += min_l(E_y[P+v:l] + [l=X])` for inter←main.
   The formula `s(A_xi[P']) = s(A_x[P]) − |Q|` of Algorithm 5 becomes `E_xi[P'] = E_x[P]`, and Merge is a plain addition.
   All of these are equivalent to the formulas in the paper.

2. **Correction of the table-reuse formula in Algorithm 5 (important)**: The paper reuses A_x, after the editable children
   have been subtracted, for every copy x_i via `s(A_x[P]) − |Q|`. However, this A_x still contains **the contributions of
   the non-editable children z that move to another copy x_j**. When several labels are assigned to the X part of B_z,
   such a contribution is not 0. The table of x_i is then overestimated, and a split may be missed (tw > k may be reported wrongly).
   This implementation uses the reuse formula only for the copy j\* that receives the most non-editable children,
   and first subtracts from A_x the non-editable children that move to the other copies. The table of each copy j ≠ j\*
   is built by Initialization followed by adding its own non-editable children. `--check` confirms that all tables match the recomputed ones.

3. **Top-down choice of the split**: At the main child y of an editable intersection node z, only the label of
   the new vertex v can be chosen. The candidates are restricted to {X} ∪ a(z) (Def 11(b)).
   Merging a label outside a(z) into a label in a(z) keeps the partition legal and does not increase its size (the argument of Lemma 12).
   Hence the restriction does not change the "minimum size" in Def 10. Among choices of the same size, X is preferred (maximum |X|).
   At the root, among the partitions that satisfy `|W_i| + size < |B_r|` and `size ≤ k+1`, the one with the minimum size, and then the maximum |X|, is chosen.

4. **Assignment of non-editable children**: A non-editable child with exactly one non-empty W is assigned to the copy of that label.
   A non-editable child whose B_z is entirely X is assigned to j\* as "an arbitrary one in a(p(x))" (Def 11(a)).

5. **Copies and Merge**: A copy x_i of a non-root node x has a home in the new tree only when the label of the home vertex of x is i.
   The other copies are contained in the copies on the grandparent side, so they are not created but merged from the start.
   Their children and their table contributions are passed to the nearest existing ancestor copy.
   This gives the same result as the paper's "create the copies, then MERGE", but avoids unnecessary table computations.

6. **Restarting the DFS**: The restart rule of the paper ("from the first pre-visited node that was not edited") is ambiguous.
   Therefore, after a split, the DFS restarts from the root and skips post-visited (closed) nodes and their subtrees.
   A closed subtree is neither on a Move path nor in the edited region, so the closed mark remains valid after the split.

7. **Initial tree decomposition and the search for k**: The initial tree decomposition is built from an elimination ordering (min-degree).
   Instead of Bodlaender's reduction in Corollary 2, Theorem 1 is applied repeatedly as follows.
   The program keeps a lower bound lb (initially the degeneracy) and the current width cw, and runs with k = max(lb, ⌈(cw−3)/4⌉).
   On success, cw ≤ 2k+1. On failure, tw > k is known, so lb is set to k+1. The program stops when cw ≤ 2·lb+1.
   Therefore the 2-approximation is guaranteed even if the initial width exceeds 4k+3.

8. **Processing each connected component separately**: The grouped TD is defined for connected graphs (the root bag S must not be empty).
   Therefore the algorithm runs on each component, and the tree decompositions are joined at the end.

## Complexity and memory in practice

Every node stores a DP table of L^{|B|} int32 entries. Memory is therefore about n·L^{w+1}·4 bytes,
and one table operation takes L^{w+1} time (corresponding to O\*(81^k)n in the paper).
In practice the limit is an initial width of about 16 or less (k ≤ 3–4). Beyond that, the program checks the `--mem-gb` limit and exits with an error before computing.

## Test results (on the development machine, all under the 32GB / 10 min limit)

- `tests/battery.sh`: 75 cases in total — random partial K-trees (K=1..3), grids, and G(n,m) — run with `--check --init random --no-lb`.
  **0 failures**; 10 runs were stopped by the memory guard (random initialization gave initial widths of 19–33).
  Checked: the output is a valid tree decomposition, width ≤ 2k+1, width ≤ 2·tw+1, lower bound ≤ tw, and all DP tables match the recomputed ones.
- 25 cases with a build that enables `_GLIBCXX_DEBUG` and UBSan: no problems.
- Larger instances: a partial 3-tree with n=16000, initial width 15 (mixed initialization): width 7 with `-k 3`, 2.1 s, max RSS 354MB.
  Partial 2-trees with n=2000–8000, initial width 9–11: width 5 with `-k 2`, under 0.05 s each.

## License

MIT License. See [`LICENSE`](LICENSE) for the full text.
The paper (arXiv:2411.16918) is copyrighted by its authors and is not covered by this license.
