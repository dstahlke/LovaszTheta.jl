Usage
=====

| Function                  | Description                                                                                                 |
| :--------                 | :-----------                                                                                                |
| `θ(g)`                    | Lovasz theta function of the graph `g`.                                                                     |
| `θ(g, w)`                 | Weighted Lovasz theta function of the graph `g`. Weights must be non-negative.                              |
| `θ⁻(g)`                   | Schrijver theta function of the graph `g`.                                                                  |
| `θ⁻(g, w)`                | Weighted Schrijver theta function of the graph `g`. Weights must be non-negative.                           |
| `θ⁺(g)`                   | Szegedy theta function of the graph `g`.                                                                    |
| `θ⁺(g, w)`                | Weighted Szegedy theta function of the graph `g`. Weights must be non-negative.                             |
| `w ∈ a * TH(g)`           | Constrain Convex.Variable `w` to be within theta body (optionally) scaled by Variable `a`.                  |
| `theta_dual_weight(g, w)` | Compute the weight `v` that saturates the inequality `θ(g, w) * θ(complement(g), v) ≥ w' v`.                |
| `theta_solution(g, w)`    | Compute `θ(g, w)` and return the matrices associated with the optimal solution.  See docstring for details. |
