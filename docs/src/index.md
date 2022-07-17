LovaszTheta.jl - Lovasz theta function and theta body for graphs
================================================================

LovaszTheta.jl provides functions for computing the Lovasz θ, Schrijver θ⁻, and Szegedy θ⁺
functions of a graph.  These provide upper bounds on the independence number of a graph and lower
bounds on the chromatic number of the complement graph.  They are homomorphism monotones, and
so can be used to provide necessary conditions on the existence of a homomorphism between a
pair of graphs.

Variations of these functions are available which accept a vector of vertex weights.  The
theta body is available as a semidefinite programming subroutine - it is possible to
constrain a Convex.jl variable to be within the theta body.

## Examples

```julia
using Graphs, LovaszTheta
@assert abs(θ(cycle_graph(5)) - √5) < 1e-7
```

Test that ``\max\left\{\sum w_i \middle| w \in \textrm{TH}(g) \right\} = \theta(g)``.

```julia
using Graphs, LovaszTheta, Convex, SCS
g = erdos_renyi(20, 0.5);
w = Variable(nv(g));
problem = maximize(sum(w), [w ∈ TH(g)]);
solve!(problem, () -> SCS.Optimizer(verbose=0, eps=1e-8))
@assert abs(problem.optval - θ(g)) < 1e-7
```

Test entropy splitting
([Entropy splitting for antiblocking corners and perfect graphs](https://link.springer.com/article/10.1007/BF02122693)).

```julia
using Graphs, LovaszTheta, Convex, SCS

function corner_entropy(p, corner)
    w = Variable(nv(g))
    problem = minimize(-p' * log(w), [w ∈ corner])
    solve!(problem, () -> SCS.Optimizer(verbose=0, eps=1e-8))
    return problem.optval
end

g = erdos_renyi(20, 0.5)
p = normalize(rand(nv(g)), 1)
ent = -p'*log.(p)
ce1 = corner_entropy(p, TH(g))
ce2 = corner_entropy(p, TH(complement(g)))
@assert abs(ent - (ce1 + ce2)) < 1e-7
```

More examples can be found in the unit tests.
