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

```jldoctest
julia> using Graphs, LovaszTheta

julia> abs(θ(cycle_graph(5)) - √5) < 1e-7
true
```

Test that ``\max\left\{\sum w_i \middle| w \in \textrm{TH}(g) \right\} = \theta(g)``.

```jldoctest
julia> using Graphs, LovaszTheta, Convex, SCS

julia> g = erdos_renyi(20, 0.5);

julia> w = Variable(nv(g));

julia> problem = maximize(sum(w), [w ∈ TH(g)]);

julia> solve!(problem, () -> SCS.Optimizer(verbose=0, eps=1e-8))

julia> abs(problem.optval - θ(g)) < 1e-7
true
```

More examples can be found in the unit tests.
