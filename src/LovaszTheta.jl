module LovaszTheta

export θ, θ⁻, θ⁺, TH, theta_dual_weight, theta_solution

using Convex, SCS
using Graphs
using LinearAlgebra
using MathOptInterface
import Base.*, Base.in

MOI = MathOptInterface

default_eps = 1e-8

@enum ThetaConeType Lovasz Schrijver Szegedy

function make_optimizer(eps)
    optimizer = SCS.Optimizer()
    if isdefined(MOI, :RawOptimizerAttribute) # as of MathOptInterface v0.10.0
        MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), 0)
        if isdefined(SCS, :ScsSettings) && hasfield(SCS.ScsSettings, :eps_rel) # as of SCS v0.9
            MOI.set(optimizer, MOI.RawOptimizerAttribute("eps_rel"), eps)
            MOI.set(optimizer, MOI.RawOptimizerAttribute("eps_abs"), eps)
        else
            MOI.set(optimizer, MOI.RawOptimizerAttribute("eps"), eps)
        end
    else
        MOI.set(optimizer, MOI.RawParameter("verbose"), 0)
        if isdefined(SCS, :ScsSettings) && hasfield(SCS.ScsSettings, :eps_rel) # as of SCS v0.9
            MOI.set(optimizer, MOI.RawParameter("eps_rel"), eps)
            MOI.set(optimizer, MOI.RawParameter("eps_abs"), eps)
        else
            MOI.set(optimizer, MOI.RawParameter("eps"), eps)
        end
    end
    return optimizer
end

"""
    θ(g::AbstractGraph, w::Convex.AbstractExpr)

Returns λ bounded from below by θ(g, w).  The constraint θ(g, w) ≤ 1 is equivalent to
w ∈ TH(complement(g)), the theta body of the complement graph.
"""
function θ(g::AbstractGraph, w::Convex.AbstractExpr; cone::ThetaConeType=Lovasz)
    # FIXME allow matrix w
    if size(w)[2] != 1
        throw(DimensionMismatch("Weight must be a vector (size(w) = $(size(w)))."))
    end
    if size(w) != (nv(g), 1)
        throw(DimensionMismatch("Weight vector length, $(size(w)[1]), didn't match number of vertices, $(nv(g))."))
    end

    λ = Variable()
    Z = Variable(nv(g), nv(g))

    add_constraint!(λ, w == diag(Z))
    if cone == Schrijver
        add_constraint!(λ, Z .* adjacency_matrix(complement(g)) <= 0)
    else
        add_constraint!(λ, Z .* adjacency_matrix(complement(g)) == 0)
    end
    add_constraint!(λ, [λ w' ; w Z] ⪰ 0)
    if cone == Szegedy
        add_constraint!(λ, Z >= 0)
    end

    return λ
end

"""
Represents the scaled (by non-negative variable λ) theta body of the graph g.
"""
struct ThetaBody
    g::AbstractGraph
    cone::ThetaConeType
    λ::Convex.AbstractExpr

    function ThetaBody(
        g::AbstractGraph,
        cone::ThetaConeType = Lovasz,
        λ::Convex.AbstractExpr = Constant(1.0))

        if size(λ) != (1,1)
            throw(DimensionMismatch("λ must be a scalar (size was $(size(λ)))"))
        end
        return new(g, cone, λ)
    end
end

"""
    TH(g::AbstractGraph)

Creates the theta body of the supplied graph.  The `*` operator can be used to scale this theta
body by a non-negative constant or by a scalar `Convex.Variable`.

julia> problem = maximize(sum(w), [w in TH(g)])

julia> problem = minimize(λ, [w in λ*TH(g)])
"""
TH(g::AbstractGraph, cone::ThetaConeType = Lovasz) = ThetaBody(g, cone)

*(λ::Convex.AbstractExpr, th::ThetaBody) = ThetaBody(th.g, th.cone, λ*th.λ)
*(λ::Real, th::ThetaBody) = Constant(λ) * th

*(th::ThetaBody, λ::Convex.AbstractExpr) = λ * th
*(th::ThetaBody, λ::Real) = λ * th

in(w::Convex.AbstractExpr, th::ThetaBody) = θ(complement(th.g), w, cone=th.cone) <= th.λ

"""
    θ(g::AbstractGraph)

Compute the Lovasz theta function for the given graph.
"""
function θ(g::AbstractGraph; eps=default_eps, cone::ThetaConeType=Lovasz)
    return θ(g, ones(nv(g)), eps=eps, cone=cone)
end

"""
    θ(g::AbstractGraph, w::AbstractArray{<:Number, 1})

Compute the weighted Lovasz theta function for the given graph.
"""
function θ(g::AbstractGraph, w::AbstractArray{<:Number, 1}; eps=default_eps, cone::ThetaConeType=Lovasz)
    if length(w) != nv(g)
        throw(DimensionMismatch("Weight vector length, $(length(w)), didn't match number of vertices, $(nv(g))."))
    end
    if minimum(w) < 0
        throw(DomainError("Weight vector must be non-negative."))
    end
    λ = θ(g, Constant(w); cone=cone)
    #λ = Variable()
    #add_constraint!(λ, w in λ * TH(complement(g), cone))
    problem = minimize(λ)
    solve!(problem, () -> make_optimizer(eps))
    if problem.status == MOI.INTERRUPTED
        throw(InterruptException())
    end
    return problem.optval
end

"""
    θ(g::AbstractGraph, w::AbstractArray{<:Number, 2})

Compute the matrix weighted Lovasz theta function for the given graph.
"""
function θ(g::AbstractGraph, w::AbstractArray{<:Number, 2}; eps=default_eps)
    x = Variable(nv(g))
    λ = θ(g, x)
    problem = minimize(λ, [ Diagonal(x) ⪰ w ])
    solve!(problem, () -> make_optimizer(eps))
    if problem.status == MOI.INTERRUPTED
        throw(InterruptException())
    end
    return problem.optval
end

"""
    θ(g::AbstractGraph)

Compute the Schrijver theta function for the given graph.
"""
θ⁻(g::AbstractGraph; eps=default_eps) = θ(g, eps=eps, cone=Schrijver)

"""
    θ(g::AbstractGraph, w::AbstractArray{<:Number, 1})

Compute the weighted Schrijver theta function for the given graph.
"""
θ⁻(g::AbstractGraph, w::AbstractArray{<:Number, 1}; eps=default_eps) = θ(g, w, eps=eps, cone=Schrijver)

"""
    θ(g::AbstractGraph)

Compute the Szegedy theta function for the given graph.
"""
θ⁺(g::AbstractGraph; eps=default_eps) = θ(g, eps=eps, cone=Szegedy)

"""
    θ(g::AbstractGraph, w::AbstractArray{<:Number, 1})

Compute the weighted Szegedy theta function for the given graph.
"""
θ⁺(g::AbstractGraph, w::AbstractArray{<:Number, 1}; eps=default_eps) = θ(g, w, eps=eps, cone=Szegedy)

"""
    theta_dual_weight(g::AbstractGraph, w::AbstractArray{<:Number, 1})

Compute the weight `v` that saturates the inequality
θ(g, w) * θ(complement(g), v) ≥ w' v.

The returned weight vector will be normalized such that θ(complement(g), v) = 1.
"""
function theta_dual_weight(g::AbstractGraph, w::AbstractArray{<:Number, 1}; eps=default_eps)
    v = Variable(nv(g))
    problem = maximize(w' * v, [ θ(complement(g), v) <= 1 ])
    solve!(problem, () -> make_optimizer(eps))
    if problem.status == MOI.INTERRUPTED
        throw(InterruptException())
    end
    # Clamp to positive orthant.  Solver can produce numbers slightly below zero.
    return max.(0, evaluate(v))
end

"""
    theta_solution(g::AbstractGraph, w::AbstractArray{<:Number, 1})

Compute the optimal matrices for the various formulations of θ.  Returns a named tuple
(λ, R, Y, Z, B, T) satisfying the following:

```
λ = θ(g, w)

R[i,j] = sqrt(w[i]*w[j])

Y ⪰ R
diag(Y) = λ*I
Y[i,j] = 0 unless i=j or i∼j

diag(Z) = w
[λ w'; w Z] ⪰ 0
Z[i,j] = 0 unless i=j or i∼j

B ⪰ 0
tr(B) = 1
tr(B*R) = λ
B[i,j] = 0 for i∼j

T ⪰ 0
diag(T) = w
|T| = λ
T[i,j] = 0 for i∼j
```
"""
function theta_solution(g::AbstractGraph, w::AbstractArray{<:Number, 1}; eps=default_eps)
    if length(w) != nv(g)
        throw(DimensionMismatch("Weight vector length, $(length(w)), didn't match number of vertices, $(nv(g))."))
    end
    if minimum(w) < 0
        throw(DomainError("Weight vector must be non-negative."))
    end

    tol = eps * 10
    λ = Variable()
    Y = Variable(nv(g), nv(g))
    R = Hermitian(sqrt.(w) * sqrt.(w)')

    constraints = Constraint[]
    push!(constraints, diag(Y) == λ*ones(nv(g)))
    push!(constraints, Y .* adjacency_matrix(complement(g)) == 0)
    push!(constraints, Y ⪰ R)

    problem = minimize(λ, constraints)
    solve!(problem, () -> make_optimizer(eps))
    if problem.status == MOI.INTERRUPTED
        throw(InterruptException())
    end

    λ = evaluate(λ)
    # Schur product with adjacency matrix is used to clean up near-zero values that should be
    # exactly zero.
    Y = Hermitian(Matrix(evaluate(Y) .* adjacency_matrix(g) + λ*I))

    Z = Hermitian(Y ./ λ .* R)
    @assert isposdef(Z + tol*I)
    @assert norm(diag(Z) - w) < tol
    @assert isposdef([λ w'; w Z] + tol*I)
    @assert norm(Z .* adjacency_matrix(complement(g))) < tol

    # Schur product with adjacency matrix is used to clean up near-zero values that should be
    # exactly zero.
    B = -diagm(constraints[1].dual[:,1]) - constraints[2].dual .* adjacency_matrix(complement(g))
    B = Hermitian((B + B') ./ 2)
    @assert isposdef(B + tol*I)
    @assert abs(tr(B) - 1) < tol
    @assert abs(tr(B * R) - λ) < tol
    @assert norm(B .* adjacency_matrix(g)) < tol

    invB = map(x -> x > tol ? 1/x : 0, diag(B))
    T = Hermitian(sqrt.(invB * invB') .* B .* R)
    T[diagind(T)] .= w
    @assert isposdef(T + tol*I)
    @assert abs(maximum(eigvals(T)) - λ) > -tol
    @assert norm(diag(T) - w) < tol
    @assert norm(T .* adjacency_matrix(g)) < tol

    return (λ=λ, R=R, Y=Y, Z=Z, B=B, T=T)
end

end # module
