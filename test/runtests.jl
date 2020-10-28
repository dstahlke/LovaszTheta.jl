using LovaszTheta
using LinearAlgebra
using Graphs
using GraphIO
using Convex, SCS
using Test
using Documenter
using Random

tol = 1e-7

doctest(LovaszTheta)

@testset "Basic tests" begin
    @test abs(θ(cycle_graph(5)) - √5) < tol

    g = loadgraph(IOBuffer("GkYYd?"), "graph1", GraphIO.Graph6.Graph6Format())

    @test abs(θ(g; cone=LovaszTheta.Schrijver) - 3.2360679766708524) < tol
    @test abs(θ(g; cone=LovaszTheta.Lovasz   ) - 3.3016943418683780) < tol
    @test abs(θ(g; cone=LovaszTheta.Szegedy  ) - 3.3380444547539740) < tol

    @test abs(θ⁻(g) - 3.2360679766708524) < tol
    @test abs(θ(g)  - 3.3016943418683780) < tol
    @test abs(θ⁺(g) - 3.3380444547539740) < tol
end

@testset "Duality" begin
    Random.seed!(0)
    g = erdos_renyi(20, 0.5)
    w = rand(nv(g))
    v = theta_dual_weight(g, w)
    @test abs(θ(complement(g), v) - 1) < tol
    @test abs(w'*v - θ(g, w) * θ(complement(g), v)) < tol
end

@testset "Trace TH" begin
    Random.seed!(0)
    g = erdos_renyi(20, 0.5)
    w = Variable(nv(g))
    problem = maximize(sum(w), [w ∈ TH(g)])
    solve!(problem, () -> SCS.Optimizer(verbose=0, eps=1e-8))
    @test abs(problem.optval - θ(g)) < tol
end

@testset "Product" begin
    Random.seed!(0)
    g = cycle_graph(5)
    g2 = union(cartesian_product(g, g), tensor_product(g, g)) # strong product
    a = rand(nv(g))
    b = rand(nv(g))
    @test abs(θ(g, a) * θ(g, b) - θ(g2, kron(a, b))) < tol
end

# Entropy splitting for antiblocking corners and perfect graphs
# Csiszar, Körner, Lovász, Marton
# DOI:10.1007/BF02122693
@testset "Entropy splitting" begin
    Random.seed!(0)

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
    @test abs(ent - (ce1 + ce2)) < tol
end

# Approximation of the Stability Number of a Graph Via Copositive Programming
# De Klerk, Pasechnik
# Also:
# A Sum of Squares Characterization of Perfect Graphs
# https://scirate.com/arxiv/2110.08950
@testset "De Klerk, Pasechnik" begin
    g = loadgraph(IOBuffer("GkYYd?"), "graph1", GraphIO.Graph6.Graph6Format())
    n = nv(g)
    J = ones(n, n)
    A = adjacency_matrix(g)
    f = k -> k*(A + I) - J
    λ = Variable()
    N = Variable(n, n)
    add_constraint!(λ, N >= 0)
    add_constraint!(λ, f(λ) - N ⪰ 0) # optimize over PSD+NONNEG, a subset of completely positive
    problem = minimize(λ)
    solve!(problem, () -> SCS.Optimizer(verbose=0, eps=1e-8))

    @test abs(problem.optval - θ⁻(g)) < tol
end

@testset "Knuth compatible matrices" begin
    Random.seed!(0)
    g = erdos_renyi(20, 0.5)
    eps = 1e-8
    tol = eps * 10
    w = rand(nv(g))
    v = theta_dual_weight(g, w, eps=eps)
    nt1 = theta_solution(g, w, eps=eps)
    nt2 = theta_solution(complement(g), v, eps=eps)
    @test abs(nt1.λ * nt2.λ - w' * v) < tol
    @test norm(nt2.Z * w - nt1.λ * v) < tol
    @test norm(nt1.Z * v - nt2.λ * w) < tol
    @test norm(nt1.Z * nt2.Z - w * v') < tol
end
