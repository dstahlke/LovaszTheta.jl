using LovaszTheta
using LinearAlgebra
using Graphs
using GraphIO
using Convex, SCS
using Test
using Documenter
using Random

tol = 1e-7

make_optimizer = LovaszTheta.make_optimizer

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
    solve!(problem, () -> make_optimizer(1e-8))
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
        solve!(problem, () -> make_optimizer(1e-8))
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
    solve!(problem, () -> make_optimizer(1e-8))

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

@testset "Dual vector from T" begin
    # Compute theta_dual_weight(g, w) from theta_solution(g, w).T.
    # We take a maximization feasible point of θ(g, w) = λ and turn it into a minimization feasible
    # point of θ(complement(g), v) = 1 with v' w = λ.
    # This `v` and `Yc` should be returned as part of theta_solution, but I'd like to work out a
    # proof of correctness before doing so.
    Random.seed!(0)
    g = erdos_renyi(20, 0.5)
    w = rand(nv(g))
    sol = theta_solution(g, w)
    T = sol.T
    λ = sol.λ
    @test opnorm(T) ≈ λ  atol=tol

    a = eigvecs(T)[:,end] # largest eigenvector
    @test a' * T * a ≈ λ  atol=tol
    b = real(sqrt(T) * a)
    @test b' * b ≈ λ  atol=tol
    v = diagm(pinv.(sqrt.(w))) * b
    v = v .* v
    @test θ(complement(g), v) ≈ 1  atol=tol
    @test v' * w ≈ λ  atol=tol
    @test v ≈ theta_dual_weight(g, w)  atol = tol

    # Create a feasible solution for θ(complement(g), v) = 1.
    Yc = diagm(pinv.(sqrt.(w))) * T * diagm(pinv.(sqrt.(w)))
    @assert norm(adjacency_matrix(g) .* Yc) < tol
    @assert minimum(eigvals(Hermitian(Yc - sqrt.(v) * sqrt.(v)'))) > -tol  # Y ⪰ v v'
    @assert diag(Yc) ≈ ones(nv(g))  atol=tol
end

@testset "Homomorphism hoax" begin
    # Theorem 6 of arXiv:1310.7120, "Bounds on entanglement assisted
    # source-channel coding via the Lovasz theta number and its
    # variants."
    #
    # A semidefinite relaxation of homomorphism g → h exists iff
    # θ(complement(g)) <= θ(complement(h)).
    #
    # See also R. Bacik and S. Mahajan, "Semidefinite programming and
    # its applications to NP problems," 1995 where they call this
    # semidefinite relaxation a "hoax" and show it exists iff
    # θ(complement(g ∘ h)) = |V(g)| where g ∘ h is the hom-product.

    g = cycle_graph(5)
    h = complement(loadgraph(IOBuffer("GkYYd?"), "graph1", GraphIO.Graph6.Graph6Format()))

    # Places where Cxyst must be zero.
    # Theorem 3 of arXiv:1310.7120.
    #Agh = adjacency_matrix(homomorphic_product(g, h)) # FIXME use this when Graphs.jl gets homomorphic_product
    Agh = reshape([
            (has_edge(g, x, y) && !has_edge(h, s, t)) || (x == y && s != t)
            for
            s in vertices(h),
            x in vertices(g),
            t in vertices(h),
            y in vertices(g)
        ], nv(g)*nv(h), nv(g)*nv(h))

    solg = theta_solution(complement(g), ones(nv(g)))
    solh = theta_solution(complement(h), ones(nv(h)))

    # Solution exists when θ(complement(g)) <= θ(complement(h)).
    @test solg.λ < solh.λ

    eye(n) = Matrix(1.0*I, (n,n))

    λ = solh.λ
    B = solh.B
    D = Hermitian(B .* eye(nv(h)))
    @test minimum(eigvals(λ*D - B)) > -tol

    J = Hermitian(ones(nv(g), nv(g)))
    Z = Hermitian(solg.Y + (solh.λ - solg.λ)*eye(nv(g)) - J)
    @test minimum(eigvals(Z)) > -tol

    # Construct the hoax matrix.
    ⊗ = kron
    C = Hermitian(inv(λ) * (J ⊗ B + inv(λ-1) * Z ⊗ (λ * D - B)))

    # C ⪰ 0
    @test minimum(eigvals(C)) > -tol
    # Zeros for the forbidden input/output pairs.
    @test norm(C .* Agh) < tol
    # \sum_{st} C_{xyst} = J
    @test norm(reshape(sum(reshape(C, nv(h), nv(g), nv(h), nv(g)), dims=(1,3)), nv(g), nv(g)) - J) < tol
end
