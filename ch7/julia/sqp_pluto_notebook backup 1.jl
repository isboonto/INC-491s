
### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 1a2b3c4d-5e6f-7a8b-9c0d-1e2f3a4b5c6d
begin
    using LinearAlgebra
    using ForwardDiff
end

# ╔═╡ 2b3c4d5e-6f7a-8b9c-0d1e-2f3a4b5c6d7e
begin
    ff = (x) -> -x[1] - x[2]
    gg1 = (x) -> x[1]^2 - x[2]
    gg2 = (x) -> x[1]^2 + x[2]^2 - 1
    gg = (x) -> [gg1(x), gg2(x)]
end

# ╔═╡ 3c4d5e6f-7a8b-9c0d-1e2f-3a4b5c6d7e8f
begin
    gradient_ff = x -> ForwardDiff.gradient(ff, x)
    gradient_gg1 = x -> ForwardDiff.gradient(gg1, x)
    gradient_gg2 = x -> ForwardDiff.gradient(gg2, x)
    hessian_ff = x -> ForwardDiff.hessian(ff, x)
    hessian_gg1 = x -> ForwardDiff.hessian(gg1, x)
    hessian_gg2 = x -> ForwardDiff.hessian(gg2, x)
end

# ╔═╡ 4d5e6f7a-8b9c-0d1e-2f3a-4b5c6d7e8f9a
begin
    function bfgs_update(H, s, y)
        rho = 1.0 / (y' * s)
        I = Matrix{Float64}(I, length(s), length(s))
        V = I - rho * s * y'
        H = V' * H * V + rho * s * s'
        return H
    end
end

# ╔═╡ 5e6f7a8b-9c0d-1e2f-3a4b-5c6d7e8f9a0b
begin
    function sqp(ff, gg1, gg2, x0; k=10)
        x = zeros(2, k)
        λ = zeros(2, k)
        d_opt = zeros(2, k)

        x[:, 1] = x0
        λ[:, 1] = [0.0, 0.0]

        for i in 2:k
            ∇f = gradient_ff(x[:, i-1])
            J = ForwardDiff.jacobian(gg, x[:, i-1])

            Hg1 = hessian_gg1(x[:, i-1])
            Hg2 = hessian_gg2(x[:, i-1])
            B = hessian_ff(x[:, i-1]) + λ[:, i-1]' * [Hg1, Hg2]

            p1 = gg1(x[:, i-1])
            p2 = gg2(x[:, i-1])
            ∇p1 = gradient_gg1(x[:, i-1])
            ∇p2 = gradient_gg2(x[:, i-1])

            # Solve the quadratic subproblem
            d = -inv(B) * ∇f  # Simplified step direction

            # Check constraints
            if p1 + ∇p1' * d > 0 || p2 + ∇p2' * d > 0
                println("Constraints violated at iteration $i")
                break
            end

            d_opt[:, i] = d
            x[:, i] = x[:, i-1] + d_opt[:, i]
            λ[:, i] = -[p1 + ∇p1' * d, p2 + ∇p2' * d]
        end

        return x, λ, d_opt
    end
end

# ╔═╡ 6f7a8b9c-0d1e-2f3a-4b5c-6d7e8f9a0b1c
begin
    x0 = [0.5, 1.0]
    x, λ, d_opt = sqp(ff, gg1, gg2, x0)
end

# ╔═╡ 7a8b9c0d-1e2f-3a4b-5c6d-7e8f9a0b1c2d
begin
    println("Optimal solution after 10 iterations: ", x[:, end])
    println("History of x values: ", x)
end
