### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ d42f2a4e-e0ac-4eb3-9b24-cf66f0df4bd6
using Pkg; Pkg.activate("/home/isudonto/.julia/environments")

# ╔═╡ acd5768b-185f-46e5-8d15-2694809d0a90
begin
	using LinearAlgebra, ForwardDiff
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using Symbolics
	using LaTeXStrings
end

# ╔═╡ fd9add82-4ab2-4ef7-a2d6-b14e978e77d2
begin
	f =  (x) -> x[1]^2 - x[1]*x[2] + x[2]^2 - 3x[1]
	a1 = (x) -> x[1] + x[2] - 2
	a2 = (x) -> -x[1]
	a3 = (x) -> -x[2]
	a4 = (x) -> x[1] - 3/2
	ak = x -> [a1(x), a2(x), a3(x), a4(x)]
end

# ╔═╡ 250253c0-eb66-4114-b0fd-ea4c64b5d5c1
begin
	@variables x₁ x₂
	
	x = [x₁, x₂]
	ak([x₁, x₂])
	x0 = [0, 0]  				
	@show gk = x -> ForwardDiff.gradient(f, x)
	@show A0 = ForwardDiff.jacobian(ak, x0)
	@show b0 = -ak([0,0])
	@show Q0 = ForwardDiff.hessian(f, x0)
	@show p0 = gk(x0)
	rk = rank([A0' gk(x0)]) == rank(A0')
end

# ╔═╡ 413b064c-e093-4176-a71e-2926b8c0eeff
function bfgs_update(H, s, y)
    rho = 1.0 / (y' * s)
    I = Matrix{Float64}(I, length(s), length(s))
    V = I - rho * s * y'
    H = V' * H * V + rho * s * s'
    return H
end

# ╔═╡ 093070cc-0966-4d3c-a726-b09458851d25
function sqp(obj_func, gradient, x0; max_iter=100, tol=1e-6)
    x = x0
    H = Matrix{Float64}(I, length(x0), length(x0))  # Initial Hessian approximation
    for iter in 1:max_iter
        grad = gradient(x)
        p = -H \ grad  # Search direction
        alpha = 1.0  # Step length (can be optimized further)
        x_new = x + alpha * p
        s = x_new - x
        y = gradient(x_new) - grad
        H = bfgs_update(H, s, y)
        x = x_new
        if norm(grad) < tol
            break
        end
    end
    return x
end

# ╔═╡ db344806-d707-45ac-bdf2-29f98b45ed08
begin
	fig1 = Figure(size=(600,400), fontsize=20)
	ax1 = Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", aspect=1.4)
	#limits!(ax1, -1.5, 1.5, -0.25, 1.5)

	xs1 = -4:0.02:4; 
	fp = (x1, x2) -> f([x1, x2])
	a1p = (x1, x2) -> a1([x1, x2])
	a2p = (x1, x2) -> a2([x1, x2])
	contour!(ax1, xs1, xs1, fp, levels=-20:0.5:20)
	contour!(ax1, xs1, xs1, a1p, levels=[0], linewidth=2, color=:red)
	contour!(ax1, xs1, xs1, a2p, levels=[0], linewidth=2, color=:blue)

	xb = -0.8:0.1:0.8
	g1b = (x) -> x^2 
	g2b = (x) -> sqrt(1-x^2)
	#band!(ax1, xb, g2b.(xb), g1b.(xb), color=(:orange, 0.4))
	
	scatterlines!(ax1, xp[1,:], xp[2,:], markersize=15, color=:blue, linestyle=:dash)
	scatter!(ax1, 1/sqrt(2), 1/sqrt(2), markersize=12, color=:red, 
		strokecolor=:black, strokewidth=1)
	text!(ax1, 0.52, 0.6, text=L"\mathbf{x}^\ast")
	text!(ax1, xp[1,1], xp[2,1], text=L"\mathbf{x}_0")
	fig1
end

# ╔═╡ 3db54003-59b1-4e72-ad41-922d8f6fb0f8
begin
#= 
Title: QR-decomposition based algorithm for convex QP problems with equality 
		constraints.
Description: Implements a QR-decomposition algorithm for convex QP problems with 
		inequality constraints.
Theory: See Practical Optimization Sec. 13.2
Input: 
	H -- positive semidefinite Hessian matrix
	A -- full row-rank constraint matrix 
	p, b -- input vectors
Output:
	xs -- solution vector 
Example:
% Find a minimizer of the QP problem
% minimize 0.5*x'*H*x +x'*p
% subject to A*x = b
% where
% H = [1 0 0; 0 1 0; 0 0 0.01]
% p = [2 1 -1]'
% A = [0 1 1]
% b = 1
% Solution:
% Execute the following commands:
% H = [1 0 0; 0 1 0; 0 0 0.01]
% p = [2 1 -1]'
% A = [0 1 1]
% b = 1
% xs = qp_e(H,p,A,b)
=#	
	function qp_e(H, p, A, b)
		if minimum(size(A)) == 0
			xs = -H\p;          # no active constraints
		else
			n = length(p);
			p1 = size(A)[1];
			Q, R = qr(A')
			R = R[1:p1, :]
			Q1 = Q[:, 1:p1]
			Q2 = Q[:, p1+1:n]
			xs = Q1 * (inv(R')*b)
			Hh = Q2' * H * Q2
			ph = Q2' * (H * xs + p)
			phi_s = -Hh\ph
			xs = Q2 * phi_s + xs
		end
	end
end

# ╔═╡ ba615e8d-2de5-421f-813e-baa401cf00f3
begin
	H = [1 0 0; 0 1 0; 0 0 0.01];
	p = [2 1 -1]';
	A = [0 1 1];
	b = 1
	xs = qp_e(H, p, A, b)
end

# ╔═╡ dc334185-09b7-421a-99e5-117006e8aafc
function qp_ie0(H, p, A, b, x0)
    println("Program qp_ie0.jl")
	n = size(x0,1); xs = x0
    ∇f = x -> H*x + p
	k = 1
	xs = x0
	#-------------------------------------------------------------------
	# Step 0
	rd = A*xs - b
	Aidx = findall(x -> x == 0, rd)
	Aa = A[Aidx,:] 					# active constraints
	ba = b[Aidx]
	y = [H Aa'; Aa zeros(size(Aa,1), size(Aa',2))]\[-p ; ba]
	# Check Feasible
	rd = A*y[1:n] - b
	xf = findall(x -> x > 0, rd)
	if isempty(xf) 
		println("y is feasible.")
		xs = y[1:n]
		Aidx = findall(x -> x .== 0, rd)
		μ = Aa'\(-∇f(xs))
		Sidx = findall(x-> x < 0, μ)
		if !isempty(Sidx)
			deleteat!(Aidx, Sidx)
		end
	else 
		println("y is infeasible.")
		d = y[1:n] - xs
		Ad = A*d
		αidx = findall(x -> x .> 0, Ad)
		α, max_idx = findmax((b[αidx] - A[αidx,:]*xs) ./ (A[αidx,:]*d))
		
		xs = xs + α*d
		sort!(push!(Aidx, max_idx))
		
	end
	
	#-----------------------------------------------------------------
	# Step 1
	Aa = A[Aidx,:]
	ba = b[Aidx]
	y = [H Aa'; Aa zeros(size(Aa,1), size(Aa',2))]\[-p ; ba]
	# Check Feasible
	rd = A*y[1:n] - b
	xf = findall(x -> x > 0, rd)
	if isempty(xf) 
		println("y is feasible.")
		xs = y[1:n]
		Aidx = findall(x -> x .== 0, rd)
		μ = Aa'\(-∇f(xs))
		Sidx = findall(x-> x < 0, μ)
		if !isempty(Sidx)
			deleteat!(Aidx, Sidx)
		end
	else 
		println("y is infeasible.")
		d = y[1:n] - xs
		Ad = A*d
		αidx = findall(x -> x .> 0, Ad)
		α, max_idx = findmax((b[αidx] - A[αidx,:]*xs) ./ (A[αidx,:]*d))
		
		xs = xs + α*d
		sort!(push!(Aidx, max_idx))
		
	end
	k += 1
	#-----------------------------------------------------------------
	# Step 2
	Aa = A[Aidx,:]
	ba = b[Aidx]
	y = [H Aa'; Aa zeros(size(Aa,1), size(Aa',2))]\[-p ; ba]
	# Check Feasible
	rd = A*y[1:n] - b
	xf = findall(x -> x > 0, rd)
	if isempty(xf) 
		println("y is feasible.")
		xs = y[1:n]
		Aidx = findall(x -> x .== 0, rd)
		μ = Aa'\(-∇f(xs))
		Sidx = findall(x-> x < 0, μ)
		deleteat!(Aidx, Sidx)
	else 
		println("y is infeasible.")
		d = y[1:n] - xs
		Ad = A*d
		αidx = findall(x -> x .> 0, Ad)
		α, max_idx = findmax((b[αidx] - A[αidx,:]*xs) ./ (A[αidx,:]*d))
		
		xs = xs + α*d
		sort!(push!(Aidx, max_idx))
		
		#@show sort!(push!(Aidx, αidx))
	end
	k += 1
	#-----------------------------------------------------------------
	# Step 3
	Aa = A[Aidx,:]
	ba = b[Aidx]
	y = [H Aa'; Aa zeros(size(Aa,1), size(Aa',2))]\[-p ; ba]
	# Check Feasible
	rd = A*y[1:n] - b
	xf = findall(x -> x > 0, rd)
	if isempty(xf) 
		println("y is feasible.")
		xs = y[1:n]
		Aidx = findall(x -> x .== 0, rd)
		μ = Aa'\(-∇f(xs))
		Sidx = findall(x-> x < 0, μ)
		deleteat!(Aidx, Sidx)
	else 
		println("y is infeasible.")
		d = y[1:n] - xs
		Ad = A*d
		αidx = findall(x -> x .> 0, Ad)
		α, max_idx = findmax((b[αidx] - A[αidx,:]*xs) ./ (A[αidx,:]*d))
		
		xs = xs + α*d
		sort!(push!(Aidx, max_idx))
	end
	k += 1
	#------------------------------------------------------------------
   
    fs = 0.5 * xs' * (H * xs + 2 * p)
    
    return xs, fs, k-1
end

# ╔═╡ 0814a569-3265-488d-9b0e-6b88aa57b977
begin
	# Example usage
	H2 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.01]
	p2 = [2.0, 1.0, -1.0]
	A2 = [0.0 1.0 1.0]
	b2 = [1.0]
	x02 = [1.0, 1.0, 1.0]
	#----------------------------
	H3 = [2 -1; -1 2]
	p3 = [-3, 0]
	A3 = [1 1; -1 0; 0 -1; 1 0]
	b3 = [2, 0, 0, 3/2]
	x03 = [0, 0]
	xs2, fs2, k = qp_ie0(H3, p3, A0, b0, x03)
	#xs2, fs2, k = qp_ie0(H2, p2, A2, b2, x02)
	println("Solution: ", xs2)
	println("Objective function value: ", fs2)
	println("Iterations: ", k)
end

# ╔═╡ fd0ca00c-4731-497c-83e5-a9baa0b156cb
2/3

# ╔═╡ cbe44268-5a49-43fb-ae01-4c5e51d735a1


# ╔═╡ cceda2a3-1528-40ed-a6d3-ec59792a6671
begin
	alx = []
	findmin(alx)
end

# ╔═╡ Cell order:
# ╠═d42f2a4e-e0ac-4eb3-9b24-cf66f0df4bd6
# ╠═acd5768b-185f-46e5-8d15-2694809d0a90
# ╠═fd9add82-4ab2-4ef7-a2d6-b14e978e77d2
# ╠═250253c0-eb66-4114-b0fd-ea4c64b5d5c1
# ╠═093070cc-0966-4d3c-a726-b09458851d25
# ╠═413b064c-e093-4176-a71e-2926b8c0eeff
# ╠═db344806-d707-45ac-bdf2-29f98b45ed08
# ╠═3db54003-59b1-4e72-ad41-922d8f6fb0f8
# ╠═ba615e8d-2de5-421f-813e-baa401cf00f3
# ╠═dc334185-09b7-421a-99e5-117006e8aafc
# ╠═0814a569-3265-488d-9b0e-6b88aa57b977
# ╠═fd0ca00c-4731-497c-83e5-a9baa0b156cb
# ╠═cbe44268-5a49-43fb-ae01-4c5e51d735a1
# ╠═cceda2a3-1528-40ed-a6d3-ec59792a6671
