### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 8703cb3e-b461-11ef-0bf4-9b990da3ce6d
using Pkg; Pkg.activate()#Pkg.activate("/home/isudonto2000/.julia/")

# ╔═╡ 0f44b0b6-b695-414e-966d-53866494ff19
begin
	using LinearAlgebra, ForwardDiff
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using Symbolics
end

# ╔═╡ 1f88d788-4f49-4258-a7c1-6325ac292949
begin
	f =  (x) -> x[1]^2 - x[1]*x[2] + x[2]^2 - 3x[1]
	#f =  (x) -> -x[1] - x[2]
	a1 = (x) -> x[1] + x[2] - 2
	a2 = (x) -> -x[1]
	a3 = (x) -> -x[2]
	a4 = (x) -> x[1] - 3/2
	#c1 = (x) -> x[1]^2 - x[2]
	#c2 = (x) -> x[1]^2 + x[2]^2 - 1
	ak = x -> [a1(x), a2(x), a3(x), a4(x)]
	#ak = x -> [c1(x), c2(x)]
end

# ╔═╡ 1d33374b-e764-4e38-9a07-88c89b522f67
begin
	@variables x₁ x₂
	
	x = [x₁, x₂]
	ak([x₁, x₂])
	x0 = [0, 0]  				
	gk = x -> ForwardDiff.gradient(f, x)
	A0 = ForwardDiff.jacobian(ak, x0)
	b0 = -ak([0,0])
	Q0 = ForwardDiff.hessian(f, x0)
	p0 = gk(x0)
	rk = rank([A0' gk(x0)]) == rank(A0')
end

# ╔═╡ cbc0a052-8db1-427f-9f40-242b8ce76c18
function qp_ie(H, p, A, b, x0)
	k = 1; N = 1000; n = size(x0,1);
	xs = zeros(n,N)
	xs[:,1] = x0;   
	∇f = x -> H*x + p                       # gradient 
	rd = A*xs[:,1] - b
	A_idx = findall( x -> x == 0, rd)
	#@show A_idx
	while (k <= 10) && (norm(∇f(xs[:,k])) >= 1)
		Aa = A[A_idx,:] 					# active constraints
		ba = b[A_idx]
		y = ([H Aa'; Aa zeros(size(Aa,1), size(Aa',2))])\[-p ; ba]
		rd = A*y[1:n] - b 					# residual vector
		xf = findall(x -> x > 0, rd) 		# feasible x ?
		
		if isempty(xf)
		# if x is feasible find μ using (A_{ak}^T)^{-1}(-∇f(x))
		# Check μ is satisfied with the KKT conditions. μ must not be negative.
		# Delete the constraint that is not satisfied. 
		# x_{k+1} = y
			println("In $k step x is feasible")
			rd = A*y[1:n] - b 	
			A_idx = findall(x -> x == 0, rd)
			μ = Aa'\(-∇f(y[1:n]))
			#@show μ
			μ_idx = findall(x -> x < 0, μ)
			#@show μ_idx
			if !isempty(μ_idx)
				deleteat!(A_idx, μ_idx)
			end
			xs[:,k+1] = y[1:n]
		else 
		# If x is infeasible, find descent direction d_k and also step size α 
		# x_{k+1} = x_k + αd_k
			println("In $k step x is infeasible")
			d = y[1:n] - xs[:,k]
			α_idx = findall(x -> x > 0, A*d)
			α, max_idx = findmax((b[α_idx] - A[α_idx, :]*xs[:,k]) ./ (A[α_idx,:]*d))

			sort!(push!(A_idx, max_idx))
			xs[:,k+1] = xs[:,k] + α*d	
		end
		#@show A_idx
		#@show xs[:,k+1]	
		k += 1 								# dont' forget this
	end
	
	fs = 0.5 * xs[:,k]' * (H * xs[:,k] + 2 * p)
	xso = xs[:,1:k]
	return xso, fs, k
end

# ╔═╡ 949f3baf-4473-4b55-a4b3-9364bcfeb867
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
	#xs2, fs2, k = qp_ie0(H3, p3, A0, b0, x03)
	#xs2, fs2, k = qp_ie0(H2, p2, A2, b2, x02)
	xs2, fs2, k = qp_ie(Q0, p0, A0, b0, x0) 
	println("Solution: ", xs2[:,k])
	println("Objective function value: ", fs2)
	println("Iterations: ", k-1)
end

# ╔═╡ e14283df-77f7-4245-a5ed-cc8b6592ad3f
begin
	fig1 = Figure(size=(600,400), fontsize=18)
	ax1 = Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", aspect=1.4)
	limits!(ax1, -1, 2.5, -1.0, 2.5)

	xs1 = -4:0.02:4;
	fp = (x1, x2) -> f([x1, x2])
	a1p = (x1, x2) -> a1([x1, x2]); a2p = (x1, x2) -> a2([x1, x2]);
	a3p = (x1, x2) -> a3([x1, x2]); a4p = (x1, x2) -> a4([x1, x2]);
	
	contour!(ax1, xs1, xs1, fp, levels = -20:1:100)
	contour!(ax1, xs1, xs1, a1p, levels = [0], linewidth=2, color=:red)
	contour!(ax1, xs1, xs1, a2p, levels = [0], linewidth=2, color=:red)
	contour!(ax1, xs1, xs1, a3p, levels = [0], linewidth=2, color=:red)
	contour!(ax1, xs1, xs1, a4p, levels = [0], linewidth=2, color=:red)

	xb = -0.0:0.1:1.5
	a1b = (x) -> 2 - x
	a2b = (x) -> 0
	band!(ax1, xb, a2b.(xb), a1b.(xb), color=(:orange, 0.4))

	scatterlines!(ax1, xs2, markersize=10, strokewidth=1, strokecolor=:black)
	scatter!(ax1, xs2[1,k], xs2[2,k], color=:green, markersize=10)

	text!(ax1, xs2[1,1]+0.1, xs2[2,1], text=L"\textbf{x}_0")
	text!(ax1, xs2[1,k]+0.1, xs2[2,k], text=L"\textbf{x}^\ast")
	
	fig1
end

# ╔═╡ Cell order:
# ╠═8703cb3e-b461-11ef-0bf4-9b990da3ce6d
# ╠═0f44b0b6-b695-414e-966d-53866494ff19
# ╠═1f88d788-4f49-4258-a7c1-6325ac292949
# ╠═1d33374b-e764-4e38-9a07-88c89b522f67
# ╠═949f3baf-4473-4b55-a4b3-9364bcfeb867
# ╠═e14283df-77f7-4245-a5ed-cc8b6592ad3f
# ╠═cbc0a052-8db1-427f-9f40-242b8ce76c18
