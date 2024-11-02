### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ f241bc24-969f-11ef-0b21-13a99b7fe19f
using Pkg; Pkg.activate()

# ╔═╡ df60ff49-b53c-4e18-847b-5de9fde20c47
begin
	using Optim
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=20)
	using PlutoUI
	import Printf: @sprintf
end

# ╔═╡ a406b379-c8ed-4434-a32b-d5016cb75923
md"""
$f(x) = \sin(x), \text{ and }  g(x) = x^2 - 1 ≤ 0$
"""

# ╔═╡ 8d662522-e87e-4840-af24-5a8af40c9562
begin
	f = x -> sin(x)
	# x^2 <= 1
	L = (x, μ) -> f(x) + μ*(x^2 - 1)
	dom = -4:0.1:4

	x_star = -1
end

# ╔═╡ 6d1961a6-7b2d-49b6-8ff5-f2967cd2ea20
begin
	fig = Figure(size = (600,400))
	ax = CairoMakie.Axis(fig[1,1], xlabel=L"x", ylabel=L"y")
	
	
	for (i, μ) in enumerate(collect(0:8)*0.25)
		g = x-> L(x, μ)
		smu = @sprintf("μ = %.2f", μ)
		lines!(ax, dom, g, color=(:royalblue2, i*0.3), label=smu)
		global x = optimize(g, -4, 4).minimizer
		scatter!(ax, x, g(x), markersize=11, strokecolor=:white, strokewidth=1, 
			color=:red)
		println(x)
		
	end
	limits!(ax, -4,4,-3,10)
	# feasible region x^2 ≤ 1
	lines!(ax, dom, f, color=:black)
	lines!(ax, -1:0.1:1, f, color=(:red, 1.0))
	text!(ax, x_star, f(x_star), text=L"x^\ast", fontsize=20)
	
	fig[1,2] = Legend(fig, ax)
	fig
end

# ╔═╡ 36173118-6a52-4b2b-84a6-8097f519e458
function penalty_method(f, p, x, k_max; ρ=1, γ=2)
	for k in 1:k_max
		x = optimize(x -> f(x) + rho*p(x), x)
		ρ *= γ
		if p(x) == 0
			return x
		end
	end
	return x
end

# ╔═╡ 9d7237a7-0ad9-4ca6-9670-e49c186db5e9
begin
	f1 = x -> x[1] + x[2] + x[1]*x[2]
	h1 = x -> x[1]^2 + x[2]^2 - 1
	P = x -> h1(x)^2

	xs_P = zeros(10,2)
	x1 = x_start = [3, 2.5]
	λ = λ_start = 1.0
	γ = 1.5

	is = collect(1:10)
	for i in is
		obj = x -> f1(x) + λ*P(x)
		global x1 = optimize(obj, x1 + randn(2)/5).minimizer
		xs_P[i,:] = x1
		global λ *= γ
	end
	xs_P
	lines(xs_P)
end

# ╔═╡ Cell order:
# ╠═f241bc24-969f-11ef-0b21-13a99b7fe19f
# ╠═df60ff49-b53c-4e18-847b-5de9fde20c47
# ╠═a406b379-c8ed-4434-a32b-d5016cb75923
# ╠═8d662522-e87e-4840-af24-5a8af40c9562
# ╠═6d1961a6-7b2d-49b6-8ff5-f2967cd2ea20
# ╠═36173118-6a52-4b2b-84a6-8097f519e458
# ╠═9d7237a7-0ad9-4ca6-9670-e49c186db5e9
