### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 5a5ddd95-e785-49c6-af4c-b91b239fc027
begin
	using Pkg; Pkg.activate()
end

# ╔═╡ d0477b7a-732f-11ee-104d-271b3cbcae8d
begin
	using LinearAlgebra
	#using Plots; default(fontfamily="Computer Modern", guidefontsize=12,
	#	tickfontsize=9, framestyle=:box) #LaTeX Style
	using CairoMakie
		set_theme!(theme_latexfonts(), fontsize=18)
	using ForwardDiff, Symbolics
	using PlutoUI
end

# ╔═╡ 9f972d7c-8515-47b4-b43c-195d4680e8f5
# from  algorithms for Optimization book pp 36 
# Algorithm 3.1
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
	a, ya = x, f(x)
	b, yb = a + s, f(a + s)
	
	if yb > ya 
		a, b = b, a 		# swap
		ya, yb = yb , ya
		s = -s
	end
	while true
		c, yc = b + s, f(b + s)
		if yc > yb
			return a < c ? (a, c) : (c, a)
		end
		a, ya, b, yb = b, yb, c, yc
		s *= k
	end
end

# ╔═╡ 6131e2ca-4ca9-49f1-86bb-29e4a3137bd5
function backtracking_line_search(f, ∇f, x, d, α; p=0.5, β=1e-4)
	y, g = f(x), ∇f(x) 
	while f(x + α*d) > y + β*α*(g'*d)
		@show α *= p
	end
	
	return α, x + α*d
end

# ╔═╡ cd301a1b-6501-4aca-b214-b35545193764
begin
	f(x) = 0.1x[1]^6 - 1.5x[1]^4 + 5x[1]^2 + 0.1x[2]^4 + 3x[2]^2 - 9x[2] + 0.5x[1]x[2]
	∇f(x) = ForwardDiff.gradient(f,x)
	pf(x,y) = f([x,y])
end

# ╔═╡ d25049c2-9a9e-462f-85c5-fa9a11092ca3
begin
	x1 = -3:0.05:3;
	x2 = 0:0.05:2.5;
end

# ╔═╡ 3815d4f2-5e9f-490d-b758-785f8fdc978c
begin
	xp1 = [-1.25, 1.25]
	d = [4, 0.75]
	αi = 1.2
	α, xf1 = backtracking_line_search(f, ∇f, xp1, d, αi; p = 0.7, β = 1e-4)
end

# ╔═╡ c650b3a0-39e0-4700-9100-7dc37632b485
begin
	f(xp1)
	∇f(xp1)
end

# ╔═╡ 9b39a828-23c2-4e48-ab73-794b9beaf9be
let
	fig1 = Figure(size = (800,400))
	ax1 = Axis(fig1[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = AxisAspect(1.5),
		xlabelsize=24, ylabelsize=24)
	limits!(ax1, -3, 3, 0, 2.5)

	text!(ax1, xp1[1]-0.3, xp1[2]+0.1, text=L"x_0", fontsize=24, color=:red)

	contour!(ax1, x1, x2, pf, levels = -20:1:40, linewidth=2)
	scatterlines!(ax1,[xp1[1],  xf1[1]],[xp1[2],  xf1[2]], linewidth=3, 
		markercolor=:LightBlue, markersize=15, strokecolor=:Blue, strokewidth=2)
	fig1
end

# ╔═╡ b26fcaaf-ddd3-4fa7-a67a-7efdc029cdcd
function example_prob(x::Vector{Float64})
	0.3x[1]^6 - 1.5x[1]^4 + 5x[1]^2 + 0.1x[2]^4 + 3x[2]^2 - 9x[2] + 0.5x[1]*x[2]
end

# ╔═╡ Cell order:
# ╠═5a5ddd95-e785-49c6-af4c-b91b239fc027
# ╠═d0477b7a-732f-11ee-104d-271b3cbcae8d
# ╠═9f972d7c-8515-47b4-b43c-195d4680e8f5
# ╠═6131e2ca-4ca9-49f1-86bb-29e4a3137bd5
# ╠═cd301a1b-6501-4aca-b214-b35545193764
# ╠═d25049c2-9a9e-462f-85c5-fa9a11092ca3
# ╠═c650b3a0-39e0-4700-9100-7dc37632b485
# ╠═3815d4f2-5e9f-490d-b758-785f8fdc978c
# ╠═9b39a828-23c2-4e48-ab73-794b9beaf9be
# ╠═b26fcaaf-ddd3-4fa7-a67a-7efdc029cdcd
