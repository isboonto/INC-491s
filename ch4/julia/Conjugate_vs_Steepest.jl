### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6ea97aac-74b5-11ef-10e1-fdeed32ca35f
using Pkg; Pkg.activate()

# ╔═╡ 9578b091-85b5-4185-a318-2d57ce65be17
begin
	using CairoMakie
	set_theme!(theme_latexfonts(), fontsize=22)
	using Optim, LinearAlgebra, ForwardDiff
	using Symbolics, LaTeXStrings
	using PlutoUI
end

# ╔═╡ 0076ca61-224c-4867-a508-c7ac9df7abd5
md"""
### Bracket_minimum 
"""

# ╔═╡ f06e8f8f-33c7-4ed8-836d-6de47a88582d
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

# ╔═╡ d478eabe-501e-4421-be8b-b7d001dcffe4
md"""
### Exact Line Search
"""

# ╔═╡ 455b4dc4-fb60-4a1f-8bf7-bb2fd2f40aca
# Exact line search
function line_search(f, x, d)
	objective = α -> f(x + α*d)
	a, b = bracket_minimum(objective)
	
	# using brent method from Optim.jl
	res = Optim.optimize(objective, a, b)     # minimize
	α = Optim.minimizer(res)
	return α
end

# ╔═╡ 191ab048-6dff-421e-9ee1-cbcb0c4169f7
md"""
### Conjugate Gradient Method
"""

# ╔═╡ a797e874-95e1-479c-8cee-ce3dd968106a
begin
	abstract type DescentMethod end
	
	mutable struct ConjugateGradientDescent <: DescentMethod
		d
		g
	end

	function init!(M::ConjugateGradientDescent, f, ∇f, x)
		M.g = ∇f(x)
		M.d = -M.g

		return M
	end
	
	function step!(M::ConjugateGradientDescent, f, ∇f, x)
		d, g = M.d, M.g
		g′ = ∇f(x)		#g\prime <TAB>
		β = max(0, dot(g′, g′ - g)/dot(g,g)) 	# Polok-Ribiere
		#β = dot(g′, g′)/dot(g,g)
		d′ = -g′/norm(g′) + β*d 
		
		α = line_search(f, x, d′)
		x′ = x + α*d′
		M.d, M.g = d′, g′

		return x′
	end
		
end

# ╔═╡ 607f325f-77ee-4421-ba31-8806e161d2e3
md"""
β = $(@bind β PlutoUI.Slider(1:1:20; show_value=true, default=3))
"""

# ╔═╡ dda03848-5c8f-4bf9-9bcc-3fa27cdc0ab2
begin
	# quadratic function
	f(x) = x[1]^2 + β*x[2]^2 
	∇f(x) = ForwardDiff.gradient(f,x)
	pf(x,y) = f([x,y])
	
	# bean function
	fb(x; a=1, b=0.5) = (a-x[1])^2 + (1-x[2])^2 + (b)*(2x[2] - x[1]^2)^2
	∇fb(x) = ForwardDiff.gradient(fb,x)		
	pfb(x,y) = fb([x, y])
end

# ╔═╡ a9a572e6-2d50-4d24-8877-81d9647a95c1
begin
	md"""
	x = $(@bind x01 PlutoUI.Slider(-2:0.1:2, show_value=true, default = -1.2)), $(@bind x02 PlutoUI.Slider(-1:0.1:3, show_value=true, default=2))
	"""
end

# ╔═╡ 60894cf7-333b-4db5-8735-ea1019ba729f
md"""
### Objective function
"""

# ╔═╡ 60dd22eb-e428-4d08-804d-90ede695295d
begin
	@variables x₁ x₂

	ft = latexstring("f(x_1, x_2) = $(pf(x₁, x₂))")
end

# ╔═╡ e3c73ba0-048d-49e4-9d9c-9903873c990e
begin
	# Commond use variable
	N = 500 			# maximum iterations
	x0 = [x01, x02] 	# initial point
	ε_G = 1e-10; ε_A = 1e-10; ε_R = 1e-10
	# Mesh Grid
	xxs1 = -5:0.05:3;
	xxs2 = -5:0.05:3; nothing
end

# ╔═╡ 30aae91c-16e5-4840-af61-e6401517599d
begin
	# Quadratic with degree 2, the CG can converse in two steps.
	
	g0 = ∇fb(x0) 							# Gradient of bean function
	d0 = -g0  								# Steepest descent	
	xcon = zeros(2, N)						# a vector of solution 
	xcon[:,1] = x0
	xrcon = x0
		
	# Conjugate Gradient initial step 
	Mc = ConjugateGradientDescent(g0, d0)
	Mc = init!(Mc, fb, ∇fb, x0)

	# Next Step
	for i = 2:N
		xcon[:,i] = step!(Mc, fb, ∇fb, xcon[:,i-1])	
		
		if norm(fb(xcon[:,i]) .- fb(xcon[:,i-1])) < ε_A + ε_R*norm(fb(xcon[:,i-1]))
			global xrcon = xcon[:,1:i]	
			break;
		end
	end
	
end

# ╔═╡ 86974899-c880-46fb-8874-6349f4fd4e57
begin
	fig1 = Figure(size=(800,400))
	bx1 = CairoMakie.Axis(fig1[1,1], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.3))
	limits!(bx1, -2, 3, -1, 3)
	text!(bx1, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	
	scatter!(bx1, x0[1], x0[2], markersize=15, strokewidth=2, strokecolor=:blue, 
		color=:lightblue)
	lv1 = -0:5:200 				# level curve
	contour!(bx1, xxs1, xxs2, pfb, levels = lv1,  color=(:blue, 0.3), linewidth=2)
	contour!(bx1, xxs1, xxs2, pfb, levels = [2, 0.3, 0.11], color=(:blue, 0.3), 
		linewidth=2)

	Np = size(xrcon, 2)
	scatterlines!(bx1, xrcon[1,1:Np-1], xrcon[2,1:Np-1], color=:blue, 
		markercolor=:lightblue, markersize=15, strokewidth = 1, strokecolor=:blue, linewidth=2)
	scatter!(bx1, xrcon[1,end], xrcon[2,end], markersize=15, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(bx1, xrcon[1,end]-0.1, xrcon[2,end]+0.1, text=L"x^\ast", fontsize=22)
	
	fig1
end

# ╔═╡ Cell order:
# ╠═6ea97aac-74b5-11ef-10e1-fdeed32ca35f
# ╠═9578b091-85b5-4185-a318-2d57ce65be17
# ╟─0076ca61-224c-4867-a508-c7ac9df7abd5
# ╠═f06e8f8f-33c7-4ed8-836d-6de47a88582d
# ╟─d478eabe-501e-4421-be8b-b7d001dcffe4
# ╠═455b4dc4-fb60-4a1f-8bf7-bb2fd2f40aca
# ╟─191ab048-6dff-421e-9ee1-cbcb0c4169f7
# ╠═a797e874-95e1-479c-8cee-ce3dd968106a
# ╠═dda03848-5c8f-4bf9-9bcc-3fa27cdc0ab2
# ╟─607f325f-77ee-4421-ba31-8806e161d2e3
# ╟─a9a572e6-2d50-4d24-8877-81d9647a95c1
# ╟─60894cf7-333b-4db5-8735-ea1019ba729f
# ╟─60dd22eb-e428-4d08-804d-90ede695295d
# ╟─e3c73ba0-048d-49e4-9d9c-9903873c990e
# ╟─30aae91c-16e5-4840-af61-e6401517599d
# ╠═86974899-c880-46fb-8874-6349f4fd4e57
