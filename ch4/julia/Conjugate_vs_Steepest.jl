### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
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

# ╔═╡ 607f325f-77ee-4421-ba31-8806e161d2e3
md"""
Quadratic function:

vertical axis β₁ = $(@bind β₁ PlutoUI.Slider(1:1:20; show_value=true, default=3)), horizontal β₂  $(@bind β₂ PlutoUI.Slider(1:1:20; show_value=true, default=3))
"""

# ╔═╡ 81cd9137-945f-4aa9-aa4a-82fdca549b61
begin
	# quadratic function
	Quadratic(x) = β₁*x[1]^2 + β₂*x[2]^2;
	
	# bean function
	Bean(x; a=1, b=0.5) = (a-x[1])^2 + (1-x[2])^2 + (b)*(2x[2] - x[1]^2)^2; nothing

	# Rosenbrock function
	Rosenbrock(x; a=1, b=0.5) = (a - x[1])^2 + b*(x[2] -x[1]^2)^2;
end

# ╔═╡ a9a572e6-2d50-4d24-8877-81d9647a95c1
begin
	md"""
	x = $(@bind x01 PlutoUI.Slider(-2:0.1:2, show_value=true, default = -1.2)), $(@bind x02 PlutoUI.Slider(-2:0.1:2, show_value=true, default=2))
	"""
end

# ╔═╡ 60894cf7-333b-4db5-8735-ea1019ba729f
md"""
### Objective function
"""

# ╔═╡ e6f40080-d609-47a9-9339-dcf6ae48c8a4
@bind f Select([Quadratic, Bean, Rosenbrock])

# ╔═╡ dda03848-5c8f-4bf9-9bcc-3fa27cdc0ab2
begin
	#∇f(x) = ForwardDiff.gradient(f,x)
	#pf(x,y) = f([x,y])
	
	∇f(x) = ForwardDiff.gradient(f,x)		
	pf(x,y) = f([x, y]); nothing
end

# ╔═╡ 60dd22eb-e428-4d08-804d-90ede695295d
begin
	@variables x₁ x₂

	ft = latexstring("f(x_1, x_2) = $(pf(x₁, x₂))")
end

# ╔═╡ 5db03d7c-a7b5-4280-aeb6-dcf1a023b9b9
function bisection(f, a₀, b₀, ϵ)

    function D(f, a)
        # Approximate the first derivative using central differences
        h = 0.001
        return (f(a + h) - f(a - h)) / (2 * h)
    end

    a = a₀
    b = b₀

    while (b - a) > ϵ
        c = (a + b) / 2.0

        if D(f, c) > 0
            b = c
        else
            a = c
        end
    end

    return (a + b) / 2 # this was changed
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
	res = Optim.optimize(objective, a, b)   # minimize
	α = Optim.minimizer(res)
	#α = bisection(objective, a, b, 1e-3)
	
	return α
end

# ╔═╡ 625cdb29-8f8d-45f6-b4e7-ff1613000295
md"""
### Backtracking line search
"""

# ╔═╡ 8f75d04c-0e57-4d26-9369-beef72b718be
function backtracking_line_search(f, ∇f, x, d, α; p=0.5, β=1e-4)
	y, g = f(x), ∇f(x) 
	while f(x + α*d) > y + β*α*(g'*d)
		α *= p
	end
	
	return α
end

# ╔═╡ e3c73ba0-048d-49e4-9d9c-9903873c990e
begin
	# Commond use variable
	N = 500 			# maximum iterations
	x0 = [x01, x02] 	# initial point
	ε_G = 1e-4; ε_A = 1e-14; ε_R = 1e-14
	# Mesh Grid
	xxs1 = -5:0.05:4;
	xxs2 = -5:0.05:4; nothing
end

# ╔═╡ 8af0fbd0-8cb0-4c25-a1f6-c4cf3d1abab7
md"""
### Construct the Gradient Descent
"""

# ╔═╡ 10039d6d-edc9-4eda-a2d3-d5e9d681a6f4
begin
	xgra = zeros(2, N)
	xgra[:, 1] = x0;  nothing
end

# ╔═╡ 4ecde5b8-aa69-4197-a06f-dd1bd67c4cba
md"""
### Construct the Conjugate Gradient
"""

# ╔═╡ 3eac4bfa-926d-499c-bd71-0205e20772a0
begin
	g0 = ∇f(x0) 							# Gradient of bean function
	d0 = g0  								# Steepest descent	
	xcon = zeros(2, N)						# a vector of solution 
	xcon[:,1] = x0; nothing
end

# ╔═╡ c9bffae7-d262-4c51-968a-891a2e5be02d
md"""
### Gradient Descent
"""

# ╔═╡ 3b119295-6694-419d-a792-7736ac39f795
begin
	abstract type DescentMethod end
	mutable struct GradientDescent <: DescentMethod
		α
	end

	function init!(M::GradientDescent, f, ∇f, x)
		return M
	end

	function step!(M::GradientDescent, f, ∇f, x)
		g = ∇f(x)
		#α =  backtracking_line_search(f, ∇f, x, -g, 10; p=0.5, β=1e-4)
		α =  line_search(f, x, -g)
		return x - α*g
	end
end

# ╔═╡ 191ab048-6dff-421e-9ee1-cbcb0c4169f7
md"""
### Conjugate Gradient Method
"""

# ╔═╡ a797e874-95e1-479c-8cee-ce3dd968106a
begin
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
		d′ = -g′ + β*d 
		#α =  backtracking_line_search(f, ∇f, x, -g′, 0.50; p=0.5, β=1e-4)
		α = line_search(f, x, d′)
		x′ = x + α*d′
		M.d, M.g = d′, g′

		return x′
	end
		
end

# ╔═╡ 45eb692d-11c4-4bbb-9b5c-490f1607538c
begin
	xrgra = [];
	Mg = GradientDescent(0) 			# initial α is 0
	Mg = init!(Mg, f, ∇f, x0)
	
	# Next Step
	for i = 2:N
		xgra[:,i] = step!(Mg, f, ∇f, xgra[:,i-1])
				
		#if norm(f(xgra[:,i]) .- f(xgra[:,i-1])) <= ε_A + ε_R*f(xcon[:,i-1])
		if norm(∇f(xgra[:,i])) <= ε_G
			global xrgra = xgra[:,1:i]	
			break;
		else
			global xrgra = xgra[:,1:i-1]
		end
	end
end

# ╔═╡ 30aae91c-16e5-4840-af61-e6401517599d
begin
	# Quadratic with degree 2, the CG can converse in two steps.
	xrcon = [];
	# Conjugate Gradient initial step 
	Mc = ConjugateGradientDescent(g0, d0)
	Mc = init!(Mc, f, ∇f, x0)

	# Next Step
	for i = 2:N
		xcon[:,i] = step!(Mc, f, ∇f, xcon[:,i-1])	
		
		#if norm(f(xcon[:,i]) .- f(xcon[:,i-1])) < ε_A + ε_R*f(xcon[:,i])
		if norm(∇f(xcon[:,i])) <= ε_G	
			global xrcon = xcon[:,1:i]	
			break;
		else
			global xrcon = xcon[:,1:i-1]
		end
	end
end

# ╔═╡ 86974899-c880-46fb-8874-6349f4fd4e57
begin
	fig1 = Figure(size=(1000,400))
	#--------------------------------------------------------------------------
	# Steepest Descent
	bx2 = CairoMakie.Axis(fig1[1,1], xlabel= L"x_1", ylabel = L"x_2",
		aspect = AxisAspect(1.3), backgroundcolor=(:blue, 0.05))
	limits!(bx2, -2.2, 2.2, -2.2, 3.2)

	# Contour and initial point
	text!(bx2, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	scatter!(bx2, x0[1], x0[2], markersize=15, strokewidth=2, strokecolor=:blue, 
		color=:lightblue)
	lv1 = -0:5:200 				# level curve
	contour!(bx2, xxs1, xxs2, pf, levels = lv1, color=(:blue, 0.4), linewidth=2)
	contour!(bx2, xxs1, xxs2, pf, levels = [2, 0.3, 0.11], color=(:blue, 0.4), 
		linewidth=2)

	# Scatter line and the optimal point.
	Np1 = size(xrgra, 2)
	scatterlines!(bx2, xrgra[1,1:Np1-1], xrgra[2,1:Np1-1], color=:blue, 
		markercolor=:lightblue, markersize=15, strokewidth = 1, strokecolor=:blue, linewidth=2, label=("SD with $(Np1-2) iterations"))
	scatter!(bx2, xrgra[1,end], xrgra[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(bx2, xrgra[1,end]-0.1, xrgra[2,end]+0.1, text=L"x^\ast", fontsize=22)

	#------------------------------------------------------------------------
	# Conjugate gradient
	bx1 = CairoMakie.Axis(fig1[1,2], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.3), backgroundcolor=(:blue, 0.05))
	limits!(bx1, -2.2, 2.2, -2.2, 3.2)

	# Contour and initial point
	text!(bx1, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	scatter!(bx1, x0[1], x0[2], markersize=15, strokewidth=2, strokecolor=:blue, 
		color=:lightblue)
	lv1 = -0:5:200 				# level curve
	contour!(bx1, xxs1, xxs2, pf, levels = lv1,  color=(:blue, 0.4), linewidth=2)
	contour!(bx1, xxs1, xxs2, pf, levels = [2, 0.3, 0.11], color=(:blue, 0.4), 
		linewidth=2)

	# Scatter line and the optimal point.
	Np = size(xrcon, 2)
	scatterlines!(bx1, xrcon[1,1:Np-1], xrcon[2,1:Np-1], color=:blue, 
		markercolor=:lightblue, markersize=15, strokewidth = 1, strokecolor=:blue, linewidth=2, label=("CG with $(Np-2) iterations"))
	scatter!(bx1, xrcon[1,end], xrcon[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(bx1, xrcon[1,end]-0.1, xrcon[2,end]+0.1, text=L"x^\ast", fontsize=22)
	
	axislegend(bx2, fontsize=18)
	axislegend(bx1, fontsize=18)
	fig1
end

# ╔═╡ Cell order:
# ╠═6ea97aac-74b5-11ef-10e1-fdeed32ca35f
# ╠═9578b091-85b5-4185-a318-2d57ce65be17
# ╟─0076ca61-224c-4867-a508-c7ac9df7abd5
# ╠═f06e8f8f-33c7-4ed8-836d-6de47a88582d
# ╠═81cd9137-945f-4aa9-aa4a-82fdca549b61
# ╠═dda03848-5c8f-4bf9-9bcc-3fa27cdc0ab2
# ╟─607f325f-77ee-4421-ba31-8806e161d2e3
# ╟─a9a572e6-2d50-4d24-8877-81d9647a95c1
# ╠═60894cf7-333b-4db5-8735-ea1019ba729f
# ╟─60dd22eb-e428-4d08-804d-90ede695295d
# ╠═e6f40080-d609-47a9-9339-dcf6ae48c8a4
# ╠═86974899-c880-46fb-8874-6349f4fd4e57
# ╠═5db03d7c-a7b5-4280-aeb6-dcf1a023b9b9
# ╟─d478eabe-501e-4421-be8b-b7d001dcffe4
# ╠═455b4dc4-fb60-4a1f-8bf7-bb2fd2f40aca
# ╟─625cdb29-8f8d-45f6-b4e7-ff1613000295
# ╠═8f75d04c-0e57-4d26-9369-beef72b718be
# ╠═e3c73ba0-048d-49e4-9d9c-9903873c990e
# ╟─8af0fbd0-8cb0-4c25-a1f6-c4cf3d1abab7
# ╠═10039d6d-edc9-4eda-a2d3-d5e9d681a6f4
# ╠═45eb692d-11c4-4bbb-9b5c-490f1607538c
# ╟─4ecde5b8-aa69-4197-a06f-dd1bd67c4cba
# ╠═3eac4bfa-926d-499c-bd71-0205e20772a0
# ╠═30aae91c-16e5-4840-af61-e6401517599d
# ╟─c9bffae7-d262-4c51-968a-891a2e5be02d
# ╠═3b119295-6694-419d-a792-7736ac39f795
# ╟─191ab048-6dff-421e-9ee1-cbcb0c4169f7
# ╠═a797e874-95e1-479c-8cee-ce3dd968106a
