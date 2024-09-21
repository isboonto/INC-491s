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

# ╔═╡ 76a55074-772f-11ef-02f5-4da9c17e00ed
using Pkg; Pkg.activate()

# ╔═╡ 0e935cc5-e81e-4cfe-9e2d-c7edd47a2d5c
begin
	using CairoMakie
	set_theme!(theme_latexfonts(), fontsize=22)
	using Optim, LinearAlgebra, ForwardDiff
	using Symbolics, LaTeXStrings
	using PlutoUI
end

# ╔═╡ 04fe0465-8a08-4974-8494-22d362cfc0d9
md"""
### Bracket_minimum
"""

# ╔═╡ b9999cd3-3c12-481b-b2f0-419ba3ab3268
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

# ╔═╡ a6b5d39f-dc54-4a6f-8dff-8d5923931983
md"""
Quadratic function:

vertical axis β₁ = $(@bind β₁ PlutoUI.Slider(1:1:20; show_value=true, default=3)), horizontal β₂  $(@bind β₂ PlutoUI.Slider(1:1:20; show_value=true, default=3))
"""

# ╔═╡ d147e839-18ff-40fe-8a61-15d4fd4209b3
begin
	# quadratic function
	Quadratic(x) = β₁*x[1]^2 + β₂*x[2]^2;
	
	# bean function
	Bean(x; a=1, b=0.5) = (a-x[1])^2 + (1-x[2])^2 + (b)*(2x[2] - x[1]^2)^2; nothing

	# Rosenbrock function
	Rosenbrock(x; a=1, b=0.5) = (a - x[1])^2 + b*(x[2] -x[1]^2)^2;
end

# ╔═╡ e105054f-b798-4b3e-85fe-0fafd4519ffe
begin
	md"""
	x = $(@bind x01 PlutoUI.Slider(-2:0.1:2, show_value=true, default = -1.2)), $(@bind x02 PlutoUI.Slider(-2:0.1:2, show_value=true, default=2))
	"""
end

# ╔═╡ 9c9328e2-78a4-4361-a18c-6b8dec5de817
md"""
### Objective function
"""

# ╔═╡ d3d14f77-86a1-4060-814c-aa27b28e3981
@bind f Select([Quadratic, Bean, Rosenbrock])

# ╔═╡ 1bb09348-4ec9-46ff-83bf-dcbdbff5a1d0
begin
	#∇f(x) = ForwardDiff.gradient(f,x)
	#pf(x,y) = f([x,y])
	
	∇f(x) = ForwardDiff.gradient(f,x)		
	∇2f(x) = ForwardDiff.hessian(f,x)
	pf(x,y) = f([x, y]); nothing
end

# ╔═╡ 944e78b3-29be-4b73-9bf4-250996963029
begin
	@variables x₁ x₂

	ft = latexstring("f(x_1, x_2) = $(pf(x₁, x₂))")
end

# ╔═╡ 12a097ad-b16a-4f78-aaf6-d416099a5d31
md"""
### Steepest Descent Plot
"""

# ╔═╡ 054a40ed-263a-4b9a-a8e2-2849268aa02c
md"""
### Newton's method plot
"""

# ╔═╡ fe1d59f0-2283-4912-a13c-f45dd42bcac8
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

# ╔═╡ 9a46e083-89a9-40f0-88e2-2999aa9e54c1
md"""
### Exact Line Search
"""

# ╔═╡ f500ed5b-d066-49db-8ac7-52fd22a3e62f
# Exact line search
function line_search(f, x, d)
	objective = α -> f(x + α*d)
	a, b = bracket_minimum(objective)
	
	# using brent method from Optim.jl
	#res = Optim.optimize(objective, a, b)   # minimize
	#α = Optim.minimizer(res)
	α = bisection(objective, a, b, 1e-3)
	
	return α
end

# ╔═╡ 9983670f-2b49-47ee-84b1-f4e866988768
md"""
### Backtracking line search
"""

# ╔═╡ 4dbf0207-7d14-444a-9ecd-3e504325ddc6
function backtracking_line_search(f, ∇f, x, d, α; p=0.5, β=1e-4)
	y, g = f(x), ∇f(x) 
	while f(x + α*d) > y + β*α*(g'*d)
		α *= p
	end
	
	return α
end

# ╔═╡ 0b7c5ddf-2348-4d20-8bb8-bf426fcea8f2
begin
	# Commond use variable
	N = 5000 			# maximum iterations
	x0 = [x01, x02] 	# initial point
	ε_G = 1e-4; ε_A = 1e-14; ε_R = 1e-14
	# Mesh Grid
	xxs1 = -5:0.05:4;
	xxs2 = -5:0.05:4; nothing
end

# ╔═╡ a562ffc6-56e0-4039-8618-4c213c1bff71
md"""
### Construct the Gradient Descent
"""

# ╔═╡ 90145150-4bc6-4c4c-861e-e6acff41075b
begin
	xgra = zeros(2, N)
	xgra[:, 1] = x0;  nothing
end

# ╔═╡ 700df667-f398-476c-9616-e8c916251e5e
md"""
### Construct the Conjugate Gradient
"""

# ╔═╡ 5cd82c00-26ef-42cb-ae7e-d58f365a919a
begin
	g0 = ∇f(x0) 							# Gradient of bean function
	d0 = g0  								# Steepest descent	
	xcon = zeros(2, N)						# a vector of solution 
	xcon[:,1] = x0; nothing
end

# ╔═╡ e927823d-f10b-491a-a7ec-87e2b2dbcc4a
md"""
### Construct Newton Method
"""

# ╔═╡ bf6c38b9-eda2-4ae0-8e5e-e1578b32666c
begin
	gnew0 = ∇f(x0) 							# Gradient of bean function
	Hnew0 = ∇2f(x0)  								# Steepest descent	
	xnew = zeros(2, N)						# a vector of solution 
	xnew[:,1] = x0; 
	xnewl = zeros(2, N)
	xnewl[:,1] = x0; 
	xnewapl = zeros(2, N)
	xnewapl[:,1] = x0; nothing
end

# ╔═╡ f2cf1870-cdc2-40b2-8229-b1f47a552508
md"""
#### Newton's method with exact linesearch
"""

# ╔═╡ 30c7a4ba-fc91-4890-8316-91398f37e0a1
md"""
#### Newton's method with approximate linesearch
"""

# ╔═╡ 20476b16-8936-4c0f-96fa-b4ed72c7cb4e
md"""
### Descent Method
"""

# ╔═╡ 9eb3739c-37a1-4063-9493-f4fa73fc8294
abstract type DescentMethod end

# ╔═╡ 658984bb-5e9b-4589-9329-1a703ab3cca1
md"""
### Gradient Descent
"""

# ╔═╡ fb150ab1-7723-4c89-8683-1eac024c2361
begin
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

# ╔═╡ 6827c066-028a-4f41-8642-40598c72c566
md"""
### Conjugate Gradient Method
"""

# ╔═╡ 1a4ddc9c-317c-4e0c-8b04-3cf83978ef1b
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

# ╔═╡ 4915e1ae-ac36-4f69-bd13-de21d5c89088
md"""
### Newton's Method
"""

# ╔═╡ 3f42e662-db76-474b-90c4-fe8e5c252890
begin
	mutable struct Newton <: DescentMethod
		g
		H
	end

	function init!(M::Newton, f, ∇f, ∇2f, x)
		M.g = ∇f(x)
		M.H = ∇2f(x)
		
		return M
	end

	# for multivariable function
	function step!(M::Newton, f, ∇f, ∇2f, x)
		g, H = M.g, M.H
		s = -H\g
		x′ = x + s
		M.g = ∇f(x′)
		M.H = ∇2f(x′)
		return x′
	end

	function step_ls!(M::Newton, f, ∇f, ∇2f, x)
		g, H = M.g, M.H
		s = -H\g
		α = line_search(f, x, s)
		x′ = x + α*s
		M.g = ∇f(x′)
		M.H = ∇2f(x′)

		return x′
	end

	function step_apls!(M::Newton, f, ∇f, ∇2f, x)
		g, H = M.g, M.H
		s = -H\g
		α = backtracking_line_search(f, ∇f, x, s, 1.0; p=0.4, β=1e-4)
		x′ = x + α*s
		M.g = ∇f(x′)
		M.H = ∇2f(x′)

		return x′
	end
	
end

# ╔═╡ b0f44627-20ba-4758-b089-5c02dd33fd38
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

# ╔═╡ adeb3f6d-72b3-42b6-93b6-e53f6b17da56
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

# ╔═╡ 3774f258-6fdb-4a1a-acef-efba80f937ee
begin
	fig1 = Figure(size=(1000,800))
	empty!(fig1)
	#--------------------------------------------------------------------------
	# Steepest Descent
	bx2 = CairoMakie.Axis(fig1[1,1], xlabel= L"x_1", ylabel = L"x_2",
		aspect = AxisAspect(1.3), backgroundcolor=(:blue, 0.01))
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
		aspect = AxisAspect(1.3), backgroundcolor=(:blue, 0.01))
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
	
end

# ╔═╡ efc966ba-6fee-48a9-ba64-4aff75ea6bc5
begin
	xrnew = [];
	# Newton initial step
	Mn = Newton(gnew0, Hnew0)
	Mn = init!(Mn, f, ∇f, ∇2f, x0)

	# Next Step
	for i = 2:N
		xnew[:,i] = step!(Mn, f, ∇f, ∇2f, xnew[:,i-1])

		if norm(∇f(xnew[:,i])) <= ε_G
			global xrnew = xnew[:,1:i]
			break;
		else
			global xrnew = xnew[:, 1:i-1]
		end
	end
	
end

# ╔═╡ f57ba092-5e70-44cb-9852-faeee9a48cc2
begin
	xrnewl = [];
	# Newton initial step
	Mnl = Newton(gnew0, Hnew0)
	Mnl = init!(Mnl, f, ∇f, ∇2f, x0)

	# Next Step
	for i = 2:N
		xnewl[:,i] = step_ls!(Mnl, f, ∇f, ∇2f, xnewl[:,i-1])

		if norm(∇f(xnewl[:,i])) <= ε_G
			global xrnewl = xnewl[:,1:i]
			break;
		else
			global xrnewl = xnewl[:, 1:i-1]
		end
	end
	
end

# ╔═╡ 3e2cf0b8-94d7-47c6-b843-badc161df75f
begin
	xrnewapl = [];
	# Newton initial step
	Mnapl = Newton(gnew0, Hnew0)
	Mnapl = init!(Mnapl, f, ∇f, ∇2f, x0)

	# Next Step
	# approximate line search sometime is quite bad. 
	for i = 2:N
		xnewapl[:,i] = step_apls!(Mnapl, f, ∇f, ∇2f, xnewapl[:,i-1])

		if norm(∇f(xnewapl[:,i])) <= ε_G
			global xrnewapl = xnewapl[:,1:i]
			break;
		else
			global xrnewapl = xnewapl[:, 1:i-1]
		end
	end
	
end

# ╔═╡ d141a0dd-460b-43b0-9dbe-e80a1ee6a7b0
begin
	# Newton's method without linesearch
	bx3 = CairoMakie.Axis(fig1[2,1], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.3), backgroundcolor=(:blue, 0.01))
	limits!(bx3, -2.2, 2.2, -2.2, 3.2)

	# Contour and initial point
	text!(bx3, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	scatter!(bx3, x0[1], x0[2], markersize=15, strokewidth=2, strokecolor=:blue, 
		color=:lightblue)
	
	contour!(bx3, xxs1, xxs2, pf, levels = lv1,  color=(:blue, 0.4), linewidth=2)
	contour!(bx3, xxs1, xxs2, pf, levels = [2, 0.3, 0.11], color=(:blue, 0.4), 
		linewidth=2)
	
	# Scatter line and the optimal point
	Nn = size(xrnew, 2)
	
	scatterlines!(bx3, xrnew[1,1:Nn], xrnew[2,1:Nn], color=:blue, 
		markercolor=:lightblue, markersize=15, strokewidth = 1, strokecolor=:blue, linewidth=2, label=("Newton with $(Nn-1) iterations"))
	scatter!(bx3, xrnew[1,end], xrnew[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(bx3, xrnew[1,end]-0.1, xrnew[2,end]+0.1, text=L"x^\ast", fontsize=22)

	#-----------------------------------------------------------------
	# Newton's with linesearch
	bx4 = CairoMakie.Axis(fig1[2,2], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.3), backgroundcolor=(:blue, 0.01))
	limits!(bx4, -2.2, 2.2, -2.2, 3.2)
	
	# Contour and initial point
	text!(bx4, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	scatter!(bx4, x0[1], x0[2], markersize=15, strokewidth=2, strokecolor=:blue, 
		color=:lightblue)
	
	contour!(bx4, xxs1, xxs2, pf, levels = lv1,  color=(:blue, 0.4), linewidth=2)
	contour!(bx4, xxs1, xxs2, pf, levels = [2, 0.3, 0.11], color=(:blue, 0.4), 
		linewidth=2)
	
	# Scatter line and the optimal point
	Nnl = size(xrnewl, 2)
	Nnapl = size(xrnewapl, 2)
	
	# exact line search
	scatterlines!(bx4, xrnewl[1,1:Nnl-1], xrnewl[2,1:Nnl-1], color=:blue, 
		markercolor=:lightblue, markersize=15, strokewidth = 1, strokecolor=:blue, linewidth=2, label=("Newton+ex LS with $(Nnl-1) iterations"))
	
	# approximate line search
	scatterlines!(bx4, xrnewapl[1,1:Nnapl-1], xrnewapl[2,1:Nnapl-1], color=:green, 
		markercolor=:lightgreen, markersize=15, strokewidth = 1, strokecolor=:green, linewidth=2, label=("Newton+app LS with $(Nnapl-1) iterations"))
	
	scatter!(bx4, xrnewl[1,end], xrnewl[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(bx4, xrnewl[1,end]-0.1, xrnewl[2,end]+0.1, text=L"x^\ast", fontsize=22)
	
	
end

# ╔═╡ e7079934-b250-4846-827f-feee8d1ac2c6
begin
	
	axislegend(bx1, fontsize=18)
	axislegend(bx2, fontsize=18)
	axislegend(bx3, fontsize=18)
	axislegend(bx4, fontsize=18)

	fig1
end

# ╔═╡ Cell order:
# ╠═76a55074-772f-11ef-02f5-4da9c17e00ed
# ╠═0e935cc5-e81e-4cfe-9e2d-c7edd47a2d5c
# ╟─04fe0465-8a08-4974-8494-22d362cfc0d9
# ╠═b9999cd3-3c12-481b-b2f0-419ba3ab3268
# ╠═d147e839-18ff-40fe-8a61-15d4fd4209b3
# ╠═1bb09348-4ec9-46ff-83bf-dcbdbff5a1d0
# ╟─a6b5d39f-dc54-4a6f-8dff-8d5923931983
# ╟─e105054f-b798-4b3e-85fe-0fafd4519ffe
# ╟─9c9328e2-78a4-4361-a18c-6b8dec5de817
# ╟─944e78b3-29be-4b73-9bf4-250996963029
# ╟─d3d14f77-86a1-4060-814c-aa27b28e3981
# ╠═e7079934-b250-4846-827f-feee8d1ac2c6
# ╟─12a097ad-b16a-4f78-aaf6-d416099a5d31
# ╠═3774f258-6fdb-4a1a-acef-efba80f937ee
# ╟─054a40ed-263a-4b9a-a8e2-2849268aa02c
# ╠═d141a0dd-460b-43b0-9dbe-e80a1ee6a7b0
# ╠═fe1d59f0-2283-4912-a13c-f45dd42bcac8
# ╟─9a46e083-89a9-40f0-88e2-2999aa9e54c1
# ╠═f500ed5b-d066-49db-8ac7-52fd22a3e62f
# ╟─9983670f-2b49-47ee-84b1-f4e866988768
# ╠═4dbf0207-7d14-444a-9ecd-3e504325ddc6
# ╠═0b7c5ddf-2348-4d20-8bb8-bf426fcea8f2
# ╠═a562ffc6-56e0-4039-8618-4c213c1bff71
# ╠═90145150-4bc6-4c4c-861e-e6acff41075b
# ╠═b0f44627-20ba-4758-b089-5c02dd33fd38
# ╟─700df667-f398-476c-9616-e8c916251e5e
# ╠═5cd82c00-26ef-42cb-ae7e-d58f365a919a
# ╠═adeb3f6d-72b3-42b6-93b6-e53f6b17da56
# ╟─e927823d-f10b-491a-a7ec-87e2b2dbcc4a
# ╠═bf6c38b9-eda2-4ae0-8e5e-e1578b32666c
# ╠═efc966ba-6fee-48a9-ba64-4aff75ea6bc5
# ╟─f2cf1870-cdc2-40b2-8229-b1f47a552508
# ╠═f57ba092-5e70-44cb-9852-faeee9a48cc2
# ╟─30c7a4ba-fc91-4890-8316-91398f37e0a1
# ╠═3e2cf0b8-94d7-47c6-b843-badc161df75f
# ╟─20476b16-8936-4c0f-96fa-b4ed72c7cb4e
# ╠═9eb3739c-37a1-4063-9493-f4fa73fc8294
# ╟─658984bb-5e9b-4589-9329-1a703ab3cca1
# ╠═fb150ab1-7723-4c89-8683-1eac024c2361
# ╠═6827c066-028a-4f41-8642-40598c72c566
# ╠═1a4ddc9c-317c-4e0c-8b04-3cf83978ef1b
# ╟─4915e1ae-ac36-4f69-bd13-de21d5c89088
# ╠═3f42e662-db76-474b-90c4-fe8e5c252890
