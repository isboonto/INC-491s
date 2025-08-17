### A Pluto.jl notebook ###
# v0.20.16

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

# ╔═╡ 6d12590c-7b74-11f0-1629-a916ddf947dc
using Pkg; Pkg.activate()

# ╔═╡ ba33c9b6-8724-4b9e-9880-4f089772f481
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using LinearAlgebra, ForwardDiff
	using PlutoUI, LaTeXStrings
end

# ╔═╡ 5e7b2f91-17a9-4793-aff3-b9ed8d8e68e0
using Symbolics

# ╔═╡ f009577f-2265-4bb8-b1d1-84bf468a7baf
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

# ╔═╡ e859c8ff-45b2-4c71-a9d5-5b4e6cec522e
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

# ╔═╡ fd5db6d2-cdd4-4340-862d-0a193a1a2810
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

# ╔═╡ c10e4bb5-6601-4f84-9481-ae30f461b1dd
md"""
### Objectiv function
"""

# ╔═╡ d6ce048a-29a0-4fe9-8d54-07b1c6115c10
md"""
Quadratic function:

vertical axis β₁ = $(@bind β₁ PlutoUI.Slider(1:1:20; show_value=true, default=3)), horizontal β₂  $(@bind β₂ PlutoUI.Slider(1:1:20; show_value=true, default=3))
"""

# ╔═╡ a711c0f4-823b-45f5-a074-7a7e4218b7a4
begin
	# quadratic function
	Quadratic(x) = β₁*x[1]^2 + β₂*x[2]^2 + x[1]*x[2]

	# Bean function
	Bean(x; a=1, b=0.5) = (a-x[1])^2 + (1-x[2])^2 + (b)*(2x[2] - x[1]^2)^2; nothing

	# Rosenbrock function
	Rosenbrock(x; a=1, b=0.5) = (a - x[1])^2 + b*(x[2] -x[1]^2)^2;
end

# ╔═╡ 54ffce70-0ac7-4783-a61f-acf577750def
@bind f Select([Quadratic, Bean, Rosenbrock])

# ╔═╡ 7514c095-cab9-484b-8c1d-171f067fffd0
begin
	∇f(x) = ForwardDiff.gradient(f, x)
	pf(x, y) = f([x, y]); nothing
end

# ╔═╡ 43682112-7a32-48bc-9d57-d7b942bed5c4
begin
	@variables x₁ x₂

	ft = latexstring("f(x_1, x_2)) = $(pf(x₁, x₂))")
end

# ╔═╡ 2b75caea-40c8-46fd-b484-35dd12e37e10
begin
	fig1 = Figure(size=(800,700)); nothing
end

# ╔═╡ d38a6184-9470-4927-a256-cafb2ae562fa
begin
	md"""
	x = $(@bind x01 PlutoUI.Slider(-2:0.1:2, show_value=true, default = -1.2)), $(@bind x02 PlutoUI.Slider(-2:0.1:2, show_value=true, default=2))
	
	α3 = $(@bind α3 PlutoUI.Slider([0.01, 0.05, 0.1, 0.25, 0.3, 0.5], show_value=true, default= 0.25))
	"""
end

# ╔═╡ 24d72df5-3cc0-4475-b36b-3aad1c693ee0
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/fixed_step_size.pdf", fig1)

# ╔═╡ 7f2952a5-1ff4-43ec-88b8-bcfe84cb4936
md"""
### Construct the Gradient Descent
"""

# ╔═╡ 61097e5a-095d-4e2c-8c03-50a15a6e2bbb
begin
	# commond use variables
	N = 500 			# maximum iterations
	x0 = [x01, x02] 	# initial point
	ε_G = 1e-4; ε_A = 1e-14; ε_R = 1e-14;

	# Mesh Grid
	xxs1 = -5:0.05:4;
	xxs2 = -5:0.05:4; nothing
end

# ╔═╡ 55828d5d-2cc2-4bfd-b179-6a8536dbdf0a
function plot_gd(α, row, col, xrgra, color::Symbol)
	ax1 = CairoMakie.Axis(fig1[row, col], xlabel=L"x_1", ylabel=L"x_2",
						 aspect = AxisAspect(1.3), backgroundcolor=(:cyan, 0.05))
	limits!(ax1, -2.2, 2.2, -2.2, 3.2)

	# Contour and initial point
	lv1 = -0:5:200
	contour!(ax1, xxs1, xxs2, pf, levels = lv1, linewidth=1)
	contour!(ax1, xxs1, xxs2, pf, levels = [2, 0.6, 0.11], linewidth=1)

	text!(ax1, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	scatter!(ax1, x0[1], x0[2], markersize=15, strokewidth=1, strokecolor=:blue,
			color=:lightblue)
	# Scatter line and the optimal point.
	Np1 = size(xrgra, 2)
	scatterlines!(ax1, xrgra[1,1:Np1], xrgra[2, 1:Np1], color=:blue,
				  markercolor=(color, 0.5), markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=2, label=("α = $(α) with $(Np1-2) iterations"))
	scatter!(ax1, xrgra[1,end], xrgra[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(ax1, xrgra[1,end]-0, xrgra[2,end]+0.15, text=L"x^\ast", fontsize=22)


	axislegend(ax1, labelsize=14)
	fig1
end

# ╔═╡ 94513c71-e7a3-49d3-a8a0-48d14848d8e7
begin
	xgra = zeros(2, N)
	xgra[:, 1] = x0; nothing
end

# ╔═╡ b67ab6c4-0ad7-4884-84a2-269af103de4b
begin
	abstract type DescentMethod end
	struct GradientDescent <: DescentMethod
		α
	end
	
	init!(M::GradientDescent, f, ∇f, x) = M
	
	# Fixed-step
	function step!(M::GradientDescent, f, ∇f, x)
		α, g = M.α, ∇f(x)
		#α =  line_search(f, x, -g)
	
		return x - α*g
	end
end

# ╔═╡ b45feff8-105f-4fe0-b7bd-9403c00752bd
function gradient_descent(xrga, α)
	
	Mg = GradientDescent(α) 	# initial α is 0
	Mg = init!(Mg, f, ∇f, x0)
	
	# Next Step
	for i = 2:N
		xgra[:,i] = step!(Mg, f, ∇f, xgra[:,i-1])
		if norm(∇f(xgra[:,i])) <= ε_G
			global xrgra = xgra[:,1:i]
			break;
		else
			global xrgra = xgra[:, 1:i-1]
		end
	end
	return xrgra
end

# ╔═╡ b818ee57-4584-49bf-b05c-a30bd21a7610
begin
	α1 = 0.01; α2 = 0.06;# α3 = 0.25;
	xrgra1 = zeros(length(x0), N);
	xrgra1 = gradient_descent(xrgra1, α1)

	xrgra2 = zeros(length(x0), N);
	xrgra2 = gradient_descent(xrgra2, α2)
	
	xrgra3 = zeros(length(x0), N);
	xrgra3 = gradient_descent(xrgra3, α3)
end

# ╔═╡ 04be0f27-9cfd-4b98-a0a7-71c941ede3c2
function plot_trend(row, col)
	ax1 = CairoMakie.Axis(fig1[row, col], xlabel=L"$$Iteration", ylabel=L"$$O
						  bjective Function Value", aspect = AxisAspect(1.3), backgroundcolor=(:cyan, 0.05))
	limits!(ax1, 0, 30, -1, 15)
	Np1 = size(xrgra1,2)
	z = [f([xrgra1[1,x], xrgra1[2,x]]) for x in 1:Np1]
	scatterlines!(ax1, 1:Np1, z, color=:blue,
				  markercolor=:magenta, markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=1, label=("α = $(α1)"))
	Np1 = size(xrgra2,2)
	z = [f([xrgra2[1,x], xrgra2[2,x]]) for x in 1:Np1]
	scatterlines!(ax1, 1:Np1, z, color=:blue,
				  markercolor=:lightblue, markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=1, label=("α = $(α2)"))
	
	axislegend(ax1, labelsize=14)

	Np1 = size(xrgra3,2)
	z = [f([xrgra3[1,x], xrgra3[2,x]]) for x in 1:Np1]
	scatterlines!(ax1, 1:Np1, z, color=:blue,
				  markercolor=:red, markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=1, label=("α = $(α3)"))
	
	axislegend(ax1, labelsize=14)
	
	fig1
	
end

# ╔═╡ 18783fa7-e31a-4233-8846-064e8b8eb601
begin
	empty!(fig1)
	plot_gd(α1, 1, 1, xrgra1, :magenta)
	plot_gd(α2, 1, 2, xrgra2, :lightblue)
	plot_gd(α3, 2, 1, xrgra3, :red)
	plot_trend(2,2)
end

# ╔═╡ Cell order:
# ╠═6d12590c-7b74-11f0-1629-a916ddf947dc
# ╠═ba33c9b6-8724-4b9e-9880-4f089772f481
# ╠═5e7b2f91-17a9-4793-aff3-b9ed8d8e68e0
# ╠═f009577f-2265-4bb8-b1d1-84bf468a7baf
# ╠═e859c8ff-45b2-4c71-a9d5-5b4e6cec522e
# ╠═fd5db6d2-cdd4-4340-862d-0a193a1a2810
# ╟─c10e4bb5-6601-4f84-9481-ae30f461b1dd
# ╠═a711c0f4-823b-45f5-a074-7a7e4218b7a4
# ╠═7514c095-cab9-484b-8c1d-171f067fffd0
# ╟─d6ce048a-29a0-4fe9-8d54-07b1c6115c10
# ╟─43682112-7a32-48bc-9d57-d7b942bed5c4
# ╠═54ffce70-0ac7-4783-a61f-acf577750def
# ╠═2b75caea-40c8-46fd-b484-35dd12e37e10
# ╟─d38a6184-9470-4927-a256-cafb2ae562fa
# ╠═18783fa7-e31a-4233-8846-064e8b8eb601
# ╠═24d72df5-3cc0-4475-b36b-3aad1c693ee0
# ╟─55828d5d-2cc2-4bfd-b179-6a8536dbdf0a
# ╟─04be0f27-9cfd-4b98-a0a7-71c941ede3c2
# ╠═7f2952a5-1ff4-43ec-88b8-bcfe84cb4936
# ╠═61097e5a-095d-4e2c-8c03-50a15a6e2bbb
# ╠═94513c71-e7a3-49d3-a8a0-48d14848d8e7
# ╠═b45feff8-105f-4fe0-b7bd-9403c00752bd
# ╠═b818ee57-4584-49bf-b05c-a30bd21a7610
# ╠═b67ab6c4-0ad7-4884-84a2-269af103de4b
