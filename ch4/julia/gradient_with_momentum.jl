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

# ╔═╡ 70563c10-824b-11f0-2115-e93c87dcd75c
using Pkg; Pkg.activate()

# ╔═╡ 2dc11439-74c6-4fde-9f3d-18dfbb8f8712
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using ForwardDiff, LinearAlgebra
	using PlutoUI, LaTeXStrings, Symbolics, Optim, Latexify
end

# ╔═╡ 1dccf6e1-4108-4eac-afd3-8b7ba9288c05
md"""
The Rosenbrock Function

$f(x) = (a - x_1)^2  + b(x_2 - x_1^2)^2,$

where $a = 1$, and $b = 100$.
"""

# ╔═╡ 527deac5-07bd-487c-a8a0-200cbaa031b6
begin
	Rosenbrock(x) = (1 - x[1])^2 + 100*(4*x[2] - x[1]^2)^2
	
	# Boyd and Vanderberghe
	Brunton(x) =  15*exp(-0.25*(x[1] - 2)^2) + 0.5*exp(-2(x[1] - 2)^2) + x[1]^2 + x[2]^2

end

# ╔═╡ 9a78d420-1aaa-452e-b655-ca247e940a76
md"""
**1. Select a function:** $(@bind f Select([Rosenbrock, Brunton]))
"""

# ╔═╡ 98ac464d-f8a0-4b0f-9175-2cd185703b96
if f == Rosenbrock
	# Grid
    xxs1 = -3:0.01:2
    xxs2 = -0.5:0.01:2
	xlim = (-3, 2); ylim = (-0.5, 2)
	level = 0:70:1000
	x0 = [-2, 1.5]
	asp = 1.3
	α = 0.0003
elseif f == Brunton
	xxs1 = -5:0.1:5
    xxs2 = -5:0.1:5
	xlim = (-5, 5); ylim = (-5, 5)
	level = 0:2:50
	x0 = [5, 3]
	asp = 1.3
	α = 0.1
end

# ╔═╡ ae7fda1b-4389-4524-b9b6-6ebf1d4ae754
md"""
**2. Select a calculation:**
"""

# ╔═╡ cacef432-87f7-4977-a602-3f2adfe5b9a2
@variables x₁ x₂

# ╔═╡ 3decb336-9894-40ed-86d3-ef23038827c4
let
	original_str = latexify(f([x₁, x₂]))
	modified_str = "f(x_1, x_2) =" * original_str
	latexstring(modified_str)
end

# ╔═╡ 22b7adf3-7296-44e7-beeb-9a037087ba7b
∇f = x -> ForwardDiff.gradient(f, x)

# ╔═╡ 86af2161-32f3-451e-a65f-ab6e5e4d62ec
let
	original_str = "∇f(x_1, x_2) = " * latexify(∇f([x₁, x₂]))
	latexstring(original_str)
end

# ╔═╡ cd2ef601-9890-44a9-a01d-8ba09a5ad9c7
function bracket_minimum(f, x=0, s=1e-2, k = 2.0)
	# Ensure we start by going downhill
	a, ya = x, f(x)
	b, yb = a + s, f(a + s)

	if yb > ya 			# swap the direction
		a, b = b, a
		ya, yb = yb, ya
		s = -s
	end
	while true
		c, yc = b + s, f(b + s)
		if yc > yb
			return a < c ? (a, c) : (c,a)
		end
		a, ya , b, yb = b, yb, c, yc
		s *= k
	end

end

# ╔═╡ e877de87-3b79-4213-ab2b-73e44c7d6f68
# Exact line search
function line_search(f, x, d)
	objective = α -> f(x + α*d)
	a, b = bracket_minimum(objective)
	
	# using brent method from Optim.jl
	res = Optim.optimize(objective, a, b)   # minimize
	α = Optim.minimizer(res)
		
	return α
end

# ╔═╡ 43009dda-1630-4970-ba5e-e6f85dfd96fd
begin
	fig1 = Figure(size=(800, 600)); nothing
end

# ╔═╡ 98050695-1fd1-4ecf-994d-06c2c8f70279
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/momentum_ex2.pdf", fig1)
	input = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/momentum_ex2.pdf"
	output = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/momentum_ex2.pdf"
	run(`pdfcrop $input $output`)
end
  ╠═╡ =#

# ╔═╡ 03dcf143-8168-40f0-a9ac-fefb23e742e8
begin
# -----------------------------------------------------------------------------------
# Gradient Descent
# -----------------------------------------------------------------------------------
	function gradient_descent_fixed_alpha(f, ∇f, x0; α = 0.1, ε_G = 1e-8, N = 500)
		xgra = zeros(length(x0), N)
		xgra[:, 1] = x0
		for i in 2:N
			d = -∇f(xgra[:, i-1])
			xgra[:, i] = xgra[:, i-1] + α * d
			# We can use line search to improve the gradient method, but it is 
			# very solow
			#α = line_search(f, xgra[:, i-1], d)
			if norm(∇f(xgra[:, i])) ≤ ε_G
				return xgra[:, 1:i]
			end
		end
	
		return xgra
	end
# -----------------------------------------------------------------------------------
# Momentum Method
# -----------------------------------------------------------------------------------
	function gradient_heavy_ball(f, ∇f, x0; α = 0.0003, β = 0.9, ε_G = 1e-8, N = 500)
		xgra = zeros(length(x0), N)
		xgra[:, 1] = x0
		v = zeros(2,1); β = β; α = α
		for i in 2:N
			g = ∇f(xgra[:,i-1])
			v = β*v - α*g
			xgra[:, i] = xgra[:, i-1] + v
			if norm(∇f(xgra[:, i])) ≤ ε_G 
				return xgra[:, 1:i]
			end
		end

		return xgra
	end
# -----------------------------------------------------------------------------------
# Nesterov Momentum
# -----------------------------------------------------------------------------------
	function gradient_nesterov(f, ∇f, x0; α = 0.0003, β = 0.9, ε_G = 1e-8, N = 500)
		xgra = zeros(length(x0), N)
		xgra[:, 1] = x0
		v = zeros(2,1); β = β; α = α
		for i in 2:N
			v = β*v - α*∇f(xgra[:,i-1] + β*v)
			xgra[:, i] = xgra[:, i-1] + v
			if norm(∇f(xgra[:, i])) ≤ ε_G
				return xgra[:, 1:i]
			end
		end

		return xgra
	end
# -----------------------------------------------------------------------------------
# Adagrad 
# -----------------------------------------------------------------------------------
	function gradient_adagrad(f, ∇f, x0; α = 0.3, ϵ = 1e-10, ε_G = 1e-8, N = 500)
		xgra = zeros(length(x0), N)
		xgra[:, 1] = x0
		v = zeros(length(x0)); ϵ = ϵ; α = α
		s = zeros(length(x0))
		for i in 2:N
			g = ∇f(xgra[:,i-1])
			s += g.*g
			xgra[:, i] = xgra[:, i-1] - α*g ./ (sqrt.(s) .+ ϵ)
			if norm(∇f(xgra[:, i])) ≤ ε_G
				return xgra[:, 1:i]
			end
		end

		return xgra
	end
# -----------------------------------------------------------------------------------
# RMSProp
# -----------------------------------------------------------------------------------
	function gradient_RMSProp(f, ∇f, x0; α = 0.3, γ = 0.1,  ϵ = 1e-10, ε_G = 1e-8, N = 500)
		xgra = zeros(length(x0), N)
		xgra[:, 1] = x0
		v = zeros(length(x0)); ϵ = ϵ; α = α
		s = zeros(length(x0))
		for i in 2:N
			g = ∇f(xgra[:,i-1])
			s += γ*s + (1-γ)*(g.*g)
			xgra[:, i] = xgra[:, i-1] - α*g ./ (sqrt.(s) .+ ϵ)
			if norm(∇f(xgra[:, i])) ≤ ε_G
				return xgra[:, 1:i]
			end
		end

		return xgra
	end
end

# ╔═╡ 7df517d3-d8f0-4273-85e4-b81d2dd5804b
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/momentum_ex3.pdf", fig2)
	input = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/momentum_ex3.pdf"
	output = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/momentum_ex3.pdf"
	run(`pdfcrop $input $output`)
end
  ╠═╡ =#

# ╔═╡ 0d00c12b-88a8-4c5b-b031-0c164a6f886b
function plot_line(fig, f, row, col, xgra_list, colors, labels; 
				   legend_offset=(0,0), legend_pos=:rt, N = 100)
	ax = Axis(fig[row, col]; xlabel=L"$$Iteration", ylabel=L"f(x_1, x_2)",
			 aspect=AxisAspect(asp), backgroundcolor=(:cyan, 0.05))
	xlims!(ax, 0, N)
	for (xgra, color, lbl) in zip(xgra_list, colors, labels)
		Np = size(xgra, 2)
		z = [f(xgra[:,k]) for k in 1:Np]
		scatterlines!(ax, 1:Np, z, markersize=10, color=(color, 0.5), 
					  strokewidth=1, strokecolor=(color), 
					  label=lbl * " , Stop at $(round(z[end], digits=2))")
	end
	axislegend(ax; position=legend_pos, offset=legend_offset, labelsize=14)
	
	fig
end

# ╔═╡ 475f7464-9d2f-45da-82f9-a480c59ba028
function plot_con(fig, f, row, col, xgra_list, colors, labels; 
                  legend_offset=(0,0), legend_pos=:rt)

    ax = Axis(fig[row, col]; xlabel=L"x_1", ylabel=L"x_2",
              aspect=AxisAspect(asp), backgroundcolor=(:cyan, 0.05))
    limits!(ax, xlim[1], xlim[2], ylim[1], ylim[2])

    # Grid
    z = [f([x1, x2]) for x1 in xxs1, x2 in xxs2]
    contour!(ax, xxs1, xxs2, z; levels=level, linewidth=1)

    # Multiple trajectories
    for (xgra, color, lbl) in zip(xgra_list, colors, labels)
        Np = size(xgra, 2)
        scatterlines!(ax, xgra[1, 1:Np], xgra[2, 1:Np]; color=(color, 0.3), 
					  linewidth=2, label=lbl, strokewidth=1, strokecolor=color)
    end

    # One legend for all
    axislegend(ax; position=legend_pos, offset=legend_offset, labelsize=14)

    return fig
end


# ╔═╡ aa5384b0-98e3-44bd-b8e3-09d37916af59
begin	
	xgra1 = gradient_descent_fixed_alpha(f, ∇f, x0; α=α); nothing
	xgra2 = gradient_heavy_ball(f, ∇f, x0, α=α)
	if 	f == Rosenbrock
		xgra3 = gradient_nesterov(f, ∇f, x0; α=0.0002, β=0.92)
		legend_pos = :rt
	elseif f == Brunton
		xgra3 = gradient_nesterov(f, ∇f, x0; α=0.1, β=0.9)
		legend_pos = :lt
	end
	
	empty!(fig1)
	plot_con(fig1, f, 1, 1, [xgra1, xgra2, xgra3], [:red, :blue, :magenta],
         ["Gradient descent", "Heavy Ball", "Nestrov Momentum"]; 
			 legend_offset=(0, -30), legend_pos=legend_pos)
	plot_line(fig1, f, 1, 2, [xgra1, xgra2, xgra3], [:red, :blue, :magenta],
		["Gradient descent", "Heavy Ball", "Nestrov Momentum"]; 
			  legend_offset=(0, -30), legend_pos=:rt)
end

# ╔═╡ 692947a3-08bf-4bb2-8c9c-b393828c31d0
let
	global fig2 = Figure(size=(800,600))
	if 	f == Rosenbrock
		xgra4 = gradient_adagrad(f, ∇f, x0; α = 0.07)
		xgra5 = gradient_RMSProp(f, ∇f, x0; α = 0.1, γ = 0.0003)
		legend_pos = :rt
		N = 200
	elseif f == Brunton
		xgra4 = gradient_adagrad(f, ∇f, x0; α = 3)
		xgra5 = gradient_RMSProp(f, ∇f, x0; α = 3, γ = 0.1)
		legend_pos = :lt
		N = 100
	end
	
	empty!(fig2)

	plot_con(fig2, f, 1, 1, [xgra1, xgra4, xgra5], [:red, :blue, :magenta],
         ["Gradient descent", "Adagrad", "RMSProp"]; 
			 legend_offset=(0, -30), legend_pos=legend_pos)
	plot_line(fig2, f, 1, 2, [xgra1, xgra4, xgra5], [:red, :blue, :magenta],
		["Gradient descent", "Adagrad", "RMSProp"]; legend_offset=(0, -30), 
			  legend_pos=:rt, N=N)
end

# ╔═╡ Cell order:
# ╠═70563c10-824b-11f0-2115-e93c87dcd75c
# ╠═1dccf6e1-4108-4eac-afd3-8b7ba9288c05
# ╠═2dc11439-74c6-4fde-9f3d-18dfbb8f8712
# ╠═527deac5-07bd-487c-a8a0-200cbaa031b6
# ╟─9a78d420-1aaa-452e-b655-ca247e940a76
# ╟─98ac464d-f8a0-4b0f-9175-2cd185703b96
# ╟─ae7fda1b-4389-4524-b9b6-6ebf1d4ae754
# ╟─cacef432-87f7-4977-a602-3f2adfe5b9a2
# ╠═3decb336-9894-40ed-86d3-ef23038827c4
# ╟─22b7adf3-7296-44e7-beeb-9a037087ba7b
# ╠═86af2161-32f3-451e-a65f-ab6e5e4d62ec
# ╟─cd2ef601-9890-44a9-a01d-8ba09a5ad9c7
# ╠═e877de87-3b79-4213-ab2b-73e44c7d6f68
# ╠═43009dda-1630-4970-ba5e-e6f85dfd96fd
# ╠═aa5384b0-98e3-44bd-b8e3-09d37916af59
# ╠═98050695-1fd1-4ecf-994d-06c2c8f70279
# ╠═692947a3-08bf-4bb2-8c9c-b393828c31d0
# ╠═03dcf143-8168-40f0-a9ac-fefb23e742e8
# ╠═7df517d3-d8f0-4273-85e4-b81d2dd5804b
# ╠═0d00c12b-88a8-4c5b-b031-0c164a6f886b
# ╠═475f7464-9d2f-45da-82f9-a480c59ba028
