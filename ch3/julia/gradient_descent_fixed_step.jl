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

# ╔═╡ c10e4bb5-6601-4f84-9481-ae30f461b1dd
md"""
### Objective function
"""

# ╔═╡ d6ce048a-29a0-4fe9-8d54-07b1c6115c10
md"""
Quadratic function:

vertical axis β₁ = $(@bind β₁ PlutoUI.Slider(1:1:20; show_value=true, default=2)), horizontal β₂  $(@bind β₂ PlutoUI.Slider(1:1:20; show_value=true, default=3))
"""

# ╔═╡ a711c0f4-823b-45f5-a074-7a7e4218b7a4
begin
	# quadratic function
	Quadratic(x) = (1/2)*[x[1]; x[2]]'*[β₁ 1; 1 β₂]*[x[1]; x[2]]
	

	# Bean function
	Bean(x; a=1, b=0.5) = (a-x[1])^2 + (1-x[2])^2 + (b)*(2x[2] - x[1]^2)^2; nothing

	# Rosenbrock function
	Rosenbrock(x; a=1, b=0.5) = (a - x[1])^2 + b*(x[2] -x[1]^2)^2;

	# Boyd and Vanderberghe
	Boyd(x) = exp(x[1] + 2x[2] -2) + exp(x[1] - 2x[2] - 2) + exp(-x[1] - 2)
end

# ╔═╡ 54ffce70-0ac7-4783-a61f-acf577750def
@bind f Select([Quadratic, Bean, Rosenbrock, Boyd])

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
	println("define figure")
end

# ╔═╡ d38a6184-9470-4927-a256-cafb2ae562fa
begin
	md"""
	x = $(@bind x01 PlutoUI.Slider(-2:0.1:2, show_value=true, default = -1.2)), $(@bind x02 PlutoUI.Slider(-2:0.1:2, show_value=true, default=2))
	
	α3 = $(@bind α3 PlutoUI.Slider([0.01, 0.05, 0.1, 0.25, 0.3, 0.4, 0.543, 0.7], show_value=true, default= 0.543))
	"""
end

# ╔═╡ 24d72df5-3cc0-4475-b36b-3aad1c693ee0
# ╠═╡ disabled = true
#=╠═╡
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/fixed_step_size.pdf", fig1)
  ╠═╡ =#

# ╔═╡ adf60444-fd6a-42d2-afe0-737220c76046
md"""
### Zig-Zags Property
"""

# ╔═╡ a366c24e-cdbe-483e-b0a7-08d91616b991
md"""
α4 = $(@bind α4 PlutoUI.Slider([0.01, 0.05, 0.1, 0.25, 0.3, 0.4, 0.543, 0.6, 0.7, 1.0], show_value=true, default= 0.543))
"""

# ╔═╡ 75863c55-09cf-47ea-b4dd-5e9d7a757c18
# ╠═╡ disabled = true
#=╠═╡
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/zig_zags.pdf", fig2)
  ╠═╡ =#

# ╔═╡ 8f0b493c-29f4-4ee3-a1fc-a13a1e4859ae
function eigen_bounds(x::AbstractVector{<:Real}, f)
	H = ForwardDiff.hessian(f, x)
	vals = eigen(H).values

	return maximum(vals), minimum(vals)
end

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
function plot_gd(α, row, col, xrgra, color::Symbol, fig)
	ax1 = CairoMakie.Axis(fig[row, col], xlabel=L"x_1", ylabel=L"x_2",
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
	fig
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

	function step_linesearch!(M::GradientDescent, f, ∇f, x)
		α, g = M.α, ∇f(x)
		λ_max, λ_min = eigen_bounds(x, f)
		α = 2/(λ_max + λ_min)

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

# ╔═╡ 8d21f692-2967-4274-897a-51ac4982c658
function gradient_descent_linesearch(xrga0, α)
	
	Mg = GradientDescent(α) 	# initial α is 0
	Mg = init!(Mg, f, ∇f, x0)
	
	# Next Step
	for i = 2:N
		xgra[:,i] = step_linesearch!(Mg, f, ∇f, xgra[:,i-1])
		if norm(∇f(xgra[:,i])) <= ε_G
			global xrgra0 = xgra[:,1:i]
			break;
		else
			global xrgra0 = xgra[:, 1:i-1]
		end
	end
	return xrgra0
end

# ╔═╡ b818ee57-4584-49bf-b05c-a30bd21a7610
begin
	α1 = 0.01; α2 = 0.4;# α3 = 0.25;
	xrgra1 = zeros(length(x0), N);
	xrgra1 = gradient_descent(xrgra1, α1)

	xrgra2 = zeros(length(x0), N);
	xrgra2 = gradient_descent(xrgra2, α2)
	
	xrgra3 = zeros(length(x0), N);
	xrgra3 = gradient_descent(xrgra3, α3)

	xrgra4 = zeros(length(x0), N);
	xrgra4 = gradient_descent(xrgra4, α4)

	xrgra5 = zeros(length(x0), N);
	γ = 2
	fx = x -> x'*[2 0; 0 2γ]*x 
	λ_max, λ_min = eigen_bounds(x0, fx)
	α5 = 2/(λ_max + λ_min)
	xrgra5 = gradient_descent_linesearch(xrgra5, α5)
end

# ╔═╡ a9c666f6-704c-4be8-baa1-d4940fc004ac
begin
	fig2 = Figure(size=(600,600))
	plot_gd(α4, 1, 1, xrgra4, :red, fig2)
end

# ╔═╡ 5a80f411-f53a-4742-9a45-1f550f207a55
begin
	fig3 = Figure(size=(600,600))
	plot_gd(α5, 1, 1, xrgra5, :red, fig3)
end

# ╔═╡ 04be0f27-9cfd-4b98-a0a7-71c941ede3c2
function plot_trend(row, col, fig)
	ax1 = CairoMakie.Axis(fig[row, col], xlabel=L"$$Iteration", ylabel=L"$$O
						  bjective Function Value", aspect = AxisAspect(1.3), backgroundcolor=(:cyan, 0.05))
	
	Np1 = size(xrgra1,2)
	z1 = [f([xrgra1[1,x], xrgra1[2,x]]) for x in 1:Np1]
	scatterlines!(ax1, 1:Np1, z1, color=:blue,
				  markercolor=:magenta, markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=1, label=("α = $(α1)"))
	Np1 = size(xrgra2,2)
	z2 = [f([xrgra2[1,x], xrgra2[2,x]]) for x in 1:Np1]
	scatterlines!(ax1, 1:Np1, z2, color=:blue,
				  markercolor=:lightblue, markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=1, label=("α = $(α2)"))
	
	Np1 = size(xrgra3,2)
	z3 = [f([xrgra3[1,x], xrgra3[2,x]]) for x in 1:Np1]
	scatterlines!(ax1, 1:Np1, z3, color=:blue,
				  markercolor=:red, markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=1, label=("α = $(α3)"))
	
	
	max_value = maximum(filter(!isnan, vcat(z1, z2, z3)))
	if isinf(max_value)
		max_value = 15
	end
	limits!(ax1, 0, 30, -max_value/10, max_value + max_value/10)
	
	axislegend(ax1, labelsize=14)
	
	fig
end

# ╔═╡ 18783fa7-e31a-4233-8846-064e8b8eb601
begin
	empty!(fig1)
	plot_gd(α1, 1, 1, xrgra1, :magenta, fig1)
	plot_gd(α2, 1, 2, xrgra2, :lightblue, fig1)
	plot_gd(α3, 2, 1, xrgra3, :red, fig1)
	plot_trend(2,2, fig1)
end

# ╔═╡ Cell order:
# ╠═6d12590c-7b74-11f0-1629-a916ddf947dc
# ╠═ba33c9b6-8724-4b9e-9880-4f089772f481
# ╠═5e7b2f91-17a9-4793-aff3-b9ed8d8e68e0
# ╟─c10e4bb5-6601-4f84-9481-ae30f461b1dd
# ╠═a711c0f4-823b-45f5-a074-7a7e4218b7a4
# ╠═7514c095-cab9-484b-8c1d-171f067fffd0
# ╟─d6ce048a-29a0-4fe9-8d54-07b1c6115c10
# ╠═43682112-7a32-48bc-9d57-d7b942bed5c4
# ╠═54ffce70-0ac7-4783-a61f-acf577750def
# ╟─2b75caea-40c8-46fd-b484-35dd12e37e10
# ╟─d38a6184-9470-4927-a256-cafb2ae562fa
# ╠═18783fa7-e31a-4233-8846-064e8b8eb601
# ╠═24d72df5-3cc0-4475-b36b-3aad1c693ee0
# ╟─adf60444-fd6a-42d2-afe0-737220c76046
# ╟─a366c24e-cdbe-483e-b0a7-08d91616b991
# ╠═a9c666f6-704c-4be8-baa1-d4940fc004ac
# ╠═75863c55-09cf-47ea-b4dd-5e9d7a757c18
# ╠═5a80f411-f53a-4742-9a45-1f550f207a55
# ╠═8f0b493c-29f4-4ee3-a1fc-a13a1e4859ae
# ╟─55828d5d-2cc2-4bfd-b179-6a8536dbdf0a
# ╟─04be0f27-9cfd-4b98-a0a7-71c941ede3c2
# ╠═7f2952a5-1ff4-43ec-88b8-bcfe84cb4936
# ╠═61097e5a-095d-4e2c-8c03-50a15a6e2bbb
# ╠═94513c71-e7a3-49d3-a8a0-48d14848d8e7
# ╠═b45feff8-105f-4fe0-b7bd-9403c00752bd
# ╠═8d21f692-2967-4274-897a-51ac4982c658
# ╠═b818ee57-4584-49bf-b05c-a30bd21a7610
# ╠═b67ab6c4-0ad7-4884-84a2-269af103de4b
