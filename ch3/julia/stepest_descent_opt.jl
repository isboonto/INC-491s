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

# ╔═╡ aa9d275a-7d03-11f0-0253-6bb58f31d4cf
using Pkg; Pkg.activate()

# ╔═╡ c294ad00-0158-415b-89f5-78748097465f
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using ForwardDiff, LinearAlgebra
	using PlutoUI
end

# ╔═╡ 8b46ab10-a2d4-48c5-9441-57519275b1c1
#@bind f Select([Quadratic, Bean, Rosenbrock, Boyd])

# ╔═╡ 3d198172-5ebb-4ca5-b9d0-6801355b22d5
md"""
γ = $(@bind γ PlutoUI.Slider([1,2,3,4, 5, 6, 7, 8], show_value=true, default= 1))
"""

# ╔═╡ c6b6a134-8dcd-44b1-a00e-0a1e53816d58
begin
	# quadratic function
	Quadratic = x->(1/2) * x' * [2 0; 0 2γ] * x;

	# Bean function
	Bean = x->(1-x[1])^2 + (1-x[2])^2 + (0.5)*(2x[2] - x[1]^2)^2;

	# Rosenbrock function
	Rosenbrock = x -> (1 - x[1])^2 + 0.5*(x[2] -x[1]^2)^2;

	# Boyd and Vanderberghe
	Boyd = x -> exp(x[1] + 2x[2] -2) + exp(x[1] - 2x[2] - 2) + exp(-x[1] - 2)
end

# ╔═╡ ec1b44c6-d19d-4c04-980b-c51eaeb2ee17
begin
	f = x -> (1/2) * x' * [2 0; 0 2γ] * x
	∇f = x -> ForwardDiff.gradient(f, x)
	H = x -> ForwardDiff.hessian(f, x)
end

# ╔═╡ b277513d-80fd-4590-939a-eb6cf9ab2f73
function gradient_descent_opt(f, ∇f, x0; α = 0.1, ε_G=1e-6, N=100)
	xgra = zeros(length(x0), N)
	xgra[:, 1] = x0

	for i in 2:N
		xgra[:, i] = xgra[:, i-1] - α * ∇f(xgra[:, i-1])
		if norm(∇f(xgra[:, i])) ≤ ε_G
			return xgra[:, 1:i]
		end
	end

	return xgra
end

# ╔═╡ fc727865-959d-46e7-b556-941d545f8bf2
function plot_gd(α, pf, row, col, xrgra, color::Symbol, fig, x0)
	# Mesh Grid
	xxs1 = -5:0.05:4;
	xxs2 = -5:0.05:4; 
	
	ax1 = CairoMakie.Axis(fig[row, col], xlabel=L"x_1", ylabel=L"x_2",
						 aspect = AxisAspect(1.3), backgroundcolor=(:cyan, 0.05))
	limits!(ax1, -2.2, 2.2, -2.2, 3.2)

	# Contour and initial point
	z = [ f([x, y]) for x in xxs1, y in xxs2]
	lv1 = -0:5:200
	contour!(ax1, xxs1, xxs2, z, levels = lv1, linewidth=1)
	contour!(ax1, xxs1, xxs2, z, levels = [2, 0.6, 0.11], linewidth=1)

	text!(ax1, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	scatter!(ax1, x0[1], x0[2], markersize=15, strokewidth=1, strokecolor=:blue,
			color=:lightblue)
	# Scatter line and the optimal point.
	Np1 = size(xrgra, 2)
	scatterlines!(ax1, xrgra[1,1:Np1], xrgra[2, 1:Np1], color=:blue,
				  markercolor=(color, 0.5), markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=2, label=("α = $(round(α, digits=2)) with $(Np1-1) iterations"))
	scatter!(ax1, xrgra[1,end], xrgra[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(ax1, xrgra[1,end]-0, xrgra[2,end]+0.15, text=L"x^\ast", fontsize=22)


	axislegend(ax1, labelsize=14)
	fig
end

# ╔═╡ a7dd0e23-3c0e-494c-bd36-340a3edc76f2
function eigen_bounds(x::AbstractVector{<:Real}, f)
	H = ForwardDiff.hessian(f, x)
	vals = eigen(H).values

	return maximum(vals), minimum(vals)
end

# ╔═╡ d51157d8-d36b-4322-bcad-021b55ad2751
begin
		fig1 = Figure(size=(800,400))
	
		x0 = [-1.2, 2.0]
		λ_max, λ_min = eigen_bounds(x0, f)
		α1 = 2/(λ_max + λ_min)
		xrgra1 = gradient_descent_opt(f, ∇f, x0; α=α1)
		plot_gd(α1, f, 1, 1, xrgra1, :red, fig1, x0)
	
end

# ╔═╡ e99f5923-0815-4972-9b62-92ec8cc9b53f
# ╠═╡ disabled = true
#=╠═╡
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/optimal_step_size2.pdf", fig1)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═aa9d275a-7d03-11f0-0253-6bb58f31d4cf
# ╠═c294ad00-0158-415b-89f5-78748097465f
# ╠═c6b6a134-8dcd-44b1-a00e-0a1e53816d58
# ╠═8b46ab10-a2d4-48c5-9441-57519275b1c1
# ╠═ec1b44c6-d19d-4c04-980b-c51eaeb2ee17
# ╠═3d198172-5ebb-4ca5-b9d0-6801355b22d5
# ╠═d51157d8-d36b-4322-bcad-021b55ad2751
# ╠═e99f5923-0815-4972-9b62-92ec8cc9b53f
# ╠═b277513d-80fd-4590-939a-eb6cf9ab2f73
# ╠═fc727865-959d-46e7-b556-941d545f8bf2
# ╠═a7dd0e23-3c0e-494c-bd36-340a3edc76f2
