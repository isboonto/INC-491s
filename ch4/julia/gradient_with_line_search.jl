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

# ╔═╡ ed326088-7e86-11f0-2d71-6f573ba75fb0
using Pkg; Pkg.activate()

# ╔═╡ 8b5e3fa3-027e-4571-9f13-e53820b5800b
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using Optim, LinearAlgebra, ForwardDiff
	using Symbolics, PlutoUI, LaTeXStrings
end

# ╔═╡ 18e81805-92de-4ea5-9e89-2d6873347ff9
md"""
## Steepest Descent with Exact and Approximate Line Search
"""

# ╔═╡ e37bd02c-f474-4329-92c7-215fb3f5b2df
md"""
vertical axis β₁ = $(@bind β₁ PlutoUI.Slider(1:1:8, show_value=true, default=2)) 
horizontal axis β₂ = $(@bind β₂ PlutoUI.Slider(1:1:8, show_value=true, default=2))
"""

# ╔═╡ db418715-334b-4171-81f1-24c4e9b0533a
begin
	# quadratic function
	Quadratic(x) = (1/2) * x' * [β₁ 1; 1 β₂] * x;

	# Bean function
	Bean(x) = (1-x[1])^2 + (1-x[2])^2 + (0.5)*(2x[2] - x[1]^2)^2;

	# Rosenbrock function
	Rosenbrock(x) = (1 - x[1])^2 + 0.5*(x[2] -x[1]^2)^2;

	# Boyd and Vanderberghe
	Boyd(x) =  exp(x[1] + 2x[2] -2) + exp(x[1] - 2x[2] - 2) + exp(-x[1] - 2)
end

# ╔═╡ 4987f44d-a95b-4b31-b958-a8bf88a2b0f5
md"""
**1. Select a function:** $(@bind f Select([Quadratic, Bean, Rosenbrock,  Boyd]))
"""

# ╔═╡ 61013f4b-77a0-464e-9178-54f2845f7bda
begin
	@variables x₁ x₂
	
	original_str = "f(x_1, x_2) = $(f([x₁, x₂]))"
	modified_str = replace(original_str, r"exp\((.*?)\)" => s"e^{\1}")
	modified_str = replace(original_str, r"exp\((.*?)\)" => s"e^{\1}") |>
               str -> replace(str, "*" => " ")
	
	ft = latexstring(modified_str)
end

# ╔═╡ 73922087-d2d1-43e8-b18b-ac990340aca4
fig1 = Figure(size=(600, 400));

# ╔═╡ dc94ac48-f82c-47b1-8593-8644bbb1f33e
begin
	∇f = x -> ForwardDiff.gradient(f, x)
	H = x -> ForwardDiff.hessian(f, x)
end

# ╔═╡ ec12c674-0c88-4fa2-9411-68bacd933eb7
function plot_con(row, col, fig, pf, x0, xrgra, color::Symbol, α)
	# Mesh Grid
	xxs1 = -5:0.05:4;
	xxs2 = -5:0.05:4;

	ax1 = CairoMakie.Axis(fig[row, col], xlabel=L"x_1", ylabel=L"x_2",
						 aspect = AxisAspect(1.2), backgroundcolor=(:cyan, 0.05))
	limits!(ax1, -2.2, 2.2, -1, 3.2)

	# Contourand initial point
	z = [f([x, y]) for x in xxs1, y in xxs2]
	lv1 = -0:5:200
	contour!(ax1, xxs1, xxs2, z, levels = lv1, linewidth=1)
	contour!(ax1, xxs1, xxs2, z, levels = [ 2, 0.6, 0.11, 0.01], linewidth=1)

	text!(ax1, x0[1], x0[2]+0.05, text=L"x_0", fontsize=22)
	# Scatter line and the optimal point.
	Np1 = size(xrgra, 2)
	scatterlines!(ax1, xrgra[1,1:Np1], xrgra[2, 1:Np1], color=:blue,
				  markercolor=(color, 0.5), markersize=15, strokewidth=1, 	
				  strokecolor=:blue, linewidth=2, label=("α = $(round(α, digits=2))  with $(Np1-1) iterations"))
	scatter!(ax1, xrgra[1,end], xrgra[2,end], markersize=10, color=:red,
		strokewidth=1, strokecolor=:black)
	text!(ax1, xrgra[1,end]-0, xrgra[2,end]+0.15, text=L"x^\ast", fontsize=22)


	axislegend(ax1, labelsize=14)
	fig
end

# ╔═╡ faf0c40d-dff5-4b6c-9833-24f153df38c0
function eigen_bounds(x::AbstractVector{<:Real}, f)
	H = ForwardDiff.hessian(f, x)
	
	if any(isnan, H)
		return 0, 0
	else
		vals = eigen(H).values
		return maximum(vals), minimum(vals)
	end
end

# ╔═╡ 3ba4b4d4-a7b7-46dd-b7bc-0cd0fa24c993
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

# ╔═╡ 0fabdeec-4bfa-4b9c-98c0-23c274d9fbaf
# Exact line search
function line_search(f, x, d)
	objective = α -> f(x + α*d)
	a, b = bracket_minimum(objective)
	
	# using brent method from Optim.jl
	res = Optim.optimize(objective, a, b)   # minimize
	α = Optim.minimizer(res)
		
	return α
end

# ╔═╡ 91269bf6-4e40-4fcc-932b-49f1c05c71dd
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

# ╔═╡ 0adf5908-458e-4356-ae0b-dfc5f7d0e00f

function backtracking_line_search(f, ∇f, x, d, α; ρ=0.5, β=1e-4)
	y, g = f(x), ∇f(x) 
	while f(x + α*d) > y + β*α*(g'*d)
		α *= ρ
	end
	
	return α
end


# ╔═╡ a01618e3-472b-4d13-ac89-9549933948cb
function gradient_descent_backtracking(f, ∇f, x0; α = 0.1, ε_G=1e-6, N=500)
	xgra = zeros(length(x0), N)
	xgra[:, 1] = x0

	for i in 2:N
		d = -∇f(xgra[:, i-1])
		α =  backtracking_line_search(f, ∇f, xgra[:,i-1], d, 10; ρ=0.5, β=1e-4)
		
		xgra[:, i] = xgra[:, i-1] + α * d
		if norm(∇f(xgra[:, i])) ≤ ε_G
			return xgra[:, 1:i]
		end
		
	end
	
	return xgra
end

# ╔═╡ d9d735b3-7962-4fae-abe5-a24348e5974a
function gradient_descent_line_search(f, ∇f, x0; α = 0.1, ε_G=1e-6, N=500)
	xgra = zeros(length(x0), N)
	xgra[:, 1] = x0

	for i in 2:N
		d = -∇f(xgra[:, i-1])
		α = line_search(f, xgra[:, i-1], d)
		xgra[:, i] = xgra[:, i-1] + α * d
		if norm(∇f(xgra[:, i])) ≤ ε_G
			return xgra[:, 1:i]
		end
	end
	
	return xgra
end

# ╔═╡ fef18a1c-9850-45f2-9def-0bba541fc2af
function gradient_descent_bisection(f, ∇f, x0; α = 0.1, ε_G=1e-6, N=500)
	xgra = zeros(length(x0), N)
	xgra[:, 1] = x0
	
	
	for i in 2:N
		d = -∇f(xgra[:, i-1])
		
		objective = α -> f(xgra[:, i-1] + α*d)
		a, b = bracket_minimum(objective)
		α = bisection(objective, a, b, 1e-3)
		
		xgra[:, i] = xgra[:, i-1] + α * d
		if norm(∇f(xgra[:, i])) ≤ ε_G
			return xgra[:, 1:i]
		end
		
	end
	
	return xgra
end

# ╔═╡ 1e03e205-e80a-4739-9750-0b4bc4b87a5b
function_map = Dict(
	"1. Exact Line Search" => gradient_descent_line_search,
	"2. Backtracking Line Search" => gradient_descent_backtracking,
	"3. Bisection Line Search" => gradient_descent_bisection,
)

# ╔═╡ ea91b07b-465a-4d41-9c0c-52ca130eaf05
@bind selected_function_label Select(collect(keys(function_map)))

# ╔═╡ ae591bc9-efe0-4dfc-bf88-5abc7a606592
selected_function = function_map[selected_function_label]

# ╔═╡ d42279d7-250a-4a63-981c-7c8f299eaea7
md"""
**2. Select a calculation:** **$selected_function**
"""

# ╔═╡ a7b41e78-05e6-4e3a-9d91-44327443af7a
begin
	empty!(fig1)
	x0 = [-1.2, 2.0]
	λ_max, λ_min = eigen_bounds(x0, f)
	#α1 = 2/(λ_max + λ_min)
	α1 = 0
	xgra1 = selected_function(f, ∇f, x0; α=α1)
	plot_con(1,1, fig1, f, x0, xgra1, :red, α1)
end

# ╔═╡ Cell order:
# ╠═ed326088-7e86-11f0-2d71-6f573ba75fb0
# ╟─18e81805-92de-4ea5-9e89-2d6873347ff9
# ╠═8b5e3fa3-027e-4571-9f13-e53820b5800b
# ╠═db418715-334b-4171-81f1-24c4e9b0533a
# ╠═4987f44d-a95b-4b31-b958-a8bf88a2b0f5
# ╠═61013f4b-77a0-464e-9178-54f2845f7bda
# ╟─d42279d7-250a-4a63-981c-7c8f299eaea7
# ╠═ea91b07b-465a-4d41-9c0c-52ca130eaf05
# ╠═ae591bc9-efe0-4dfc-bf88-5abc7a606592
# ╟─e37bd02c-f474-4329-92c7-215fb3f5b2df
# ╟─73922087-d2d1-43e8-b18b-ac990340aca4
# ╠═a7b41e78-05e6-4e3a-9d91-44327443af7a
# ╟─1e03e205-e80a-4739-9750-0b4bc4b87a5b
# ╠═dc94ac48-f82c-47b1-8593-8644bbb1f33e
# ╠═ec12c674-0c88-4fa2-9411-68bacd933eb7
# ╠═faf0c40d-dff5-4b6c-9833-24f153df38c0
# ╠═3ba4b4d4-a7b7-46dd-b7bc-0cd0fa24c993
# ╠═0fabdeec-4bfa-4b9c-98c0-23c274d9fbaf
# ╠═91269bf6-4e40-4fcc-932b-49f1c05c71dd
# ╠═0adf5908-458e-4356-ae0b-dfc5f7d0e00f
# ╠═a01618e3-472b-4d13-ac89-9549933948cb
# ╠═d9d735b3-7962-4fae-abe5-a24348e5974a
# ╠═fef18a1c-9850-45f2-9def-0bba541fc2af
