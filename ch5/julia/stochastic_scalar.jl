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

# ╔═╡ b40d3cb4-86d9-11f0-3903-7946b5f673e5
using Pkg; Pkg.activate()

# ╔═╡ 741ac103-8199-4ec6-a790-47bd6c56aca3
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using ForwardDiff, LinearAlgebra
	using PlutoUI, LaTeXStrings, Symbolics, Optim, Latexify
end

# ╔═╡ 0f8fc61b-2fab-42fa-a1ec-ae4818175f3a
begin
	Rosenbrock(x) = (1 - x[1])^2 + 100*(4*x[2] - x[1]^2)^2
	
	# Boyd and Vanderberghe
	Brunton(x) =  15*exp(-0.25*(x[1] - 2)^2) + 0.5*exp(-2(x[1] - 2)^2) + x[1]^2 + x[2]^2

end

# ╔═╡ b354454e-8923-4399-8cb1-e93ed8efc48a
md"""
**1. Select a function:** $(@bind f Select([Rosenbrock, Brunton]))
"""

# ╔═╡ d2ddff30-5e6d-4291-8c68-49c917f20bd8
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

# ╔═╡ b2f496a8-7eb1-49c6-b3ed-f00f2c333583
md"""
**2. Select a calculation:**
"""

# ╔═╡ 07cbb4df-9b03-48bb-a785-f7820f5b0028
@variables x₁, x₂

# ╔═╡ 30bc0071-3c3e-45f3-aaa5-d5f9278f5cc8
let
	original_str = latexify(f([x₁, x₂]))
	modified_str = "f(x_1, x_2) =" * original_str
	latexstring(modified_str)
end

# ╔═╡ 7d118708-1aef-4c96-b7e2-3b0320c28948
∇f = x -> ForwardDiff.gradient(f, x)

# ╔═╡ 8271eb6a-51b3-4c66-9bb3-75430822e9ee
let
	original_str = "∇f(x_1, x_2) = " * latexify(∇f([x₁, x₂]))
	latexstring(original_str)
end

# ╔═╡ 007d606a-d399-4bbd-8ae0-51ff3896411e


# ╔═╡ Cell order:
# ╠═b40d3cb4-86d9-11f0-3903-7946b5f673e5
# ╠═741ac103-8199-4ec6-a790-47bd6c56aca3
# ╠═0f8fc61b-2fab-42fa-a1ec-ae4818175f3a
# ╠═b354454e-8923-4399-8cb1-e93ed8efc48a
# ╟─d2ddff30-5e6d-4291-8c68-49c917f20bd8
# ╟─b2f496a8-7eb1-49c6-b3ed-f00f2c333583
# ╟─07cbb4df-9b03-48bb-a785-f7820f5b0028
# ╟─30bc0071-3c3e-45f3-aaa5-d5f9278f5cc8
# ╠═7d118708-1aef-4c96-b7e2-3b0320c28948
# ╠═8271eb6a-51b3-4c66-9bb3-75430822e9ee
# ╠═007d606a-d399-4bbd-8ae0-51ff3896411e
