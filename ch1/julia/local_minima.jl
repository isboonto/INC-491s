### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ f436fc48-76a6-11f0-339a-3fca0e87b170
using Pkg; Pkg.activate()

# ╔═╡ 713f47ee-c9ec-4d57-9a6d-dcade7fab641
using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)

# ╔═╡ 291d00e2-d4e3-4dfb-a967-0ef82f1111df
f = x -> (-5/6)x^3 + (7/2)x^2 - (11/3)x + 3

# ╔═╡ a92f89fc-31ef-40f3-a456-3a40417a4857
let
	x = 0:0.01:3
	z = [f(x) for x in x]

	global fig1 = Figure(size=(600, 400))
	ax = CairoMakie.Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$")
	limits!(ax, 0, 3, 1.0, 3.5)
	
	lines!(ax, x, z, linewidth=2) 
	x1 = 0.6972; x2 = 2.1024
	scatter!(ax, [x1, x2], [f(x1), f(x2)], color=:red, markersize=15,
			strokecolor=:blue, strokewidth=1)
	text!(ax, x1-0.4, f(x1)-0.3, text=L"$$local minimum")
	text!(ax, x1-0, f(x1)+0.1, text=L"$\mathbf{x^\ast}$")
	text!(ax, x2-0.4, f(x2)+0.3, text=L"$$local minimum")
	text!(ax, x2-0, f(x2)+0.1, text=L"$\mathbf{x^\ast}$")
end

# ╔═╡ 7cc6c9d2-7386-4750-85cf-36ba0c128cc0
fig1

# ╔═╡ 36fb9cf7-b725-4ff3-bee3-d1f4560da0c6


# ╔═╡ Cell order:
# ╠═f436fc48-76a6-11f0-339a-3fca0e87b170
# ╠═713f47ee-c9ec-4d57-9a6d-dcade7fab641
# ╠═291d00e2-d4e3-4dfb-a967-0ef82f1111df
# ╠═a92f89fc-31ef-40f3-a456-3a40417a4857
# ╠═7cc6c9d2-7386-4750-85cf-36ba0c128cc0
# ╠═36fb9cf7-b725-4ff3-bee3-d1f4560da0c6
