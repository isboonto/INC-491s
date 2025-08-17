### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ a524ef58-7a77-11f0-2e68-835aa5683d02
using Pkg; Pkg.activate()

# ╔═╡ 1938e7ca-38c7-4948-a0f3-3fd12e9d70cb
begin
	using LinearAlgebra 
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
end

# ╔═╡ d4ea2e1c-c265-4307-888a-f8834775c31f
begin
	
	# Define the function
	f(x, y) = x^2 + y^2
	
	# Grid for contour plot
	x = -2:0.05:2
	y = -2:0.05:2
	Z = [f(xi, yi) for yi in y,  xi in x]   # NOTE: order is (x,y) in Makie
	
	# Point of interest
	x0, y0 = 1.0, 1.0
	grad = [2x0, 2y0]
	grad_norm = grad / norm(grad)
	
	# Directions
	d_grad = grad_norm
	d_neg_grad = -grad_norm
	d_descent = [-grad_norm[2], grad_norm[1]] - grad_norm
	
	# Scale factor to make arrows smaller
	scale = 1
	
	# --- Plot ---
	fig = Figure(size=(500,500))
	ax = Axis(fig[1,1], aspect=1, xlabel=L"$x_1$", ylabel=L"$x_2$")
	
	contour!(ax, x, y, Z; levels=20, colormap=:grays)
	scatter!(ax, [x0], [y0], color=:red, markersize=15, label="Point x")
	
	# Arrows (Makie uses arrows!(x, y, dx, dy))
	arrows!(ax, [x0], [y0],
	        [scale*d_grad[1]], [scale*d_grad[2]],
	        arrowsize=15, linewidth=2, color=:red, label="∇f(x)")
	
	arrows!(ax, [x0], [y0],
	        [scale*d_neg_grad[1]], [scale*d_neg_grad[2]],
	        arrowsize=15, linewidth=2, color=:blue, label="-∇f(x)")
	
	arrows!(ax, [x0], [y0],
	        [scale*d_descent[1]], [scale*d_descent[2]],
	        arrowsize=15, linewidth=2, color=:green, label="Other descent")
	
	axislegend(ax, position=:rb)
	fig
	
end

# ╔═╡ bac41322-5d77-488d-85c1-6411d6d0efe0
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/figures/descent_direction.png", fig)

# ╔═╡ Cell order:
# ╠═a524ef58-7a77-11f0-2e68-835aa5683d02
# ╠═1938e7ca-38c7-4948-a0f3-3fd12e9d70cb
# ╠═d4ea2e1c-c265-4307-888a-f8834775c31f
# ╠═bac41322-5d77-488d-85c1-6411d6d0efe0
