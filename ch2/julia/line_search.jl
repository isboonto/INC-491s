### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ 82e29556-7b3c-11f0-0f82-df6e4b39e425
using Pkg; Pkg.activate()

# ╔═╡ ad11c6d0-6db2-4009-a854-252c5b42e700
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize = 18)
	using Symbolics, ForwardDiff, LinearAlgebra, Roots
	using PlutoUI
	using SymbolicUtils
end

# ╔═╡ 8334b98a-d068-4a78-89ac-482c0fd033f7
begin
	f = x -> x[1]^2 + 5x[2]^2
	∇f = x -> ForwardDiff.gradient(f, x)
end

# ╔═╡ a05a9e9e-7950-46cd-a2fb-8147d2944803
begin
	x0 = [3, 1]
	f([3, 1])
end

# ╔═╡ 5f90fa96-a94e-4c8f-95be-f1cdde5669b7
begin
	d = [-3, -5]
	∇f(x0)'*d
end

# ╔═╡ 02a6036c-18f2-4c3f-a80e-e6d4ecb2e11e
begin
	# f(α) = f(x + αd)
	@variables α
	f1 = α -> f(x0 + α*d)

	f1(α)
	df1 = α -> ForwardDiff.derivative(f1, α)
	df1(α)
end

# ╔═╡ 7b68bc2b-48fe-4319-927c-7430c6b66be8
begin
	α0 = find_zeros(df1, [-1, 1])[1]                     # vector to scalar
	x1 = x0 + α0*d 
end

# ╔═╡ 676850c2-c316-434d-b59c-6b250f1cc4a8
let
	global fig1 = Figure(size=(400, 400))
	ax1 = Axis(fig1[1,1], xlabel=L"x_1", ylabel=L"x_2")

	xs = ys = -5:0.1:5
	z = [f([x,y]) for x in xs, y in ys]
	x = [x0[1], x0[2]]; y = [x0[2], x1[2]]
	dx = diff(x); dy = diff(y)
	
	contour!(ax1, xs, ys, z, levels = 0:5:200)
	scatterlines!(ax1, x, y, strokewidth=1, strokecolor=:black, color=:red)

	# Compute midpoints and directions
	for i in 1:length(x)-1
    	xmid = (x[i] + x[i+1]) / 2
    	ymid = (y[i] + y[i+1]) / 2
    	dx = x[i+1] - x[i]
    	dy = y[i+1] - y[i]

    	arrows!(ax1, [xmid], [ymid], [dx], [dy]; arrowsize=10, linewidth=2, 
				color=:red)
	end
	text!(ax1, x[1], y[1]+0.2, text=L"\mathbf{x_0}")
	text!(ax1, x[2]-0.4, y[2]+0.2, text=L"\mathbf{x_1}")
end

# ╔═╡ acbc3f67-5027-4b32-92ff-20adfdb3c8a0
fig1

# ╔═╡ 0df1e015-6616-456d-b27f-09d70acc3d5f
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/figures/line_search2.pdf", fig1)

# ╔═╡ Cell order:
# ╠═82e29556-7b3c-11f0-0f82-df6e4b39e425
# ╠═ad11c6d0-6db2-4009-a854-252c5b42e700
# ╠═8334b98a-d068-4a78-89ac-482c0fd033f7
# ╠═a05a9e9e-7950-46cd-a2fb-8147d2944803
# ╠═5f90fa96-a94e-4c8f-95be-f1cdde5669b7
# ╠═02a6036c-18f2-4c3f-a80e-e6d4ecb2e11e
# ╠═7b68bc2b-48fe-4319-927c-7430c6b66be8
# ╠═676850c2-c316-434d-b59c-6b250f1cc4a8
# ╠═acbc3f67-5027-4b32-92ff-20adfdb3c8a0
# ╠═0df1e015-6616-456d-b27f-09d70acc3d5f
