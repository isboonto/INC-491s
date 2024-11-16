### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 7c7d889c-a3df-11ef-1699-f99f594e0205
using Pkg; Pkg.activate(@__DIR__)

# ╔═╡ 90ddcb7c-d784-4ac5-82bb-02bb71196fc8
begin
	using CairoMakie;
	using MathTeXEngine

	
	# Path to the Latin Modern font files
	font_path = "/usr/share/fonts/opentype/latex-xft-fonts/latinmodern-math.otf"

	#set_theme!(Theme(texfont = font_path))
	set_theme!(theme_latexfonts(), fontsize=20)
end

# ╔═╡ 35be0489-ca00-4f5c-a604-7697f19740d4
begin
	f =  (x1, x2) -> -x1 - x2
	g1 = (x1, x2) -> x1^2 - x2
	g2 = (x1, x2) -> x1^2 + x2^2 - 1
	#g1 = (x1, x2) -> x1^2 + (x2-1)^2 - 1
end

# ╔═╡ 84612895-9a53-405d-b1c7-f965872bce0a
begin
	fig1 = Figure(size=(600,400), fontsize=20)
	ax1 = Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", aspect=1.4)
	limits!(ax1, -1.5, 1.5, -0.25, 1.5)

	xs1 = -4:0.02:4; 
	contour!(ax1, xs1, xs1, f, levels=-20:0.5:20)
	contour!(ax1, xs1, xs1, g1, levels=[0], linewidth=2, color=:red)
	contour!(ax1, xs1, xs1, g2, levels=[0], linewidth=2, color=:blue)

	xb = -0.8:0.1:0.8
	g1b = (x) -> x^2 
	g2b = (x) -> sqrt(1-x^2)
	band!(ax1, xb, g2b.(xb), g1b.(xb), color=(:orange, 0.4))
	scatter!(ax1, 1/sqrt(2), 1/sqrt(2), markersize=10, strokecolor=:black, 
		strokewidth=1)
	text!(ax1, 0.72, 0.8, text=L"$\mathbf{x}$")
	fig1
end

# ╔═╡ 9162ffab-1f53-4ace-8dfb-d86c861e0503
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/ex5_han.png", fig1)
end

# ╔═╡ Cell order:
# ╠═7c7d889c-a3df-11ef-1699-f99f594e0205
# ╠═90ddcb7c-d784-4ac5-82bb-02bb71196fc8
# ╠═35be0489-ca00-4f5c-a604-7697f19740d4
# ╠═84612895-9a53-405d-b1c7-f965872bce0a
# ╠═9162ffab-1f53-4ace-8dfb-d86c861e0503
