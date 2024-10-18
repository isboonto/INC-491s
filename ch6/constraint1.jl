### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 2f4683f8-8c5e-11ef-3ff0-7107adf83b3c
using Pkg; Pkg.activate()

# ╔═╡ c18a845b-8c9d-41da-a528-fa9b6d60be90
using CairoMakie

# ╔═╡ 89119540-f1ae-4c6c-a184-4a509081fa65
begin
	f(x1, x2) = x1^2 - (1/2)*x1 - x2 - 2
	g1(x1, x2) = x1^2 - 4x1 + x2 + 1
	g2(x1, x2) = (1/2)x1^2 + x2^2 - x1 - 4
end

# ╔═╡ 31fdb1ba-e671-4ae7-bf95-705dcb9a5a7c
begin
	fig1 = Figure(size=(600,400), fontsize=18)
	set_theme!(theme_latexfonts())

	ax1 = CairoMakie.Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$")

	xs = -2.5:0.01:5; ys = -2.5:0.01:5
	lev = collect(-10:1:30)
	contour!(ax1, xs, ys, f, levels = lev, colormap=:jet)#(:blue, 0.3))
	#contour!(ax1, xs, ys, g1, levels = [0], linewidth = 2, color=:red, label=L"$g_1(\mathbf{x})$", labelcolor=:red)
	#contour!(ax1, xs, ys, g2, levels = [0], linewidth = 2, color=:green, label=L"$g_2(\mathbf{x})$", labelcolor=:green)

	xs1_bound = -2.5:0.01:5
	xs2_bound = 1:0.01:3.3
	xs3_bound = -0.21:0.01:3.89
	xs4_bound = -0.31:0.01:-0.21
	xs5_bound = 3.89:0.01:4.5
	xs6 = -2:0.01:4
	
	g1_low = x -> -(x^2 - 4x + 1)  
	g2_low = x -> sqrt(-(1/2)x^2 + x + 4)
	
	lines!(ax1, xs, g1_low.(xs), color=:red, linewidth=2, label=L"$g_1(\mathbf{x})$")
	lines!(ax1, xs6, g2_low.(xs6), color=:green, linewidth=2, 
		label=L"$g_2(\mathbf{x})$")
	lines!(ax1, xs6, -g2_low.(xs6), color=:green , linewidth=2)
	band!(ax1, xs1_bound, g1_low.(xs1_bound), 5, color=(:brown, 0.2) )
	band!(ax1, xs2_bound, g2_low.(xs2_bound), g1_low.(xs2_bound), color=(:brown, 0.2) )
	band!(ax1, xs3_bound, -2.5, -g2_low.(xs3_bound), color=(:brown, 0.2))
	band!(ax1, xs4_bound, -2.5, g1_low.(xs4_bound), color=(:brown, 0.2))
	band!(ax1, xs5_bound, -2.5, g1_low.(xs5_bound), color=(:brown, 0.2))
	limits!(ax1, -2.5, 5, -2.5, 5)

	scatter!(ax1, [1.05],[2.1], marker=:circle, markersize=12, color=:blue, strokecolor=(:black, 0.5), strokewidth=2)
	text!(ax1, [0.9], [2.3], text=L"$\mathbf{x^\ast}$", color=:black)

	axislegend(ax1,fontsize=18)
	fig1
end

# ╔═╡ 14ce9e57-ca9f-46e9-a979-bcc2cd850d94
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/constrain1.pdf", fig1)
end

# ╔═╡ b607a3a1-e22b-43c6-87de-7c92306654bd
begin
	f4b(x1, x2, x3) = x1^2 + x2^2 - 2x1
	f4a(x1, x2, x3) = -x1 + x3 - 1
end

# ╔═╡ 614f6e60-98b5-4699-94d0-105553e6a62d
begin
	fig2 = Figure(fontsize=18)
	ax2 = Axis3(fig2[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", zlabel=L"$x_3$")

	xp = 0:0.01:2π
	x1 = 1 .+ cos.(xp); x2 = sin.(xp); x3 = 2 .+ cos.(xp)
	lines!(ax2, x1, x2, x3, linewidth=2)
	
	fig2
end

# ╔═╡ 65e5b885-c495-4942-a897-c1bf5fafdd65
begin
	fig3 = Figure(size=(600,400), fontsize=18)
	ax3 = CairoMakie.Axis(fig3[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
		yticks=[-2, -1, 0, 1, 2])

	f5(x1, x2) = (x2 + 100)^2 + 0.01*x1^2
	c5(x1, x2) = x2 - cos(x1)
	cx5 = x -> cos(x) 

	xs5 = -10:0.1:10
	ys5 = -4:0.1:4

	contour!(ax3, xs5, ys5, f5)
	lines!(ax3, xs5, cx5)

	fig3
end

# ╔═╡ Cell order:
# ╠═2f4683f8-8c5e-11ef-3ff0-7107adf83b3c
# ╠═c18a845b-8c9d-41da-a528-fa9b6d60be90
# ╠═89119540-f1ae-4c6c-a184-4a509081fa65
# ╠═31fdb1ba-e671-4ae7-bf95-705dcb9a5a7c
# ╠═14ce9e57-ca9f-46e9-a979-bcc2cd850d94
# ╠═b607a3a1-e22b-43c6-87de-7c92306654bd
# ╠═614f6e60-98b5-4699-94d0-105553e6a62d
# ╠═65e5b885-c495-4942-a897-c1bf5fafdd65
