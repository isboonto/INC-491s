### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 0710492c-9759-11ef-2c65-f3dd16b09457
using Pkg; Pkg.activate()

# ╔═╡ e332a6d0-fcaa-4d62-a6e6-0135461e1d22
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using PlutoUI, Printf
	using Optim, LinearAlgebra
	using ForwardDiff
end

# ╔═╡ e67d5422-b549-4a3a-a3b0-cd13d424243e
begin
	f = x -> x[1] + 2x[2]
	g = x -> (1/4)*x[1]^2 + x[2]^2 - 1
	L = (x,λ, μ) -> f(x) + λ*g(x) + (μ/2)*g(x)^2

	# for band
	gcb = x -> sqrt(-(1/4)x^2 + 1)
end

# ╔═╡ 5c7b7401-97f7-42b0-8fd4-ec610f51de63
begin

	cs = 2; level1 = -20:cs:20			# contour space
	xopt = [-1.41434,  -0.707164]
	x1s = -4:0.01:4; x2s = -4:0.01:4
	μ = 0.50; λ = 0;
	Lc = [L([x1, x2], λ, μ) for x1 in x1s, x2 in x2s]
	fc = [f([x1, x2]) for x1 in x1s, x2 in x2s]
	gc = [g([x1, x2]) for x1 in x1s, x2 in x2s]

	fig = Figure(size=(800,600), fontsize=20)
	
	for (ii, k) in enumerate([0, 2, 5, 9])
		global xs_P = zeros(k+1,2)
		global x0 = x_start = [1, 1.5]
		xs_P[1,:] = x0
		global λ1 = λ_start = 0
		ρ = 1.1
		global μ1 = μ_start = 0.5
		global i = 1
		
		obj = x -> f(x) + λ1*g(x) + (μ1/2)*(g(x))^2
		global x1 = optimize(obj, x0).minimizer

		while (norm(ForwardDiff.gradient(obj, x1)) > 1e-10) && (i < k+1) 
			#obj = x -> f(x) + λ1*g(x) + (μ1/2)*(g(x))^2
			x1 = optimize(obj, x0).minimizer
			xs_P[i+1,:] = x1
			x0 = x1
			λ1 += μ1*g(x1)
			μ1 *= ρ
			i += 1
		end

		
		if k == 0
			ax1 = CairoMakie.Axis(fig[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L"k = %$(i-1),\quad \mu = %$(round(μ1, digits=2)),\quad \lambda = %$(round(λ1, digits=2))")
		elseif k == 2 
			ax1 = CairoMakie.Axis(fig[1,2], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L"k = %$(i-1),\quad \mu = %$(round(μ1, digits=2)),\quad \lambda = %$(round(λ1, digits=2))")
		elseif k == 5
			ax1 = CairoMakie.Axis(fig[2,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L"k = %$(i-1),\quad \mu = %$(round(μ1, digits=2)),\quad \lambda = %$(round(λ1, digits=2))")
		elseif k == 9
			ax1 = CairoMakie.Axis(fig[2,2], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L"k = %$(i-1),\quad \mu = %$(round(μ1, digits=2)),\quad \lambda = %$(round(λ1, digits=2))")
		end
		

		Lc1 = [L([x1, x2], λ1, μ1) for x1 in x1s, x2 in x2s]
	
		limits!(ax1, -3, 3, -2.5, 2.5)
		hidedecorations!(ax1, ticklabels=false, ticks=false, grid=true, label=false)
	
		x1b = -2:0.1:2; x1bl = -4:0.1:-2; x1br = 2:0.1:4
		band!(ax1, x1b, gcb.(x1b), 2.5*ones(length(x1b)), color=(:cyan, 0.4))
		band!(ax1, x1b, -2.5*ones(length(x1b)), -gcb.(x1b), color=(:cyan, 0.4))
		band!(ax1, x1bl, -2.5*ones(length(x1bl)), 2.5*ones(length(x1bl)), 
			color=(:cyan, 0.4))
		band!(ax1, x1br, -2.5*ones(length(x1br)), 2.5*ones(length(x1br)), 
			color=(:cyan, 0.4))
	
		contour!(ax1, x1s, x2s, fc, levels=-100:1:100, linewidth=2, linestyle=:dash)
		contour!(ax1, x1s, x2s, gc, levels=[0], color=:red, linewidth=2)
		scatter!(ax1, xopt[1], xopt[2], markersize = 15, color = :blue, 	
			strokewidth=1, strokecolor=:black)
		text!(ax1, xopt[1]+0.1, xopt[2], text=L"\mathbf{x}^\ast", fontsize=24, 
			color=:black)

	
		xi = [1,1.5]
		text!(ax1, xi[1] + 0.1, xi[2], text=L"\mathbf{x}_0", fontsize=24, 
			color=:black)
		contour!(ax1, x1s, x2s, Lc1, levels=level1, colormap = :viridis, 
			linewidth=1)
		scatterlines!(ax1, xs_P[1:i,1], xs_P[1:i,2], color=:red, strokewidth=1, 
			strokecolor=:black, markersize=15)
		text!(ax1, xs_P[i,1]+0.1, xs_P[i,2]-0.6, text=L"\mathbf{x}_A", 
			fontsize=24, color=:brown)
	end
	fig	
end

# ╔═╡ f93c1e5e-a3ab-4d5c-a19b-751c0dac7188
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("aug_lagrange3b.pdf", fig)
	run(`pdfcrop  --clip aug_lagrange3b.pdf aug_lagrange3.pdf`)
end

# ╔═╡ Cell order:
# ╠═0710492c-9759-11ef-2c65-f3dd16b09457
# ╠═e332a6d0-fcaa-4d62-a6e6-0135461e1d22
# ╠═e67d5422-b549-4a3a-a3b0-cd13d424243e
# ╠═5c7b7401-97f7-42b0-8fd4-ec610f51de63
# ╠═f93c1e5e-a3ab-4d5c-a19b-751c0dac7188
