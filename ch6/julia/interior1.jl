### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 3d2c5314-9b98-11ef-35c7-151aa31823ed
using Pkg; Pkg.activate()

# ╔═╡ ae748eb1-16c3-4a4c-beb3-41319dce2b44
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using PlutoUI, Printf
	using Optim, LinearAlgebra
	using ForwardDiff
	#using Nonconvex
	#Nonconvex.@load Semidefinite Ipopt
end

# ╔═╡ bfc1ae1c-a9d2-480e-a412-05a539cd3895
begin
	g1 = x -> x
	inv_b = x -> -1/g1(x)
	log_b = x -> x <= 0 ?  -log(-x) : 0
end

# ╔═╡ 2c1d34f1-f771-40c6-bd8c-b2ee0ba2705e
begin
	xs1 = -3:0.01:0
	
	fig1 = Figure(size=(400,400), fontsize=18)
	ax1 = CairoMakie.Axis(fig1[1,1], ylabel=L"\pi(\mathbf{x})", 
		xlabel=L"g(\mathbf{x})")
	limits!(ax1, -3, 0.1, -4, 8)

	lines!(ax1, xs1, inv_b, label=L"$$ Inverse barrier")
	lines!(ax1, xs1, log_b, label=L"$$ Logarithmic barrier")

	axislegend(ax1, position=:rb)

	fig1
end

# ╔═╡ d7c05d36-7123-4ae7-92c6-cda1cd325432
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("interior_penal1b.pdf", fig1)
	run(`pdfcrop  --clip interior_penal1b.pdf interior_penal1.pdf`)
end

# ╔═╡ 379f08c2-4b8d-429e-8272-a1d14e8b44f8
begin
	f = x -> x[1] + 2x[2]
	g = x -> (1/4)*x[1]^2 + x[2]^2 - 1
	Qb = (x, μ) -> f(x) + (μ)*log_b(g(x)) # need to have 0im

	# for band
	gcb = x -> sqrt(-(1/4)x^2 + 1)
end

# ╔═╡ 630e5707-6f06-4cbb-bfef-c636a68b3a03
begin
	N = 40
	cs = 2; level1 = -20:cs:20			# contour space
	xopt = [-1.41434,  -0.707164]
	x1s = -4:0.01:4; x2s = -4:0.01:4
	
	
	fc = [f([x1, x2]) for x1 in x1s, x2 in x2s]
	gc = [g([x1, x2]) for x1 in x1s, x2 in x2s]

	fig = Figure(size=(800,600), fontsize=20)
	μar = [3, 1, 0.2, 0.01]
	for (ii, μ) in enumerate(μar)
		Qbc = [Qb([x1, x2], μ) for x1 in x1s, x2 in x2s]
		global xs_P = zeros(N,2)
		global x0 = x_start = [1, 1.5] 
		xs_P[1,:] = x0
		ρ = 0.9
		
		global μ1 = μ_start = μ
		global i = 1
		
		obj = x -> f(x) + (μ1)*log_b(g(x))
		global x1 = optimize(obj, x0, Newton()).minimizer
		#x1 = x0

		while (norm(ForwardDiff.gradient(obj, x1)) > 1e-10) && (i < N) 
			#obj = x -> f(x) + (μ1)*log_b(g(x))
			x1 = optimize(obj, x0, Newton()).minimizer
			xs_P[i+1,:] = x1
			x0 = x1
			μ1 = ρ * μ1
			i = i + 1
		end
		
		println(xs_P[1:i-1,:])
		if μ == μar[1]
			ax1 = CairoMakie.Axis(fig[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L" \mu = %$(round(μ, digits=2))")
		elseif μ == μar[2] 
			ax1 = CairoMakie.Axis(fig[1,2], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L" \mu = %$(round(μ, digits=2))")
		elseif μ == μar[3]
			ax1 = CairoMakie.Axis(fig[2,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L" \mu = %$(round(μ, digits=2))")
		elseif μ == μar[4]
			ax1 = CairoMakie.Axis(fig[2,2], xlabel=L"$x_1$", ylabel=L"$x_2$", 
				title=L" \mu = %$(round(μ, digits=2))")
		end
		

		Qbc1 = [Qb([x1, x2],  μ1) for x1 in x1s, x2 in x2s]
	
		limits!(ax1, -3, 3, -2.5, 2.5)
		hidedecorations!(ax1, ticklabels=false, ticks=false, grid=true, label=false)
	
		x1b = -2:0.1:2; x1bl = -4:0.1:-2; x1br = 2:0.1:4
		band!(ax1, x1b, gcb.(x1b), 2.5*ones(length(x1b)), color=(:cyan, 0.4))
		band!(ax1, x1b, -2.5*ones(length(x1b)), -gcb.(x1b), color=(:cyan, 0.4))
		band!(ax1, x1bl, -2.5*ones(length(x1bl)), 2.5*ones(length(x1bl)), 
			color=(:cyan, 0.4))
		band!(ax1, x1br, -2.5*ones(length(x1br)), 2.5*ones(length(x1br)), 
			color=(:cyan, 0.4))
	
		contour!(ax1, x1s, x2s, fc, levels=-100:2:100, linewidth=2, linestyle=:dash)
		contour!(ax1, x1s, x2s, gc, levels=[0], color=:red, linewidth=2)
		scatter!(ax1, xopt[1], xopt[2], markersize = 15, color = :blue, 	
			strokewidth=1, strokecolor=:black)
		text!(ax1, xopt[1]+0.0, xopt[2]+0.3, text=L"\mathbf{x}^\ast", fontsize=24, 
			color=:black)
		

	
		xi = [1,1.5]
		text!(ax1, xi[1] + 0.1, xi[2], text=L"\mathbf{x}_0", fontsize=24, 
			color=:black)
		contour!(ax1, x1s, x2s, Qbc, levels=-100:2:100, color = :red,
			linewidth=1)
		#scatterlines!(ax1, xs_P[1:i,1], xs_P[1:i,2], color=:red, strokewidth=1, 
		#	strokecolor=:black, markersize=15)
		#text!(ax1, xs_P[i,1]+0.1, xs_P[i,2]-0.6, text=L"\mathbf{x}_A", 
		#	fontsize=24, color=:brown)
	end
	fig	
end

# ╔═╡ e3548ef3-d8a3-42a3-8cb2-de241df7b33e
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("interior_penal2b.pdf", fig)
	run(`pdfcrop  --clip interior_penal2b.pdf interior_penal2.pdf`)
end

# ╔═╡ Cell order:
# ╠═3d2c5314-9b98-11ef-35c7-151aa31823ed
# ╠═ae748eb1-16c3-4a4c-beb3-41319dce2b44
# ╠═bfc1ae1c-a9d2-480e-a412-05a539cd3895
# ╠═2c1d34f1-f771-40c6-bd8c-b2ee0ba2705e
# ╠═d7c05d36-7123-4ae7-92c6-cda1cd325432
# ╠═379f08c2-4b8d-429e-8272-a1d14e8b44f8
# ╠═630e5707-6f06-4cbb-bfef-c636a68b3a03
# ╠═e3548ef3-d8a3-42a3-8cb2-de241df7b33e
