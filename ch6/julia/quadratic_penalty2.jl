### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

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
	Q = (x, μ) -> f(x) + (μ/2)*(max(0,g(x)))^2

	# for band
	gcb = x -> sqrt(-(1/4)x^2 + 1)
end

# ╔═╡ 79507495-51e1-4a30-bb92-e664243b23ab
md"""
The penalty value
μ = $(@bind μ PlutoUI.Slider(0:0.5:20, show_value=true, default=1))
"""

# ╔═╡ 5c7b7401-97f7-42b0-8fd4-ec610f51de63
begin
	x1s = -4:0.01:4; x2s = -4:0.01:4
	Qc = [Q([x1, x2], μ) for x1 in x1s, x2 in x2s]
	fc = [f([x1, x2]) for x1 in x1s, x2 in x2s]
	gc = [g([x1, x2]) for x1 in x1s, x2 in x2s]
	fig = Figure(size=(800,600), fontsize=20)
	cs = 2; level1 = -20:cs:20			# contour space
	xopt = [-1.41434,  -0.707164]
	tst = @sprintf("%1f", μ )
	ax1 = CairoMakie.Axis(fig[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ")
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
	contour!(ax1, x1s, x2s, Qc, levels=level1, colormap = :viridis, 
		linewidth=1)
	contour!(ax1, x1s, x2s, gc, levels=[0], color=:red, linewidth=2)
	scatter!(ax1, xopt[1], xopt[2], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax1, -1, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)

	
	
	xi = [1,0.5]
	# the optimal point of Q
	obji = x -> f(x) + (μ/2)*(max(0,g(x)))^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax1, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax1, Xsol[1]+0.1, Xsol[2]-0.75, text=L"\mathbf{x}_Q", fontsize=24, color=:black)


	# this part for μ=4
	μ1 = 4
	Qc = [Q([x1, x2], μ1) for x1 in x1s, x2 in x2s]
	ax2 = CairoMakie.Axis(fig[1,2], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ1")
	limits!(ax2, -3, 3, -2.5, 2.5)
	hidedecorations!(ax2, ticklabels=false, ticks=false, grid=true, label=false)
	
	band!(ax2, x1b, gcb.(x1b), 2.5*ones(length(x1b)), color=(:cyan, 0.4))
	band!(ax2, x1b, -2.5*ones(length(x1b)), -gcb.(x1b), color=(:cyan, 0.4))
	band!(ax2, x1bl, -2.5*ones(length(x1bl)), 2.5*ones(length(x1bl)), 
		color=(:cyan, 0.4))
	band!(ax2, x1br, -2.5*ones(length(x1br)), 2.5*ones(length(x1br)), 
		color=(:cyan, 0.4))
	
	contour!(ax2, x1s, x2s, fc, levels=-100:1:100, linewidth=2, linestyle=:dash)
	contour!(ax2, x1s, x2s, Qc, levels=level1, colormap = :viridis)
	contour!(ax2, x1s, x2s, gc, levels=[0], color=:red, linewidth=2)
	scatter!(ax2, xopt[1], xopt[2], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax2, -1, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	# the optimal point of Q
	obji = x -> f(x) + (μ1/2)*(max(0,g(x)))^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax2, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax2, Xsol[1]+0.1, Xsol[2]-0.75, text=L"\mathbf{x}_Q", fontsize=24, 
		color=:black)
	
	# this part for μ=8
	μ1 = 8
	Qc = [Q([x1, x2], μ1) for x1 in x1s, x2 in x2s]
	ax3 = CairoMakie.Axis(fig[2,1], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ1")
	limits!(ax3, -3, 3, -2.5, 2.5)
	hidedecorations!(ax3, ticklabels=false, ticks=false, grid=true, label=false)
	
	band!(ax3, x1b, gcb.(x1b), 2.5*ones(length(x1b)), color=(:cyan, 0.4))
	band!(ax3, x1b, -2.5*ones(length(x1b)), -gcb.(x1b), color=(:cyan, 0.4))
	band!(ax3, x1bl, -2.5*ones(length(x1bl)), 2.5*ones(length(x1bl)), 
		color=(:cyan, 0.4))
	band!(ax3, x1br, -2.5*ones(length(x1br)), 2.5*ones(length(x1br)), 
		color=(:cyan, 0.4))
	
	contour!(ax3, x1s, x2s, fc, levels=-100:1:100, linewidth=2, linestyle=:dash)
	contour!(ax3, x1s, x2s, Qc, levels=level1, colormap = :viridis)
	contour!(ax3, x1s, x2s, gc, levels=[0], color=:red, linewidth=2)
	scatter!(ax3, xopt[1], xopt[2], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax3, -1, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	# the optimal point of Q
	obji = x -> f(x) + (μ1/2)*(max(0,g(x)))^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax3, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax3, Xsol[1]+0.1, Xsol[2]-0.75, text=L"\mathbf{x}_Q", fontsize=24, 
		color=:black)

	# this part for μ = 10
	μ1 = 10
	Qc = [Q([x1, x2], μ1) for x1 in x1s, x2 in x2s]
	ax4 = CairoMakie.Axis(fig[2,2], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ1")
	limits!(ax4, -3, 3, -2.5, 2.5)
	hidedecorations!(ax4, ticklabels=false, ticks=false, grid=true, label=false)
	
	band!(ax4, x1b, gcb.(x1b), 2.5*ones(length(x1b)), color=(:cyan, 0.4))
	band!(ax4, x1b, -2.5*ones(length(x1b)), -gcb.(x1b), color=(:cyan, 0.4))
	band!(ax4, x1bl, -2.5*ones(length(x1bl)), 2.5*ones(length(x1bl)), 
		color=(:cyan, 0.4))
	band!(ax4, x1br, -2.5*ones(length(x1br)), 2.5*ones(length(x1br)), 
		color=(:cyan, 0.4))
	
	contour!(ax4, x1s, x2s, fc, levels=-100:1:100, linewidth=2, linestyle=:dash)
	contour!(ax4, x1s, x2s, Qc, levels=level1, colormap = :viridis)
	contour!(ax4, x1s, x2s, gc, levels=[0], color=:red, linewidth=2)
	scatter!(ax4,xopt[1], xopt[2], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax4, -1, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	# the optimal point of Q
	obji = x -> f(x) + (μ1/2)*(max(0,g(x)))^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax4, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax4,Xsol[1]+0.1, Xsol[2]-0.75, text=L"\mathbf{x}_Q", fontsize=24, 
		color=:black)

	fig
	
end

# ╔═╡ f93c1e5e-a3ab-4d5c-a19b-751c0dac7188
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("quad_penal2b.pdf", fig)
	run(`pdfcrop  --clip quad_penal2b.pdf quad_penal2.pdf`)
end

# ╔═╡ e184f0f6-54f2-48a4-9ba7-da277b20bd25
begin
	N = 2000
	xs_P = zeros(N,2)
	x0 = x_start = [1, 0.5]
	xs_P[1,:] = x0
	λ = λ_start = 1
	γ = 1.7

	i = 2
	obj = x -> f(x) + (λ/2)*(max(0,g(x)))^2
	x1 = optimize(obj, x0).minimizer
	
	while (norm(ForwardDiff.gradient(obj, x1)) > 1e-10) && (i < N) 
		#obj = x -> f(x) + (λ/2)*h(x)^2
		global x1 = optimize(obj, x0 + randn(2)/5).minimizer
		xs_P[i,:] = x1
		global x0 = x1
		global λ *= γ
		global i += 1
	end
	
	#scatterlines!(ax2, xs_P[1:i-1,:], color=:cyan, strokewidth=1, 
	#	strokecolor=:black, markersize=15)
	
	#xs_P[1:i-1,:]
end

# ╔═╡ 6a9166f5-6692-4ba4-96dd-f22a262852d9
xs_P[1:i-1,:]

# ╔═╡ Cell order:
# ╠═0710492c-9759-11ef-2c65-f3dd16b09457
# ╠═e332a6d0-fcaa-4d62-a6e6-0135461e1d22
# ╠═e67d5422-b549-4a3a-a3b0-cd13d424243e
# ╠═79507495-51e1-4a30-bb92-e664243b23ab
# ╠═5c7b7401-97f7-42b0-8fd4-ec610f51de63
# ╠═f93c1e5e-a3ab-4d5c-a19b-751c0dac7188
# ╠═e184f0f6-54f2-48a4-9ba7-da277b20bd25
# ╠═6a9166f5-6692-4ba4-96dd-f22a262852d9
