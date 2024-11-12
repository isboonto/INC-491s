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
	h = x -> (1/4)x[1]^2 + x[2]^2 - 1
	Q = (x, μ) -> f(x) + (μ/2)*(h(x).^2)
	#Q = (x, μ) -> f(x) + (μ)*norm(h(x),1)
	Ql1 = (x, μ) -> f(x) + (μ)*norm(h(x),1)
end

# ╔═╡ 79507495-51e1-4a30-bb92-e664243b23ab
md"""
The penalty value
μ = $(@bind μ PlutoUI.Slider(0:0.5:10, show_value=true, default=1))
"""

# ╔═╡ 5c7b7401-97f7-42b0-8fd4-ec610f51de63
begin
	x1s = -4:0.01:4; x2s = -4:0.01:4
	Qc = [Q([x1, x2], μ) for x1 in x1s, x2 in x2s]
	fc = [f([x1, x2]) for x1 in x1s, x2 in x2s]
	hc = [h([x1, x2]) for x1 in x1s, x2 in x2s]
	
	fig = Figure(size=(800,700), fontsize=20)

	ax1 = CairoMakie.Axis(fig[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ")
	limits!(ax1, -3, 3, -3, 3)
	hidedecorations!(ax1, ticklabels=false, ticks=false, grid=true, label=false)
	
	contour!(ax1, x1s, x2s, fc, levels=-100:2:100, linewidth=2, linestyle=:dash)
	contour!(ax1, x1s, x2s, Qc, levels=-30:2:30, colormap = :viridis, 
		linewidth=1)
	contour!(ax1, x1s, x2s, hc, levels=[0], color=:red, linewidth=2)
	scatter!(ax1,  [-1.414], [-0.707], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax1, -1.4, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	
	xi = [1,0.5] # initial point
	# the optimal point of Q
	obji = x -> f(x) + (μ/2)*h(x)^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax1, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax1, -2.1, -1, text=L"\mathbf{x}_Q", fontsize=24, color=:black)


	# this part for μ=4
	μ1 = 4
	Qc = [Q([x1, x2], μ1) for x1 in x1s, x2 in x2s]
	ax2 = CairoMakie.Axis(fig[1,2], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ1")
	limits!(ax2, -3, 3, -3, 3)
	hidedecorations!(ax2, ticklabels=false, ticks=false, grid=true, label=false)
	
	contour!(ax2, x1s, x2s, fc, levels=-100:2:100, linewidth=2, linestyle=:dash)
	contour!(ax2, x1s, x2s, Qc, levels=-30:2:30, colormap = :viridis)
	contour!(ax2, x1s, x2s, hc, levels=[0], color=:red, linewidth=2)
	scatter!(ax2,  [-1.414], [-0.707], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax2, -1.4, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	# the optimal point of Q
	obji = x -> f(x) + (μ1/2)*h(x)^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax2, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax2, -1.4, -1.4, text=L"\mathbf{x}_Q", fontsize=24, color=:black)
	
	# this part for μ=8
	μ1 = 8
	Qc = [Q([x1, x2], μ1) for x1 in x1s, x2 in x2s]
	ax3 = CairoMakie.Axis(fig[2,1], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ1")
	limits!(ax3, -3, 3, -3, 3)
	hidedecorations!(ax3, ticklabels=false, ticks=false, grid=true, label=false)
	
	contour!(ax3, x1s, x2s, fc, levels=-100:2:100, linewidth=2, linestyle=:dash)
	contour!(ax3, x1s, x2s, Qc, levels=-20:2:20, colormap = :viridis)
	contour!(ax3, x1s, x2s, hc, levels=[0], color=:red, linewidth=2)
	scatter!(ax3,  [-1.414], [-0.707], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax3, -1.4, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	# the optimal point of Q
	obji = x -> f(x) + (μ1/2)*h(x)^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax3, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax3, -1.4, -1.4, text=L"\mathbf{x}_Q", fontsize=24, color=:black)

	# this part for μ = 10
	μ1 = 10
	Qc = [Q([x1, x2], μ1) for x1 in x1s, x2 in x2s]
	ax4 = CairoMakie.Axis(fig[2,2], xlabel=L"$x_1$", ylabel=L"$x_2$", title=L"$\mu$ = %$μ1")
	limits!(ax4, -3, 3, -3, 3)
	hidedecorations!(ax4, ticklabels=false, ticks=false, grid=true, label=false)
	
	contour!(ax4, x1s, x2s, fc, levels=-100:2:100, linewidth=2, linestyle=:dash)
	contour!(ax4, x1s, x2s, Qc, levels=-20:2:20, colormap = :viridis)
	contour!(ax4, x1s, x2s, hc, levels=[0], color=:red, linewidth=2)
	scatter!(ax4, [-1.414], [-0.707], markersize = 10, color = :blue, strokewidth=1, 
		strokecolor=:black)
	text!(ax4, -1.4, -0.9, text=L"\mathbf{x}^\ast", fontsize=24, color=:black)
	# the optimal point of Q
	obji = x -> f(x) + (μ1/2)*h(x)^2
	Xsol = optimize(obji, xi).minimizer
	scatter!(ax4, Xsol[1], Xsol[2], markersize = 10, color = :red, strokewidth=1, 
		strokecolor=:black)
	text!(ax4, -1.4, -1.4, text=L"\mathbf{x}_Q", fontsize=24, color=:black)

	fig
	
end

# ╔═╡ f93c1e5e-a3ab-4d5c-a19b-751c0dac7188
# ╠═╡ disabled = true
#=╠═╡
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("quad_penal1b.pdf", fig)
	run(`pdfcrop  --clip quad_penal1b.pdf quad_penal1.pdf`)
end
  ╠═╡ =#

# ╔═╡ 5ad31383-2d68-42f9-9e8f-10cd1e0b9fc3
md"""
The penalty value
μ = $(@bind μ_itd PlutoUI.Slider(0:0.5:10, show_value=true, default=1))
"""

# ╔═╡ d1e45dd2-fa0e-4152-91a3-a12bfdb6e1cd
begin
	N = 20
	xs_P = zeros(N,2)
	x0 = x_start = [1.7, 2.0]
	xs_P[1,:] = x0
	μ_it = μ_itd
	γ = 1.7
	Qc1 = [Q([x1, x2], μ_it) for x1 in x1s, x2 in x2s]
	fc1 = [f([x1, x2]) for x1 in x1s, x2 in x2s]
	hc1 = [h([x1, x2]) for x1 in x1s, x2 in x2s]

	fig7 = Figure(size=(600,400), fontsize=18)
	ax7 = CairoMakie.Axis(fig7[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
		title=L"$\mu$ = %$μ_it")
	limits!(ax7, -3, 3, -3, 3)
	
	i = 2
	obj = x -> f(x) + (μ_it/2)*h(x)^2
	x1 = optimize(obj, x0).minimizer
	
	while (norm(ForwardDiff.gradient(obj, x1)) > 1e-6) && (i < 20) 
		global x1 = optimize(obj, x0).minimizer
		xs_P[i,:] = x1
		global x0 = x1
		global μ_it *= γ
		global i += 1
	end

	

	contour!(ax7, x1s, x2s, fc1, levels=-100:1:100, linewidth=2, linestyle=:dash)
	contour!(ax7, x1s, x2s, Qc1, levels=-80:3:80, colormap = :viridis)
	contour!(ax7, x1s, x2s, hc1, levels=[0], color=:red, linewidth=2)
	scatterlines!(ax7, xs_P[1:i-1,:], color=:blue, strokewidth=1, 
			strokecolor=:black, markersize=15, markercolor=:cyan)
	fig7
end

# ╔═╡ a0ed3eb0-49ab-48ee-b8a5-8105dcb51f8d
xs_P

# ╔═╡ 534e0019-ae40-42c3-a1e1-1644abcdca71
xs_P

# ╔═╡ 1125e36e-3295-4217-8b13-4ed40801acd9
# ╠═╡ disabled = true
#=╠═╡
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("quad_penal_iter1b.pdf", fig7)
	run(`pdfcrop  --clip quad_penal_iter1b.pdf quad_penal_iter1.pdf`)
end
  ╠═╡ =#

# ╔═╡ 15bed972-ff53-4bcc-ac8f-f97663c939a6
begin
	xs_Pl1 = zeros(N,2)
	xs_Pl1[1,:] = x0
	μ_itl1 = μ_itd
	Qcl1 = [Ql1([x1, x2], μ_itl1) for x1 in x1s, x2 in x2s]
	#fc1 = [f([x1, x2]) for x1 in x1s, x2 in x2s]
	#hc1 = [h([x1, x2]) for x1 in x1s, x2 in x2s]
	
	fig8 = Figure(size=(600,400), fontsize=18)
	ax8 = CairoMakie.Axis(fig8[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", 
		title=L"$\mu$ = %$μ_itl1")
	limits!(ax8, -2, 2, -2, 2)
	
	il1 = 2
	x0l1 = [1.7, 0.0]
	objl1 = x -> f(x) + (μ_it)*norm(h(x),1)
	x1l1 = optimize(objl1, x0l1).minimizer
	
	while (norm(ForwardDiff.gradient(objl1, x1l1)) > 1e-6) && (il1 < 20) 
		global x1l1 = optimize(objl1, x0l1).minimizer
		xs_Pl1[i,:] = x1l1
		global x0l1 = x1l1
		global μ_itl1 *= γ
		global il1 += 1
	end

	

	contour!(ax8, x1s, x2s, fc1, levels=-100:1:100, linewidth=2, linestyle=:dash)
	contour!(ax8, x1s, x2s, Qcl1, levels=-80:3:80, colormap = :viridis)
	contour!(ax8, x1s, x2s, hc1, levels=[0], color=:red, linewidth=2)
	scatterlines!(ax8, xs_Pl1[1:i-1,:], color=:blue, strokewidth=1, 
			strokecolor=:black, markersize=15, markercolor=:cyan)
	fig8
end

# ╔═╡ Cell order:
# ╠═0710492c-9759-11ef-2c65-f3dd16b09457
# ╠═e332a6d0-fcaa-4d62-a6e6-0135461e1d22
# ╠═e67d5422-b549-4a3a-a3b0-cd13d424243e
# ╠═79507495-51e1-4a30-bb92-e664243b23ab
# ╠═5c7b7401-97f7-42b0-8fd4-ec610f51de63
# ╠═a0ed3eb0-49ab-48ee-b8a5-8105dcb51f8d
# ╠═f93c1e5e-a3ab-4d5c-a19b-751c0dac7188
# ╟─5ad31383-2d68-42f9-9e8f-10cd1e0b9fc3
# ╠═d1e45dd2-fa0e-4152-91a3-a12bfdb6e1cd
# ╠═534e0019-ae40-42c3-a1e1-1644abcdca71
# ╠═1125e36e-3295-4217-8b13-4ed40801acd9
# ╠═15bed972-ff53-4bcc-ac8f-f97663c939a6
