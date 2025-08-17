### A Pluto.jl notebook ###
# v0.20.16

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

# ╔═╡ 5e5879b2-79fc-11f0-2893-c3bb752fc2b3
using Pkg; Pkg.activate()

# ╔═╡ 58ee1144-8078-4ad6-ab1a-e4097cf589d3
begin
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using ForwardDiff
	using PlutoUI
end

# ╔═╡ 11a0dadd-a843-4113-85c3-4e2244d807e9
begin
	f(x) = 1/5*(x[1]^2 + x[2]^2) + 1
	fig1 = Figure(size=(1000,800))
	empty!(fig1)
end

# ╔═╡ a69dac68-ded1-4a25-b69f-fc195f663402
begin
	md"""
	Elevator = $(@bind e PlutoUI.Slider(π/9:0.1:π/2, show_value=true, default=π/8))

	Perspectiveness = $(@bind p PlutoUI.Slider(0:0.1:1, show_value=true, default=0.5))

	Azimuth = $(@bind a PlutoUI.Slider(π/9:0.1:π/2, show_value=true, default=π/9))
	"""
end

# ╔═╡ bc47126b-c6b5-45d5-94c0-5bace05e1988
function plot_surface1()
	ax1 = Axis3(fig1[1,1], aspect = (1,1,1), perspectiveness = p, 
			elevation = e, azimuth = a, xspinesvisible=false, yspinesvisible=false,
			zspinesvisible=false, xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false, xgridvisible = false,
			ygridvisible = false, zgridvisible = false, xticksvisible=false,
			yticksvisible=false, zticksvisible=false,
			xlabel="", ylabel="", zlabel="", title=L"$ $ Positive Definite")
	
	xs = ys = -2:0.05:2
	zs = [f([x, y]) for x in xs, y in ys]	
	limits!(ax1, -2,2, -2, 2, -0.5, 2)
	
	surface!(ax1, xs, ys, zs)

	xws = -2:0.5:2; yws = -2:0.5:2;
	zws = [f([x, y]) for x in xws, y in yws]
	wireframe!(ax1, xws, yws, zws, overdraw = true, transparency = true, 
			   color=(:black, 0.7), linewidth=0.5)
	contour!(ax1, xs, ys, zs, linewidth=1, color=:blue)
	scatter!(ax1, 0, 0, 0, color=(:red, 1.0))
	text!(ax1, 0.75, -0.25, -0.3, text=L"$ $ Minimum", color=(:black, 1.0))
	fig1
end

# ╔═╡ 5d245b93-4f1f-4604-8ff5-3ff6e6150fdc
function plot_surface2()
	ax1 = Axis3(fig1[1,2], aspect = (1,1,1), perspectiveness = p, 
			elevation = e, azimuth = a, xspinesvisible=false, yspinesvisible=false,
			zspinesvisible=false, xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false, xgridvisible = false,
			ygridvisible = false, zgridvisible = false, xticksvisible=false,
			yticksvisible=false, zticksvisible=false,
			xlabel="", ylabel="", zlabel="", title=L"$ $ Negative Definite")
	
	
	xs = ys = -2:0.05:2
	zs = [-f([x, y]) for x in xs, y in ys]	
	limits!(ax1, -2,2, -2, 2, -2, 0.3)
	
	surface!(ax1, xs, ys, zs)

	xws = -2:0.5:2; yws = -2:0.5:2;
	zws = [-f([x, y]) for x in xws, y in yws]
	wireframe!(ax1, xws, yws, zws, overdraw = true, transparency = true, 
			   color=(:black, 0.7), linewidth=0.5)
	scatter!(ax1, 0, 0, 0, color=(:red, 1.0))
	contour!(ax1, xs, ys, zs, linewidth=1, color=:blue)
	text!(ax1, 0.75, -0.25, 0.2, text=L"$ $ Maximum", color=(:black, 1.0))
	fig1
end

# ╔═╡ 3797ab71-c9ef-4b77-9545-28f22e5fb8d9
function plot_surface3()
	f(x) = 1 * (x[1]^2 + 1)
	xs = ys = -2:0.05:2
	z = [f([x, y]) for x in xs, y in ys]	
	
	
	ax1 = Axis3(fig1[2,1], aspect = (1,1,1), perspectiveness = p, 
			elevation = e, azimuth = a*3, xspinesvisible=false, yspinesvisible=false,
			zspinesvisible=false, xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false, xgridvisible = false,
			ygridvisible = false, zgridvisible = false, xticksvisible=false,
			yticksvisible=false, zticksvisible=false,
			xlabel="", ylabel="", zlabel="", title=L"$ $ Positive Semidefinite")
	
	surface!(ax1, xs, xs, z, color=:blue, colorrange=(-20,10))
	contour!(ax1, xs, xs, z, levels= -100:0.2:100, color=(:blue, 0.7), linewidth =1,
		transformation = (:xy, -0.5))
	contour!(ax1, xs, xs, z, levels= -1:0.1:1, color=(:red, 0.7), linewidth =1,
		transformation = (:xy, -0.5))
		
	text!(ax1, 1, -0.25, -2.5, text=L"$ $ Weak Minimum Line", color=(:black, 1.0))
	x1qw = -2:0.5:2; x2qw = -2:0.5:2;
	zqw3 = [f([x, y]) for x in x1qw, y in x2qw]
	wireframe!(ax1, x1qw, x2qw, zqw3, overdraw = true, transparency = true,
	        color = (:black, 0.7), linewidth=0.5)
	limits!(ax1, -2,2,-2, 2, -3, 4)
	fig1
end

# ╔═╡ 648d6d39-5116-4868-a84f-69438e1d3098
function plot_surface4()
	f(x) = 1/5 * (x[1]^2 - x[2]^2)
	xs = ys = -2:0.05:2
	z = [f([x, y]) for x in xs, y in ys]
	
	ax1 = Axis3(fig1[2,2], aspect = (1,1,1), perspectiveness = p, 
			elevation = e, azimuth = a*2.5, xspinesvisible=false, yspinesvisible=false,
			zspinesvisible=false, xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false, xgridvisible = false,
			ygridvisible = false, zgridvisible = false, xticksvisible=false,
			yticksvisible=false, zticksvisible=false,
			xlabel="", ylabel="", zlabel="", title=L"$ $ Indefinite")
	surface!(ax1, xs, ys, z, color=(:blue, 1), colorrange = (-30, 20))
	contour!(ax1, xs, ys, z, levels= -100:0.2:100, color=(:blue, 0.7), linewidth =1,
		transformation = (:xy, -1))

	scatter!(ax1, 0, 0, -1, color=(:red, 1.0))
	
	text!(ax1, 1, -0.25, -2.25, text=L"$ $ Saddle Point", color=(:red, 1.0))
	
	xws = yws = -2:0.5:2
	zw = [f([x, y]) for x in xws, y in yws]
	wireframe!(ax1, xws, yws, zw, overdraw = true, transparency = true,
	        color = (:black, 0.7), linewidth=0.5)
	
	limits!(ax1, -2,2,-2, 2, -2, 0.4)
	text!(ax1, 0.75, -0.25, -1.5, text=L"$ $ Saddle Point", color=(:black, 1.0))
	fig1
end

# ╔═╡ 503ca806-ee9a-4f12-a246-b7ff1dfd874f
begin
	empty!(fig1)
	plot_surface1()
	plot_surface2()
	plot_surface3()
	plot_surface4()
end

# ╔═╡ e05ae4f4-7613-4023-beab-ae93f7827892
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/figures/condition_h.png", fig1)
# don't save as PDF, the file is too complicated.

# ╔═╡ Cell order:
# ╠═5e5879b2-79fc-11f0-2893-c3bb752fc2b3
# ╠═58ee1144-8078-4ad6-ab1a-e4097cf589d3
# ╠═11a0dadd-a843-4113-85c3-4e2244d807e9
# ╟─a69dac68-ded1-4a25-b69f-fc195f663402
# ╠═bc47126b-c6b5-45d5-94c0-5bace05e1988
# ╠═5d245b93-4f1f-4604-8ff5-3ff6e6150fdc
# ╟─3797ab71-c9ef-4b77-9545-28f22e5fb8d9
# ╠═648d6d39-5116-4868-a84f-69438e1d3098
# ╠═503ca806-ee9a-4f12-a246-b7ff1dfd874f
# ╠═e05ae4f4-7613-4023-beab-ae93f7827892
