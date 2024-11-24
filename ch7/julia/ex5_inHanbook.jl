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

# ╔═╡ 7c7d889c-a3df-11ef-1699-f99f594e0205
using Pkg; Pkg.activate()

# ╔═╡ 90ddcb7c-d784-4ac5-82bb-02bb71196fc8
begin
	using CairoMakie; CairoMakie.activate!()
	using LaTeXStrings, PlutoUI

	
	# Path to the Latin Modern font files
	font_path = "/usr/share/fonts/opentype/latex-xft-fonts/latinmodern-math.otf"

	set_theme!(theme_latexfonts(), fontsize=20)


end

# ╔═╡ b171beff-1510-4320-baea-cc2e92aaa229
begin
	using ForwardDiff, LinearAlgebra
	using Symbolics
	using JuMP, Ipopt
end

# ╔═╡ 35be0489-ca00-4f5c-a604-7697f19740d4
begin
	f =  (x1, x2) -> -x1 - x2
	g1 = (x1, x2) -> x1^2 - x2
	g2 = (x1, x2) -> x1^2 + x2^2 - 1
	#g1 = (x1, x2) -> x1^2 + (x2-1)^2 - 1
end

# ╔═╡ 7681215a-728f-4001-9af3-8c485694fc8e
begin
	#@variables x1, x2, λ1, λ2
	
	ff = (x) -> -x[1] - x[2]
	gg = (x) ->  [x[1]^2 - x[2], x[1]^2 + x[2]^2 - 1]
	gg1 = (x) -> gg(x)[1]
	gg2 = (x) -> gg(x)[2]
	ggv = [gg1; gg2]
	
end

# ╔═╡ 708cbab8-a258-4adb-9f74-47fc1466139d
md"""
k = $(@bind k PlutoUI.Slider(1:10, show_value = true, default = 0))
"""

# ╔═╡ 93604d82-3a0e-410b-b161-86e5d1dbe9ee
begin
	#k = 10
	x = zeros(2,k); 
	λ = zeros(2,k)
	d_opt = zeros(2,k) 

	x[:,1] = [0.5, 1]; λ[:,1] = [0, 0]

	for i = 2:k
		#global i = i+1
		∇f = ForwardDiff.gradient(ff, x[:,i-1])
		J = ForwardDiff.jacobian(gg, x[:,i-1])
		
		Hg1 = ForwardDiff.hessian(gg1, x[:,i-1])
		Hg2 = ForwardDiff.hessian(gg2, x[:,i-1])
		B = ForwardDiff.hessian(ff, x[:,i-1]) + λ[:,i-1]' * [Hg1, Hg2]
		 p1 = gg1(x[:,i-1]); p2 = gg2(x[:,i-1]);
		∇p1 = ForwardDiff.gradient(gg1, x[:,i-1]); 
		∇p2 = ForwardDiff.gradient(gg2, x[:,i-1])

		model = Model(Ipopt.Optimizer)
		set_silent(model)
		@variable(model, d[1:2])
		@objective(model, Min, (0.5)*d'*B*d + ∇f'*d)	
		con1 = @constraint(model, p1 .+ ∇p1'*d <= 0)
		con2 = @constraint(model, p2 .+ ∇p2'*d <= 0)
		set_start_value(d[1], x[1,i-1])
		set_start_value(d[2], x[2,i-1])
		optimize!(model); 
		
		d_opt[:,i] = value.(d);
		
		x[:,i] = x[:,i-1] + d_opt[:,i]
		λ[:,i] = -[dual(con1); dual(con2)];
	end
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
	
	scatterlines!(ax1, x[1,:], x[2,:], markersize=15, color=:blue, linestyle=:dash)
	scatter!(ax1, 1/sqrt(2), 1/sqrt(2), markersize=12, color=:red, 
		strokecolor=:black, strokewidth=1)
	text!(ax1, 0.52, 0.6, text=L"\mathbf{x}^\ast")
	text!(ax1, x[1,1], x[2,1], text=L"\mathbf{x}_0")
	fig1
end

# ╔═╡ 8d6abc20-b1c1-4345-ba46-381d21dce36a
x

# ╔═╡ 04f88027-fa90-48df-9b56-452a0a05a627
λ

# ╔═╡ 57bda8fc-8750-4d3a-9371-194c5f13dcde
d_opt

# ╔═╡ 9162ffab-1f53-4ace-8dfb-d86c861e0503
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/ex5_han.png", fig1)
end
  ╠═╡ =#

# ╔═╡ 08018a43-3350-411b-981b-92c60ef00d74


# ╔═╡ 779db52c-b0fe-4a0d-83f4-4f51c39c3f0c
function bfgs_update(H, s, y)
    rho = 1.0 / (y' * s)
    I = Matrix{Float64}(I, length(s), length(s))
    V = I - rho * s * y'
    H = V' * H * V + rho * s * s'
    return H
end

# ╔═╡ 5ec86475-3c5a-4cb2-aad5-f0578e25c8fc
function sqp(obj_func, gradient, x0; max_iter=100, tol=1e-6)
    x = x0
    H = Matrix{Float64}(I, length(x0), length(x0))  # Initial Hessian approximation
    for iter in 1:max_iter
        grad = gradient(x)
        p = -H \ grad  # Search direction
        alpha = 1.0  # Step length (can be optimized further)
        x_new = x + alpha * p
        s = x_new - x
        y = gradient(x_new) - grad
        H = bfgs_update(H, s, y)
        x = x_new
        if norm(grad) < tol
            break
        end
    end
    return x
end

# ╔═╡ Cell order:
# ╠═7c7d889c-a3df-11ef-1699-f99f594e0205
# ╠═90ddcb7c-d784-4ac5-82bb-02bb71196fc8
# ╠═b171beff-1510-4320-baea-cc2e92aaa229
# ╠═35be0489-ca00-4f5c-a604-7697f19740d4
# ╟─84612895-9a53-405d-b1c7-f965872bce0a
# ╟─7681215a-728f-4001-9af3-8c485694fc8e
# ╠═93604d82-3a0e-410b-b161-86e5d1dbe9ee
# ╠═708cbab8-a258-4adb-9f74-47fc1466139d
# ╠═8d6abc20-b1c1-4345-ba46-381d21dce36a
# ╠═04f88027-fa90-48df-9b56-452a0a05a627
# ╠═57bda8fc-8750-4d3a-9371-194c5f13dcde
# ╠═9162ffab-1f53-4ace-8dfb-d86c861e0503
# ╠═5ec86475-3c5a-4cb2-aad5-f0578e25c8fc
# ╠═08018a43-3350-411b-981b-92c60ef00d74
# ╠═779db52c-b0fe-4a0d-83f4-4f51c39c3f0c
