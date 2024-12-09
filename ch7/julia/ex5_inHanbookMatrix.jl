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
	using LaTeXStrings
	
	# Path to the Latin Modern font files
	font_path = "/usr/share/fonts/opentype/latex-xft-fonts/latinmodern-math.otf"

	set_theme!(theme_latexfonts(), fontsize=20)


end

# ╔═╡ b171beff-1510-4320-baea-cc2e92aaa229
begin
	using ForwardDiff, LinearAlgebra
	using Symbolics
	using JuMP, Ipopt
	using PlutoUI
end

# ╔═╡ 35be0489-ca00-4f5c-a604-7697f19740d4
begin
	#f =  (x1, x2) -> (1-x1)^2 + (1-x2)^2 + 0.5*(2x2 - x1^2)^2
	f = (x1, x2) -> -x1 - x2
	#g1 = (x1, x2) -> 2x1 + 3x2 + 1  # -2x1 - 3x2 - 1 <= 0
	g1 = (x1, x2) -> x1^2 - x2 
	#g2 = (x1, x2) -> (1/4)*x1^2 + x2^2 - 1
	g2 = (x1, x2) -> x1^2 + x2^2 - 1
end

# ╔═╡ 5f4eabf2-6949-42ee-b458-1e76742f4b47
md"""
 $x_1$ = $(@bind xv1 PlutoUI.Slider(-2.5:0.1:2.5, default = -1, show_value=true))

 $x_2$ = $(@bind xv2 PlutoUI.Slider(-2.5:0.1:2.5, default = 2, show_value=true))
"""

# ╔═╡ 7681215a-728f-4001-9af3-8c485694fc8e
begin
	xx = [0,0]
	#ff = (x) -> (1-x[1])^2 + (1-x[2])^2 + 0.5*(2x[2] - x[1]^2)^2
	ff = (x) -> -x[1] - x[2]
	#gg = (x) -> [2x[1] + 3x[2] + 1, (1/4)*x[1]^2 + x[2]^2 - 1 ]
	gg = (x) -> [x[1]^2 - x[2] , x[1]^2 + x[2]^2 - 1]
	gg1 = (x) -> gg(x)[1]
	gg2 = (x) -> gg(x)[2]
	∇gg = x -> ForwardDiff.jacobian(gg, x)
	∇ff = x -> ForwardDiff.gradient(ff, x)
end

# ╔═╡ 332e4e2a-8938-49ac-9e89-9bee128a34b6
ff([1,2])

# ╔═╡ 0d246343-119f-4e93-ad79-66a5c60f73cb
begin
	gradient_ff = x -> ForwardDiff.gradient(ff, x)
	gradient_gg1 = x -> ForwardDiff.gradient(gg1, x)
	gradient_gg2 = x -> ForwardDiff.gradient(gg2, x)
end

# ╔═╡ 30a1b1d5-5e91-4c38-87a2-109b46745d47
function bfgs_update(H, s, y)
    rho = 1.0 / (y' * s)
    #I = Matrix{Float64}(I, length(s), length(s))
    V = I - rho * s * y'
    H = V' * H * V + rho * s * s'
    return H
end

# ╔═╡ 93604d82-3a0e-410b-b161-86e5d1dbe9ee
begin
	k = 30; 
	x = zeros(2,k); λ = zeros(2,k); d_opt = zeros(2,k) 
	x[:,1] = [xv1, xv2]; λ[:,1] = [1, 1]

	for i = 2:k
		∇f = ForwardDiff.gradient(ff, x[:,i-1])
		J = ForwardDiff.jacobian(gg, x[:,i-1])
		
		Hg1 = ForwardDiff.hessian(gg1, x[:,i-1])
		Hg2 = ForwardDiff.hessian(gg2, x[:,i-1])
		
		B = ForwardDiff.hessian(ff, x[:,i-1]) + λ[:,i-1]' * [Hg1, Hg2]
		
		
		p1 = gg1(x[:,i-1]); 
		p2 = gg2(x[:,i-1]);
		∇p1 = ForwardDiff.gradient(gg1, x[:,i-1]); 
		∇p2 = ForwardDiff.gradient(gg2, x[:,i-1])
		
		
		
		if i == 1
			if !isposdef(B)
				B = I
			end
		else 
			if !isposdef(B)
				s = x[:,i] .- x[:,i-1]
				y = ForwardDiff.gradient(ff, x[:,i]) .- ForwardDiff.gradient(ff, x[:,i-1])
				B = bfgs_update(B, s, y)
			end
		end
		
		@show BH = [B J'; J zeros(size(J,1), size(J,2))]
	
		bb = [∇f + J'*λ[:,i-1] ; -gg(x[:,i-1])];
		res = BH\bb;
		
			
		x[:,i] = x[:,i-1] + res[1:2,1];
		λ[:,i] = res[3:4,1];
	end
end

# ╔═╡ 84612895-9a53-405d-b1c7-f965872bce0a
begin
	fig1 = Figure(size=(600,400), fontsize=20)
	ax1 = Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", aspect=1.4)
	limits!(ax1, -2, 2, -2, 2)

	xs1 = -5:0.02:5; 
	contour!(ax1, xs1, xs1, f, levels=-200:4:300)
	contour!(ax1, xs1, xs1, f, levels=-3:1:3)
	contour!(ax1, xs1, xs1, g1, levels=[0, 0], linewidth=2, color=:red)
	contour!(ax1, xs1, xs1, g2, levels=[0], linewidth=2, color=:blue)

	xb = -1.5:0.05:0.89
	xb2 = -0.8:0.05:0.8
	#g1b = (x) -> -(1/3) - (2/3)x
	g1b = (x) -> x[1]^2
	#g2b = (x) -> sqrt(1-(1/4)*x^2)
	g2b = (x) -> sqrt(1 - x[1]^2)
	band!(ax1, xb2, g1b.(xb2), g2b.(xb2), color=(:orange, 0.4))
	#band!(ax1, xb2, -g2b.(xb2), g2b.(xb2), color=(:orange, 0.4))
	
	scatterlines!(ax1, x[1,:], x[2,:], markersize=15, color=:blue, linestyle=:dash)
	scatter!(ax1, 0.707, 0.707, markersize=12, color=:red, 
		strokecolor=:black, strokewidth=1)
	text!(ax1, x[1,end], x[2,end]+0.1, text=L"\mathbf{x}^\ast")
	text!(ax1, x[1,1], x[2,1], text=L"\mathbf{x}_0")
	fig1
end

# ╔═╡ 9162ffab-1f53-4ace-8dfb-d86c861e0503
# ╠═╡ disabled = true
#=╠═╡
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("sqp_1b.pdf", fig1)
	
	run(`pdfcrop  --clip sqp_1b.pdf sqp_1.pdf`)
	run(`rm sqp_1b.pdf`)
end
  ╠═╡ =#

# ╔═╡ 8d6abc20-b1c1-4345-ba46-381d21dce36a
x

# ╔═╡ 04f88027-fa90-48df-9b56-452a0a05a627
λ

# ╔═╡ c84bc9ce-77f9-41e1-aed0-3a65becb0e25
begin
	function sqp(ff, gradient_ff, gg1, gg2, gradient_gg1, gradient_gg2, x0; max_iter=100, tol=1e-6)
	    x = x0
	    H = Matrix{Float64}(I, length(x0), length(x0))  # Initial Hessian approximation
	    x_history = [x0]  # Initialize history of x values
	    for iter in 1:max_iter
	        grad_ff = gradient_ff(x)
	        c1 = gg1(x)
	        c2 = gg2(x)
	        A1 = gradient_gg1(x)
	        A2 = gradient_gg2(x)
	        
	        # Check constraints
	        if c1 > 0 || c2 > 0
	            println("Constraints violated at iteration $iter")
	            break
	        end
	        
	        # Solve the quadratic subproblem with constraints
	        p = -H \ grad_ff  # Search direction (simplified, should include constraints)
	        alpha = 1.0  # Step length (can be optimized further)
	        x_new = x + alpha * p
	        s = x_new - x
	        y = gradient_ff(x_new) - grad_ff
	        H = bfgs_update(H, s, y)
	        x = x_new
	        push!(x_history, x)  # Store the new x value
	        if norm(grad_ff) < tol
	            break
	        end
	    end
	    return x, x_history
	end
	
	# Initial guess
	x0 = [0.5, 0.5]
	# Run SQP
	optimal_x, x_history = sqp(ff, gradient_ff, gg1, gg2, gradient_gg1, gradient_gg2, x0)
	println("Optimal solution: $optimal_x")
	println("History of x values: $x_history")
end

# ╔═╡ 47d91953-78a2-46a2-8e38-d28c12fe33bf
I

# ╔═╡ Cell order:
# ╠═7c7d889c-a3df-11ef-1699-f99f594e0205
# ╠═90ddcb7c-d784-4ac5-82bb-02bb71196fc8
# ╠═b171beff-1510-4320-baea-cc2e92aaa229
# ╠═35be0489-ca00-4f5c-a604-7697f19740d4
# ╟─5f4eabf2-6949-42ee-b458-1e76742f4b47
# ╠═84612895-9a53-405d-b1c7-f965872bce0a
# ╠═7681215a-728f-4001-9af3-8c485694fc8e
# ╠═9162ffab-1f53-4ace-8dfb-d86c861e0503
# ╠═93604d82-3a0e-410b-b161-86e5d1dbe9ee
# ╠═8d6abc20-b1c1-4345-ba46-381d21dce36a
# ╠═04f88027-fa90-48df-9b56-452a0a05a627
# ╠═332e4e2a-8938-49ac-9e89-9bee128a34b6
# ╠═0d246343-119f-4e93-ad79-66a5c60f73cb
# ╠═30a1b1d5-5e91-4c38-87a2-109b46745d47
# ╠═c84bc9ce-77f9-41e1-aed0-3a65becb0e25
# ╠═47d91953-78a2-46a2-8e38-d28c12fe33bf
