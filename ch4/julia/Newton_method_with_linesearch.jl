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

# ╔═╡ c745146e-823a-11ee-1c63-a995bb36c382
begin
	using Pkg; Pkg.activate()
end

# ╔═╡ 9d9cc882-0077-48cf-80e1-e9e015127ff3
begin
	using LinearAlgebra
	#using Plots; default(fontfamily="Computer Modern", guidefontsize=12,
	#	tickfontsize=9, framestyle=:box) #LaTeX Style
	using CairoMakie
		set_theme!(theme_latexfonts(),fontsize=18)
	using ForwardDiff
	using PlutoUI
	using Symbolics, LaTeXStrings
	using JuMP, Ipopt
end

# ╔═╡ bee0071a-e636-4944-918f-d80d4958d75e
md"""
### Bean Function
"""

# ╔═╡ 8f3a8415-618c-48ef-89aa-8f4f028aba2a
begin
	f1 = x -> (1-x[1])^2 + (1-x[2])^2 + 0.5*(2x[2] - x[1]^2)^2
	#f1 = x -> x[1] + 2x[2]
	#f1 = x -> x[1]^2 + x[2]^2 
	#c = x->  x[1]^2 + 2x[1] - x[2]
	
	∂c = x->[2*x[1]+2 -1]
	pf1(x,y) = f1([x, y])
end

# ╔═╡ a542fc3c-1951-490d-97ed-627fe4590738
begin
	# For scalar input x
	function ∇f(f,x)
		return ForwardDiff.gradient(f,[x])
	end

	# For a vector input x
	function ∇f(f, x::Array)
		return ForwardDiff.gradient(f,x)
	end

	function ∇f(f::Array, x::Array)
		return ForwardDiff.jacobian(x) do f end
	end

	# For a scalar input x
	function ∇2f(f,x)
		# to specific that g is a gradient of f
		return ForwardDiff.hessian(f,[x])
	end
	# For a vector input x
	function ∇2f(f, x::Array)
		return ForwardDiff.hessian(f,x)
	end
	
end

# ╔═╡ 64196d3e-b82e-4ae4-819e-05a4dac54a1e
md"""
### No constraints
"""

# ╔═╡ 7fcd8f3d-cd44-447a-8676-be8120b57604
begin
	function newton_step(f, x0)
		return [x0] - ∇2f(f, x0)\ ∇f(f, x0)
	end

	# for multivariable function
	function newton_step(f, x0::Array)
		return x0 - ∇2f(f, x0)\ ∇f(f,x0)
	end
end

# ╔═╡ 9c278d10-31b3-4011-9621-1cd3fa1836c4
begin
	function mod_chol(A)
		τ = 0
		if minimum(diag(A)) > 0
			τ = 0
		else
			τ = 0.5*norm(A, 2)  # Frobenius norm 
		end
		
		A = A + τ*I
		
		while ~isposdef(A)
			τ = max(2*τ, 0.5*norm(A, 2))
			A =  A + τ*I
		end
		C = cholesky(A)
		
		return C.L	
	end
end

# ╔═╡ ab71ecaa-08d1-4a5f-8105-7be741d40cf8
begin
	AA = [2 7; 7 8]
	C = cholesky(AA, check = false)
	issuccess(C)
	#mod_chol(AA)
	max(0, 0.5*norm(AA, 2))
	maximum(diag(AA))
end

# ╔═╡ bde59df7-37e6-426a-9588-16594b729f5d
begin
	function backtracking(f, x, d; β = 1e-4, ρ = 0.5, α = 1)
		
		objective = α -> f(x + α*d)
		while (objective(α) > f(x) + β*α*∇f(f, x)' * d)[1]
			α = ρ * α
		end

		return α
	end

	function backtracking(f, x::Array, d::Array; β = 1e-4, ρ = 0.5, α = 1)
		
		objective = α -> f(x + α*d)				
		while (objective(α) > f(x) + β*α* ∇f(f, x)' * d)
			α = ρ * α
		end

		return α
	end
	
end

# ╔═╡ ed9cea50-c521-461d-96a7-d438286e4dc5
begin
	function newton_step_with_ls(f, x0)
		g = ∇f(f, x0); H = ∇2f(f, x0);
		L = mod_chol(H)
		zk = L\g;
		dk = L'\(-zk)
	    α = backtracking(f, x0, dk)
		
		return [x0] + α*dk
	end

	# for multivariable function
	function newton_step_with_ls(f, x0::Array)
		g = ∇f(f, x0); H = ∇2f(f, x0);
		L = mod_chol(H)
		zk = L \ g;
		dk = L'\ (-zk)
	    α = backtracking(f, x0, dk)	
		
		return x0 + α*dk
	end
end

# ╔═╡ d5cb99e2-de14-400e-bfe4-0696441668ea
md"""
 $n$-Step = $(@bind ns  PlutoUI.Slider(2:100, default=100, show_value=true))

 $x_1$ = $(@bind xv01 PlutoUI.Slider(-2.5:0.1:2.5, default= -1, show_value=true)) 
 
 $x_2$ = $(@bind xv02 PlutoUI.Slider(-2.5:0.1:2.5, default= 2, show_value=true))
"""

# ╔═╡ a4d4389d-a0c2-414a-8bb0-b5f50699d4fb
begin
	xv0 = [xv01; xv02]
	εₙ = 1e-11
end

# ╔═╡ 708276ef-3cfb-4ded-9689-16b9707be47f
∇f(f1,xv0)'*[1, 2]

# ╔═╡ 22b5e8b4-e073-464e-8738-ab562965da75
begin
	xs1 = []
	xguess1 = Float64[]; xguess1 = append!(xguess1, xv0)

	for i = 2:ns
		global xnew1 = newton_step(f1, xguess1[:,end])
		if norm(∇f(f1, xguess1[:,end])) <= εₙ
	       global xs1 = xguess1[:,1:i-1]
			break
		end
		global xguess1 = [xguess1 xnew1]
	end
	xs1 = xguess1[:,:]
end

# ╔═╡ 42f7341b-3376-4dba-95ee-6848d5d97010
begin
	xs2 = Float64[]
	xguess2 = Float64[]; xguess2 = append!(xguess2, xv0)

	for i = 2:ns
		global xnew2 = newton_step_with_ls(f1, xguess2[:,end])
		g = ∇f(f1, xguess2[:,end])
		H = ∇2f(f1, xguess2[:,end])
		if norm(g'*(H\g)) <= εₙ
			global xs2 = xguess2[:,1:1]
			break
		end
		global xguess2 = [xguess2 xnew2]
	end
	xs2 = xguess2[:,:]
end

# ╔═╡ 832a046b-e241-4cad-89de-04690837c4f7
begin
	xxs1 = -2.6:0.05:2.6;
	xxs2 = -2.6:0.05:4;
	
	fig1 = Figure(size = (600,400), font = "CMU Serif")
	ax1 = CairoMakie.Axis(fig1[1,1], xlabel = L"x_1", ylabel = L"x_2",  
		aspect = AxisAspect(1.2),xlabelsize=20, ylabelsize=20)
	limits!(ax1, -2.6,2.6,-2.6, 4)
   
	text!(ax1, 1, 1.0, text= L"x^\ast", fontsize=25, color=:black) 
	text!(ax1, xv0[1,1], xv0[2,1]+0.2, text= L"x_0", fontsize=25, color=:black) 
	
	level1 = [0, 0.1, 0.15, 0.5,1, 2, 3, 4, 5, 6, 7, 10, 15, 30]
	contour!(ax1,xxs1, xxs2, pf1, levels= -100:40:500,
				 color=:blue, linewidth=1) 
	contour!(ax1,xxs1, xxs2, pf1, levels=level1,
				color=:blue, linewidth=1)
	
	
	scatterlines!(ax1,xguess1[1,:], xguess1[2,:], color=:red, linestyle=:solid, 
		linewidth = 2, marker=:circle, markercolor=:cyan, markersize=15, label=("Newton with $(size(xguess1,2)) iterations"), strokecolor=:blue, strokewidth=2)
	
	scatterlines!(ax1,xs2[1,:], xs2[2,:], color=:green, linestyle=:dash, linewidth=2,
		marker=:circle, markercolor=:lightgreen, markersize=15, label=("Newton LS with $(size(xs2,2)) iterations"), strokecolor=:green, strokewidth=2)
	
	axislegend(ax1; labelsize=14, position=:rb)
	
	#empty!(fig5)
	fig1
end

# ╔═╡ a27a78db-c3e9-48b5-accc-65fcb73c70eb
# ╠═╡ disabled = true
#=╠═╡
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("newton_ls2b.pdf", fig1)
	run(`pdfcrop  --clip newton_ls2b.pdf newton_ls2.pdf`)
	run(`rm newton_ls2b.pdf`)
end
  ╠═╡ =#

# ╔═╡ 940d736c-cfd2-4d11-b119-b6c6021a9405
md"""
### With Linear Constraints
"""

# ╔═╡ ecbb6913-d499-432b-88c2-051fbb5f9c7d
begin
	function newton_step_lc(f, A, b, c, x0, λ0)
		J = ForwardDiff.jacobian(x -> ∇f(c, x)'*λ0[end],x0)
		H = ∇2f(f, x0)  + J
		g = ∇f(f, x0)
		  
		Δz = [H A'; A zeros(size(A,1), size(A',2))] \ [-g - A'*λ0; -c(x0)]# zeros(size(A,1),1)]
		
		
		λk1 = Δz[3] #-(A*(H\A'))*A*(H\g)
		dk = Δz[1:2] #-H \ (A'*λk1 + g) #Δz[1:size(A,1)]
		
		α =  backtracking(f, x0, dk)	
		#println(α, λk1, x0+α*dk)

		return x0+ α*dk, λk1
	end
end

# ╔═╡ be1ffaf9-7cdf-4b88-9ab1-5c7645ce7431
begin
	function newton_step_lc_inf(f, A, b, c, x0, λ0)
		#J = ForwardDiff.jacobian(x -> ∇f(c, x)'*λ0[end],x0)
		H = ∇2f(f, x0) # + J
		g = ∇f(f, x0)
		
		rd = g + A'*λ0                 # dual residual
		#println(A*x0 -b)
		rp = A*x0 - b                  # primal residual
		
		  
		Δz = [H A'; A zeros(size(A,1), size(A',2))] \ [-rd; -rp]# zeros(size(A,1),1)]
		
		
		δk = Δz[3] #-(A*(H\A'))*A*(H\g)
		dk = Δz[1:2] #-H \ (A'*λk1 + g) #Δz[1:size(A,1)]

		 
		α = 1; ρ = 0.1; γ = 0.5
		 
		#while norm([g + A'*(λ0 + α*δk); A*(x0 + α*dk) - b]) > (1 - ρ*α)*norm([g + #A'*λ0; A*x0 - b])
		#	α = γ * α
		#end
		#println(norm([g + A'*(λ0 + α*δk); A*(x0 + α*dk) - b]))
		
		return x0+ α*dk, λ0 + α*δk
	end
end

# ╔═╡ 741d02ca-af70-4e5c-88cc-67f6ac2e806f
function newton_step2(f, c, x0, λ0)
	# Should be ∇f(c,x) because we need to substitute the value of x after the jacobian
	J = ForwardDiff.jacobian(x -> ∇f(c, x)'*λ0[end],x0)
    H = ∇2f(f, x0) + J
	#println(H, J) for check bug only
    C = ∇f(c, x0)'
		
    Δz = [H C'; C zeros(size(C,1), size(C',2))]\[-∇f(f,x0)-C'*λ0; -c(x0)]
    Δx = Δz[1:2]
    Δλ = Δz[3]
    return x0+Δx, λ0+Δλ
end

# ╔═╡ 91c4d639-c27a-4db2-9756-7aa3079b724c
begin
	c1 = x -> 2x[1] + 3x[2] +1 
	c2 = x -> 6*x[1] + 8x[2]
	c3 = x -> (1/4)*x[1]^2 + x[2]^2 - 1 
	c = c2
	cp(x,y) = c([x,y])
	#cp2(x,y) = c2([x,y])

	# when change c don't forget to change 
	bl = 1
	
	Al = ∇f(c,[1, 2])';
	#xc = [cp1(1,1); cp2(1,2)]
end

# ╔═╡ 7cf23816-ff86-412e-acb0-d8aaa40ebd9b
md"""
 $n$-Step = $(@bind nss  PlutoUI.Slider(2:2000, default=2000, show_value=true))

 $x_1$ = $(@bind xvc01 PlutoUI.Slider(-2.5:0.1:2.5, default= -1, show_value=true)) 
 
 $x_2$ = $(@bind xvc02 PlutoUI.Slider(-2.5:0.1:2.5, default= 2, show_value=true))
"""

# ╔═╡ d5715cbd-edf9-4789-bbad-2c31af535e91
begin
	#xv0 = [-1; -1]; xv0 = [-3, 2]
	xvc0 = [xvc01; xvc02]
	λv0 = [1];  xs4=[]

	xguess4 =  Float64[]; xguess4 = append!(xguess4, xvc0)
	λguess = Float64[]; λguess = append!(λguess, λv0)
	
	xnew, λnew = newton_step2(f1, c,xguess4[:,end], λguess[end])
	 
	xguess4 = [xguess4 xnew]
	λguess = [λguess λnew]
	
	for i=2:nss
		global xnew, λnew = newton_step2(f1, c, xguess4[:,end],λguess[end])
		
		if norm(xguess4[:,i]-xguess4[:,i-1]) <= εₙ
			global xs4 = xguess4[:,1:i-1]
			break
		end
		
		global xguess4 = [xguess4 xnew]
		global λguess = [λguess λnew]
	end
	
	xs4 = xguess4[:,:]
end

# ╔═╡ e26e6c7a-5614-4394-bf22-99387782926a
begin
	#xv0 = [-1; -1]; xv0 = [-3, 2]
	#xvc0 = [xvc01; xvc02]
	#λv0 = [1];  
	xs5=[]

	xguess5 =  Float64[]; xguess5 = append!(xguess5, xvc0)
	λguess5 = Float64[]; λguess5 = append!(λguess5, λv0)
	
	xnew5, λnew5 = newton_step_lc_inf(f1, Al, bl, c, xguess5[:,end], λguess5[end])
	 
	xguess5 = [xguess5 xnew5]
	λguess5 = [λguess5 λnew5]
	
	for i=2:nss
		global xnew5, λnew5 = newton_step_lc_inf(f1, Al, bl, c ,
			xguess5[:,end],λguess5[end])
		
		if norm(xguess5[:,i]-xguess5[:,i-1]) <= εₙ
#if norm([∇f(f1,xguess5[:,end]) - Al'*λguess5[end]; Al*xguess5[:,end] - bl]) <= εₙ
			
			global xs5 = xguess5[:,1:i-1]
			break
		end
		
		global xguess5 = [xguess5 xnew5]
		global λguess5 = [λguess5 λnew5]
	end
	
	xs5 = xguess5[:,:]
end

# ╔═╡ 5e5961a3-f90d-46b9-b66c-5432ab1cd4c4
begin
	#xs1 = -5:0.05:5;
	#xs2 = -5:0.05:5;
	
	fig2 = Figure(size = (600,400), font = "CMU Serif")
	ax2 = Axis(fig2[1,1], xlabel = L"x_1", ylabel = L"x_2",  
		aspect = AxisAspect(1.2),xlabelsize=20, ylabelsize=20)
	limits!(ax2, -2.6,2.6,-2.6, 4)
   
	text!(ax2, xs4[1,end], xs4[2,end]+ 0.1, text= L"x^\ast", 
		fontsize=25, color=:black) 
	text!(ax2, xvc0[1,1]-0.3, xvc0[2,1]+0.1, text= L"x_0", fontsize=25, color=:black) 
	#level1 = [0, 0.1, 0.15, 0.5,1, 2, 3, 4, 5, 6, 7, 10]
	contour!(ax2,xxs1, xxs2, pf1,  levels= -100:40:500
				, color=:blue, linewidth=1) #
	contour!(ax2,xxs1, xxs2, pf1,  levels=level1,
				color=:blue, linewidth=1)
#	contour!(ax2,xxs1, xxs2, cp1, levels=0:0, linewidth=3)
	contour!(ax2,xxs1, xxs2, cp, levels=0:0, linewidth=3, color=:green)
	
	scatterlines!(ax2,xs4[1,:], xs4[2,:], color=:red, linestyle=:dash, linewidth=2,
		marker=:circle, markercolor=:orange, markersize=15, label=("NL Con with $(size(xs4,2)) iterations"), strokecolor=:red, strokewidth=2)
	scatterlines!(ax2,xs5[1,:], xs5[2,:], color=:black, linestyle=:dash, linewidth=2,
		marker=:circle, markercolor=:gray, markersize=15, label=("Linear Con with $(size(xs5,2)) iterations"), strokecolor=:black, strokewidth=2)
	
	axislegend(ax2; labelsize=14, position=:lt)
	
	#empty!(fig5)
	fig2
end

# ╔═╡ 3b8cdd93-d417-42d0-8dea-f33e3544ce9c
# ╠═╡ disabled = true
#=╠═╡
begin
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
		
	save("newton_ls_NL1b.pdf", fig2)
	run(`pdfcrop  --clip newton_ls_NL1b.pdf newton_ls_NL1.pdf`)
	run(`rm newton_ls_NL1b.pdf`)
end
  ╠═╡ =#

# ╔═╡ 0f902c4e-7ad1-4032-884b-10a140233955
function solve_constrained_lp(f, A::Matrix, b::Vector, xv0)
           m, n = size(A)
		
           model = Model(Ipopt.Optimizer)
           #set_silent(model)
		    
           @variable(model, x[1:n])
	
           #@variable(model, residuals[1:m])
           @constraint(model, con, A * x == b)

		   set_start_value(con, xv0)	
           @objective(model, Min, f(x))

		  #set_optimizer_attribute(model, "print_level", 7)	
           optimize!(model)
           return value.(x)
       end

# ╔═╡ 33b53170-932d-47d8-b04e-9b18d8ffc676
begin
  x =  solve_constrained_lp(f1, Matrix(Al), [bl], xv0)	
	x
end

# ╔═╡ Cell order:
# ╠═c745146e-823a-11ee-1c63-a995bb36c382
# ╠═9d9cc882-0077-48cf-80e1-e9e015127ff3
# ╟─bee0071a-e636-4944-918f-d80d4958d75e
# ╠═8f3a8415-618c-48ef-89aa-8f4f028aba2a
# ╠═a542fc3c-1951-490d-97ed-627fe4590738
# ╟─64196d3e-b82e-4ae4-819e-05a4dac54a1e
# ╠═7fcd8f3d-cd44-447a-8676-be8120b57604
# ╠═9c278d10-31b3-4011-9621-1cd3fa1836c4
# ╠═ab71ecaa-08d1-4a5f-8105-7be741d40cf8
# ╠═bde59df7-37e6-426a-9588-16594b729f5d
# ╠═708276ef-3cfb-4ded-9689-16b9707be47f
# ╠═ed9cea50-c521-461d-96a7-d438286e4dc5
# ╠═a4d4389d-a0c2-414a-8bb0-b5f50699d4fb
# ╠═22b5e8b4-e073-464e-8738-ab562965da75
# ╠═42f7341b-3376-4dba-95ee-6848d5d97010
# ╟─d5cb99e2-de14-400e-bfe4-0696441668ea
# ╠═832a046b-e241-4cad-89de-04690837c4f7
# ╠═a27a78db-c3e9-48b5-accc-65fcb73c70eb
# ╠═940d736c-cfd2-4d11-b119-b6c6021a9405
# ╠═ecbb6913-d499-432b-88c2-051fbb5f9c7d
# ╠═be1ffaf9-7cdf-4b88-9ab1-5c7645ce7431
# ╠═741d02ca-af70-4e5c-88cc-67f6ac2e806f
# ╠═91c4d639-c27a-4db2-9756-7aa3079b724c
# ╠═d5715cbd-edf9-4789-bbad-2c31af535e91
# ╠═e26e6c7a-5614-4394-bf22-99387782926a
# ╟─7cf23816-ff86-412e-acb0-d8aaa40ebd9b
# ╠═5e5961a3-f90d-46b9-b66c-5432ab1cd4c4
# ╠═3b8cdd93-d417-42d0-8dea-f33e3544ce9c
# ╠═0f902c4e-7ad1-4032-884b-10a140233955
# ╠═33b53170-932d-47d8-b04e-9b18d8ffc676
