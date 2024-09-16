### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e8cec1d6-a477-430e-b6dc-e29bd84b4244
using Pkg; Pkg.activate()

# ╔═╡ 822befd8-36b4-11ed-3d98-7587449763ff
begin
	using PlutoUI, LaTeXStrings, Colors, ColorSchemes, StaticArrays
	using ForwardDiff, Symbolics, LinearAlgebra,Optim
	using CairoMakie
	set_theme!(theme_latexfonts(), fontsize=18)
	#lift_parent_attribute(scene, :font, "Computer Modern")
end

# ╔═╡ 33fc4bcd-044b-4ac6-a85a-a333ce6eaea1
# from  algorithms for Optimization book pp 36 
# Algorithm 3.1
function bracket_minimum(f, x=0; s=1e-2, k=2.0)
	a, ya = x, f(x)
	b, yb = a + s, f(a + s)
	
	if yb > ya 
		a, b = b, a 		# swap
		ya, yb = yb , ya
		s = -s
	end
	while true
		c, yc = b + s, f(b + s)
		if yc > yb
			return a < c ? (a, c) : (c, a)
		end
		a, ya, b, yb = b, yb, c, yc
		s *= k
	end
end

# ╔═╡ b61ee663-2ce8-4e8b-8a83-e4ea8c6c1525
# Exact line search
function line_search(f, x, d)
	objective = α -> f(x + α*d)
	a, b = bracket_minimum(objective)
	
	# using brent method from Optim.jl
	res = Optim.optimize(objective, a, b)     # minimize
	α = Optim.minimizer(res)
	return α
end

# ╔═╡ 6b5c04dd-d69e-4da9-8ad1-da96734812de
begin
	abstract type DescentMethod end
	
	mutable struct ConjugateGradientDescent <: DescentMethod
		d
		g
	end

	function init!(M::ConjugateGradientDescent, f, ∇f, x)
		M.g = ∇f(x)
		M.d = -M.g

		return M
	end
	function step!(M::ConjugateGradientDescent, f, ∇f, x)
		d, g = M.d, M.g
		g′ = ∇f(x)		#g\prime <TAB>
		β = max(0, dot(g′, g′ - g)/dot(g,g)) 	# Polok-Ribiere
		#β = dot(g′, g′)/dot(g,g)
		d′ = -g′/norm(g′) + β*d 
		
		α = line_search(f, x, d′)
		x′ = x + α*d′
		M.d, M.g = d′, g′

		return x′
	end
		
end

# ╔═╡ 10b834f0-e284-489c-8e25-17bd4c1c7708
function g_step(α, f, pk, x)
	α = line_search(f, x, pk)
	
	return α
end

# ╔═╡ cdba9b91-bf96-4e55-9de3-732f65ad8a6f
md"""
β = $(@bind β PlutoUI.Slider(1:1:20; show_value=true, default=3))
"""

# ╔═╡ b88cc623-3112-4210-9df5-b0325e6a4a1d
begin
	f(x) = x[1]^2 + β*x[2]^2 
	∇f(x) = ForwardDiff.gradient(f,x)
	pf(x,y) = f([x,y])
	
	# bean function
	fb(x; a=1, b=0.5) = (a-x[1])^2 + (1-x[2])^2 + (b)*(2x[2] - x[1]^2)^2
	#fb(x) = f(x)
	∇fb(x) = ForwardDiff.gradient(fb,x)		
	pfb(x,y) = fb([x, y])
end

# ╔═╡ 2fd32917-d47a-4ff0-aed1-ca34a8ca2586
begin
	@variables x₁ x₂
	
	ft = latexstring("f(x_1,x_2) = $(pf(x₁,x₂))")
end

# ╔═╡ 8f0ffe59-f808-4f8f-bb32-0f022a2047f5
let
	# Quadratic with degree 2, the CG can converse in two-step.
	ε_G = 1e-10; ε_A = 1e-10; ε_R = 1e-10;
	N = 500
	xx = zeros(2,N)
	# Initial Point
	x0 = [10, 1]
	#x0 = [1, 1]
	d0 = -∇f(x0)/norm(∇f(x0))
	α0 = line_search(f, x0, d0)

	# Step 1
	x1 = x0 + α0*d0
	β0 = norm(∇f(x1))^2/norm(∇f(x0))^2
	xx[:,1:2] = [x0; x1]
	
	if norm(∇f(x1)) > ε_G
		# Step 2
		d1 = -∇f(x1)/norm(∇f(x0)) + β0*d0
		α1 = line_search(f, x1, d1)
		x2 = x1 + α1*d1
	
		
		xx[:,1:3] = [x0; x1; x2]  

		global xs2 = xx[:,1:3]
		Mcon = ConjugateGradientDescent(d1,∇f(x2)/norm(∇f(x2)))
		init!(Mcon, f, ∇f, xx[:,3])
	
	
	
		for i = 4:N							
			# conjugate gradient
			g = ∇f(xx[:,i-1])

			# first criteria
			if norm(g) <= ε_G
				break
			end
			xx[:,i] = step!(Mcon, f, ∇f, xx[:, i-1])
		
		
			nX = abs(f(xx[:,i]) - f(xx[:,i-1]))
			if nX < ε_A + ε_R*abs(f(xx[:,i-1]))
				global xs2 = xx[:,1:i-1]
				break
			else
				global xs2 = xx[:,1:i-1]
			end
		
		end
	else
		global xs2 = xx[:,1:2]
	end

	#---------------------------
	# Plot
	xxs1 = -5:0.05:12;
	xxs2 = -5:0.05:12;

	fig0 = Figure(size = (800,400), font = "CMU Serif");
	ax = Axis(fig0[1,1], xlabel = L"x_1", ylabel = L"x_2",  aspect = AxisAspect(1.3))
	limits!(ax, -5,12,-5, 12)

	text!(ax, x0[1], x0[2]-0.75, text= L"x_0") 
	contour!(ax,xxs1, xxs2, pf, levels=-100:20:500, color = (:blue,0.9), 
		linewidth=2) #
	contour!(ax,xxs1, xxs2, pf, levels=[2,0.3], color= (:blue, 0.9), linewidth=2)
	
	Np = size(xs2,2)
	scatterlines!(ax, xs2[1,1:Np], xs2[2,1:Np], linewidth=2, color=:red,
		markersize=15, strokecolor=:black, strokewidth=2, 
		label=("Conjugate gradient  with $(size(xs2,2)-1) iterations"))
	
	scatter!(ax,xs2[1,end], xs2[2,end], color=:red, markersize=15, 
		strokecolor=:black, strokewidth=2)
	
	text!(ax, xs2[1,end] + 0.2 , xs2[2,end]-1, text= L"x_\ast") 
	
	axislegend(ax; labelsize=18)
	
	fig0
end

# ╔═╡ cb86ca7c-0fbe-40fd-bede-9f6c00f97d80
md"""
Using **$(size(xs2,2)-1)** iteration\
Optimal solution **$x$= [ $(round(xs2[1,size(xs2,2)], digits=4)), $(round(xs2[2,size(xs2,2)], digits=4))]**\
Optimal value **$f(x)$= $(round(f(xs2[:,size(xs2,2)]), digits=4))**
"""

# ╔═╡ 84798745-2e93-406b-8575-48db6e1b26a5
#save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/conjugate_1.pdf", fig0)

# ╔═╡ 08f03238-e0c0-4f41-bb18-862cf6d792f1
begin
	md"""
	x = [ $(@bind x01 PlutoUI.Slider(-2:0.1:2, show_value=true, default = -1.2)), $(@bind x02 PlutoUI.Slider(-1:0.1:3, show_value=true, default=2))]
	"""
end

# ╔═╡ 56f7dd10-d8eb-4002-a184-84bce9560437
let
	
	# Initial Point
	x0 = [x01, x02]
	d0 = -∇fb(x0)/norm(∇fb(x0))
	α0 = line_search(fb, x0, d0)

	# Step 1
	x1 = x0 + α0*d0
	β0 = norm(∇fb(x1))^2/norm(∇fb(x0))^2

	# Step 2
	d1 = -∇fb(x1)/norm(∇fb(x0)) + β0*d0
	α1 = line_search(fb, x1, d1)
	x2 = x1 + α1*d1

	
	N = 1000
	xx = zeros(2,N)
	xx[:,1:3] = [x0; x1; x2]

	global xs3 = xx[:,1:3]
	Mcon = ConjugateGradientDescent(∇fb(x0)/norm(∇fb(x0)),∇fb(x0))
	init!(Mcon, fb, ∇fb, xx[:,3])
	ε_G = 1e-10; ε_A = 1e-10; ε_R = 1e-10;
	 
	for i = 4:N	
		# conjugate gradient
		g = ∇fb(xx[:,i-1])

		# first criteria
		if norm(g) <= ε_G
			break
		end
		# conjugate gradient
		xx[:,i] = step!(Mcon, fb, ∇fb, xx[:, i-1])
		
		nX = abs(fb(xx[:,i]) - fb(xx[:,i-1]))
		if nX < ε_A + ε_R*abs(fb(xx[:,i-1]))
			global xs3 = xx[:,1:i-1]
			break
		else
			global xs3 = xx[:,1:i-1]
		end
	
	end
	xxs1 = -2:0.05:3;
	xxs2 = -1:0.05:3;
	
	fig1 = Figure(size = (800,400));

	bx1 = CairoMakie.Axis(fig1[1,1], xlabel = L"x_1", ylabel = L"x_2",  aspect = AxisAspect(1.3))
	limits!(bx1, -2,3,-1, 3)

	text!(bx1, x1[1]+0.1, x1[2], text= L"x_0") 
	contour!(bx1,xxs1, xxs2, pfb, levels=-0:5:200, space=:data,
				color = (:blue, 0.9), linewidth=2) #
	contour!(bx1,xxs1, xxs2, pfb, levels=[2,0.3, 0.11], color=(:blue, 0.9), 
		linewidth=2)
	
	Np = size(xs3,2)
	#scatterlines!(bx1,xs3[1,1:Np], xs3[2,1:Np], color=:red, markercolor=:blue, 
	#	markersize=10, linewidth=2, label=("CG with $(size(xs3,2)) iterations"))
	scatterlines!(bx1,xs3[1,1:Np], xs3[2,1:Np], color=:red, markercolor=:blue, 
		strokecolor=:red, strokewidth=1, linewidth=2, markersize=10, label=("SD with $(size(xs3,2)) iterations"))
	
	scatter!(bx1,xs3[1,end], xs3[2,end], color=:red, markersize=7)
	
	text!(bx1, xs3[1,end] + 0.2 , xs3[2,end]-0.1, text= L"x_\ast") 
	
	axislegend(bx1; labelsize=18)
	
	fig1
	
end

# ╔═╡ f6066c19-0b06-4d3b-9c9f-97ef25fc9103
let
	
	# Initial Point
	x1 = [x01, x02]
	α0 = 0.1
	
	N = 500
	xx = zeros(2,N)
	xx[:,1] = x1
	
	ε_G = 1e-10; ε_A = 1e-10; ε_R = 1e-10;   
	for i = 2:N			
		# α is defined using line_search 
		g = ∇fb(xx[:,i-1])

		# first criteria
		if norm(g) <= ε_G
			break
		end

		# line_search
		pk = -g/norm(g)
		α = g_step(α0, fb, pk, xx[:, i-1])
		# the learning rate doesn't change
		#global α0 = α
		
		xx[:,i] = xx[:,i-1] + α*pk	

		nX = abs(fb(xx[:,i]) - fb(xx[:,i-1]))
		if nX < ε_A + ε_R*abs(fb(xx[:,i-1]))
			global xs1 = xx[:,1:i-1]
			break
		end
		xs1 = xx[:,1:i]
	end
	
	#-------------------------------------------------------------------------------
	# Plot
	#-------------------------------------------------------------------------------
	xxs1 = -2:0.05:3;
	xxs2 = -1:0.05:3;
	
	fig1 = Figure(size = (800,400));
	
	bx = CairoMakie.Axis(fig1[1,1], xlabel = L"x_1", ylabel = L"x_2",  
		aspect = AxisAspect(1.3))
	limits!(bx,-2,3,-1,3)
	
	text!(bx, x1[1], x1[2] + 0.05, text= L"x_0") 
	contour!(bx,xxs1, xxs2, pfb, levels=-0:5:200, linewidth=2,
		color=(:blue, 0.8)) #
	contour!(bx,xxs1, xxs2, pfb, levels=[2,0.3, 0.11], color=(:blue, 0.9), 
		linewidth=2)
	Np = size(xs1,2)
	scatterlines!(bx,xs1[1,1:Np], xs1[2,1:Np], color=:red, markercolor=:blue, 
		strokecolor=:red, strokewidth=1, linewidth=2, markersize=10, label=("SD with $(size(xs1,2)-1) iterations"))
	
	# The optimal value of pf is at x = [1, 1]
	scatter!(bx,xs1[1,end], xs1[2,end], color=:red, markersize=10)
	
	text!(bx, xs1[1,end] +0.1, xs1[2,end]-0.3, text= L"x_\ast") 
	
	axislegend(bx; labelsize=18)
	
	#empty!(fig1)
	fig1
	
end

# ╔═╡ ae877346-3cf2-4db0-aa8c-5c7a86ce2ee2
#save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/compare_CG_SD.pdf", fig1)

# ╔═╡ Cell order:
# ╠═e8cec1d6-a477-430e-b6dc-e29bd84b4244
# ╠═822befd8-36b4-11ed-3d98-7587449763ff
# ╠═33fc4bcd-044b-4ac6-a85a-a333ce6eaea1
# ╠═b61ee663-2ce8-4e8b-8a83-e4ea8c6c1525
# ╠═b88cc623-3112-4210-9df5-b0325e6a4a1d
# ╠═6b5c04dd-d69e-4da9-8ad1-da96734812de
# ╠═10b834f0-e284-489c-8e25-17bd4c1c7708
# ╠═cdba9b91-bf96-4e55-9de3-732f65ad8a6f
# ╟─2fd32917-d47a-4ff0-aed1-ca34a8ca2586
# ╟─cb86ca7c-0fbe-40fd-bede-9f6c00f97d80
# ╠═8f0ffe59-f808-4f8f-bb32-0f022a2047f5
# ╠═84798745-2e93-406b-8575-48db6e1b26a5
# ╠═56f7dd10-d8eb-4002-a184-84bce9560437
# ╟─08f03238-e0c0-4f41-bb18-862cf6d792f1
# ╠═f6066c19-0b06-4d3b-9c9f-97ef25fc9103
# ╠═ae877346-3cf2-4db0-aa8c-5c7a86ce2ee2
