### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 42c440f8-7661-11f0-279f-dd279171e082
using Pkg; Pkg.activate()

# ╔═╡ 18c17031-0ca4-4a01-9c2f-a43f88465f31
using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)

# ╔═╡ 3af2e4bb-e82f-446a-9443-ca3a930f1668
using JuMP, Ipopt   # for NLP

# ╔═╡ f0d1fb80-56da-4ba0-b31b-efc257bdacf9
begin
	using NLsolve
	
	# cross section
	function f1!(F, x)
		F[1] = x[2] - 2x[1]
		F[2] = x[1]^2 - x[2]
	end

	function f2!(F, x)
		F[1] = x[1]^2 - x[2]
		F[2] = x[1] + x[2] - 2
	end

	function f3!(F, x)
		F[1] = x[2] - 2x[1]
		F[2] = x[1] + x[2] - 2
	end
	
	sol1 = nlsolve(f1!, [-1.0, 1.0])
	sol2 = nlsolve(f2!, [0.0, 1.0])
	sol3 = nlsolve(f3!, [1.0, 1.0])

	x11, y11 = sol1.zero
	x21, y21 = sol2.zero
	x31, y31 = sol3.zero
	
	println("g1, g2 = $(sol1.zero)")
	println("g2, g3 = $(sol2.zero)")
	println("g1, g3 = $(sol3.zero)")
	println(x31)
end	

# ╔═╡ 71a3472d-b602-4ec0-a3c3-8dd11e837179
md"""
# Graphical Method Example 
"""

# ╔═╡ 252a93c8-1688-49eb-a264-15ebd3d52a4d
md"""
## Example 1
"""

# ╔═╡ 4f4b01db-4cd6-408c-894e-0e18592dde96
md"""
$\begin{align*}
&\operatorname*{minimize}_{x_1, x_2} \qquad (x_1-2)^2 - (x_2-1)^2\\
&\operatorname*{subject to} \qquad x_2 - 2x_1 = 0\\
&\hspace{3cm} x^2 - x_2  \leq 0 \\ 
&\hspace{3cm} x_1 + x_2 \leq 2
\end{align*}$ 
"""

# ╔═╡ 5d901014-2444-4c3f-8fde-fe9063736763
begin
	function exfig4()
	    model = Model(Ipopt.Optimizer)
	    set_silent(model)
	    @variable(model, x1)
	    @variable(model, x2)
	
	    @objective(model, Min, (x1-2)^2 + (x2-1)^2)
		
	    @constraint(model, x2 - 2x1 == 0)
	    @constraint(model, x1^2 - x2 <= 0)
		@constraint(model, x1 + x2 <= 2)
		
	    optimize!(model)
	    print(model)
	    println("Objective value: ", objective_value(model))
	    println("x_1 = ", value(x1))
	    println("x_2 = ", value(x2))

	    return
	end
	
	exfig4()
end

# ╔═╡ 18ac92e0-d4e7-48ac-b14d-c92411d4ac07
let
	global fig2 = Figure(size=(600, 400), fontsize = 18)
	ax = CairoMakie.Axis(fig2[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", aspect=1.2)
	limits!(ax, -1, 4, -1, 3)

	f(x1, x2) = (x1-2)^2 + (x2-1)^2
	g1(x1, x2) = x2 - 2x1
	g2(x1, x2) = x1^2 - x2 
	g3(x1, x2) = x1 + x2 - 2

	
	xs = -1:0.05:4; ys = -1:0.05:4
	la = -1:0.01:-2
	lb = -2:0.2:5
	lc = [lb; la]
	
	contour!(ax, xs, ys, f, levels=lc, label=L"$$Objective function contours")
	contour!(ax, xs, ys, g1, levels=[0], linewidth=2, color=:red)
	contour!(ax, xs, ys, g2, levels=[0], linewidth=2, color=:blue)
	contour!(ax, xs, ys, g3, levels=[0], linewidth=2, color=:green)

	g1_up = x -> 2x
	g2_low = x -> x^2
	g3_up = x -> -x + 2
	xs1_bound = x11:0.001:x31
	xs2_bound = x31:0.05:x21
	
	band!(ax, xs1_bound, g1_up.(xs1_bound), g2_low.(xs1_bound), color=(:cyan, 0.7),
		label=L"$$Feasible region")
	band!(ax, xs2_bound, g2_low.(xs2_bound), g3_up.(xs2_bound), color=(:cyan, 0.7))

	scatter!(ax, [2], [1], marker=:circle, markersize=10, color=:red)
	text!(ax, [2], [1], text="Unconstraint", color=:red)
	scatter!(ax, [0.8], [1.6], marker=:circle, color=:green)
	text!(ax, [0.5], [1.8], text="Only EC", color=:green)
	scatter!(ax, [0.67], [1.33], marker=:circle, color=:blue)
	text!(ax, [-0.8], [1.33], text="EC and IC", color=:blue)
	
	axislegend(ax, labelsize=14, position=:rb)
	
end

# ╔═╡ 7c7ef920-7621-409b-a9fb-0667b4a1f115
fig2

# ╔═╡ 15c76892-daf9-403e-a95c-5ccd664e366f
md"""
## Example 2
"""

# ╔═╡ 1cef2f66-bad6-43e1-871e-1570845560ae
md"""
In this example, we consider the feasible region $\Omega$, and the objective function of the NLP problem.

$\begin{align*}
&\operatorname*{minimize}_{x_1, x_2} \qquad (x_1+2)^2 - x_2\\
&\operatorname*{subject to} \qquad \frac{x_1^2}{4} + x_2 - 1 \leq 0\\
&\hspace{3cm} 2 + x_1 -2x_2 \leq 0
\end{align*}$ 
"""

# ╔═╡ ea9f1f11-8682-4e9a-bc29-c25927a0e909
let
	global fig1 = Figure(size=(600, 400))
	ax = CairoMakie.Axis(fig1[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$")
	ylims!(ax, -2, 2)
	xlims!(ax, -3, 1)

	f(x1, x2) = (x1+2)^2 - x2;
	g1(x1, x2) = x1^2/4 + x2 - 1
	g2(x1, x2) = 2 + x1 - 2x2
	
	la = collect(-1:-0.1:-2.0)
	lb = [-0.2005, -0.5, -0.8]
	lc = [lb; la]

	xs = -3:0.01:1; ys = -3:0.01:3
	contour!(ax, xs, ys, f, levels = lc, label=L"$$Objective function contours")
	contour!(ax, xs, ys, g1, levels = [0], linewidth=2, color=:red)
	contour!(ax, xs, ys, g2, levels = [0], linewidth=2, color=:green)

	g1_up = x -> -x^2/4 + 1
	g2_low = x -> 1 + x/2
	xs_bound = -2:0.1:0
	band!(ax, xs_bound, g1_up.(xs_bound), g2_low.(xs_bound), color=(:cyan, 0.7),
		label=L"$$Feasible region")
	
	axislegend(ax, labelsize=14, position=:rb)
end 

# ╔═╡ 6a737afe-3fef-4d1e-b94f-de79c80e8d85
fig1

# ╔═╡ Cell order:
# ╠═42c440f8-7661-11f0-279f-dd279171e082
# ╠═18c17031-0ca4-4a01-9c2f-a43f88465f31
# ╠═3af2e4bb-e82f-446a-9443-ca3a930f1668
# ╟─71a3472d-b602-4ec0-a3c3-8dd11e837179
# ╟─252a93c8-1688-49eb-a264-15ebd3d52a4d
# ╟─4f4b01db-4cd6-408c-894e-0e18592dde96
# ╠═5d901014-2444-4c3f-8fde-fe9063736763
# ╠═f0d1fb80-56da-4ba0-b31b-efc257bdacf9
# ╠═18ac92e0-d4e7-48ac-b14d-c92411d4ac07
# ╠═7c7ef920-7621-409b-a9fb-0667b4a1f115
# ╟─15c76892-daf9-403e-a95c-5ccd664e366f
# ╟─1cef2f66-bad6-43e1-871e-1570845560ae
# ╠═ea9f1f11-8682-4e9a-bc29-c25927a0e909
# ╠═6a737afe-3fef-4d1e-b94f-de79c80e8d85
