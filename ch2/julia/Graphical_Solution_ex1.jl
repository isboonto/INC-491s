### A Pluto.jl notebook ###
# v0.20.13

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

# ╔═╡ 34cb4d35-319e-48ea-b532-2403d8950301
begin
	import Pkg; Pkg.activate();
end

# ╔═╡ 806be89e-167a-11ed-1f7e-87dd862ce225
begin
	using PlutoUI, LaTeXStrings, StaticArrays
	using ForwardDiff,  LinearAlgebra, Optim
	using CairoMakie
	set_theme!(theme_latexfonts(), fontsize=16)
	
end

# ╔═╡ 2da61ffe-a205-48c7-b7df-3331882d0ec7
using JuMP, Ipopt

# ╔═╡ 41c4cc45-83f9-47ba-a708-19dca0286201
md"""
## Graphical Solution of One- and Two-Variable Problems

For $n = 1$ or $n=2$, It ispossible to plot the NLP problem and obtain the solutions.
"""

# ╔═╡ 28ca0c07-840f-459b-bfb1-ca06445877d6
md"""
### Example 1.10

Plot $f(x) = (x-1)^2$ ,over the domain $x \in [0, 2]$
"""

# ╔═╡ 6c52e6a4-aa96-4c1b-8bed-9afeda5017ad
begin
	f1(x) = (x - 1)^2

	xr1 = 0:0.01:2

	fig1 = Figure(size = (600,400))
	ax = CairoMakie.Axis(fig1[1,1], xlabel = L"$x$", ylabel = L"$f(x)$", aspect = 
		AxisAspect(1.5), xlabelsize=20, ylabelsize=20)
		limits!(ax, 0,2,-0.1, 1.0)

		lines!(ax, xr1, f1, linewidth=2)
	fig1
end

# ╔═╡ 02036da6-4704-48b5-847a-8ed5d45c86ad
md"""
### Example 1.11 
Consider the "corrugated spring" function $f = -\cos(kR) + 0.1R^2$, where $R$ is the radius $R = \sqrt{(x_1 -c_1)^2 + (x_2 - c_2)^2}$, $c_1 = 5$, $c_2=5$, $k=5$

- Draw a 3-D plot of this function
- Draw a contour plot of the function in variable space or $x_1-x_2$ space.
"""

# ╔═╡ fa521ffd-f538-4115-ad94-d1b058315687
md"""
perspectiveness = $(@bind pe PlutoUI.Slider(0:0.1:1, show_value=true, default=0))

elevation = $(@bind el PlutoUI.Slider(1.5:0.5:10, show_value=true, default=8.0))

azimuth =  $(@bind az PlutoUI.Slider(0:0.1:1, show_value=true, default=0.3))
"""

# ╔═╡ 3f7f91f4-1c8a-4276-a403-9620cf086096
begin
	c1 = 5; c2 = 5; k = 5
	f2(x1, x2) = -cos(k*√((x1 - c1)^2 + (x2 - c2)^2) ) + 0.1*((x1 - c1)^2 + (x2 - c2)^2)

	#x21 = 0:0.1:2; x22 = 0:0.1:2
	xs = 8:-0.1:2
	ys = 1:0.1:8 # LinRange(8, -8, 200)
	zs = [f2(x1,x2) for x1 in xs, x2 in ys]
	
	fig2 = Figure(size = (600,800))
	bx = Axis3(fig2[1, 1], xlabel=L"x_1", ylabel=L"x_2", zlabel=L"f(x_1,x_2)", 
		aspect = (1.5, 1.5, 1), perspectiveness = pe, elevation = π / el, azimuth = az*π,)
	
		xw = 8:-0.2:2; yw = 1:0.2:8;
		zw = [f2(x1,x2) for x1 in xw, x2 in yw]
	
		surface!(bx, xs, ys, zs, colormap=:jet, transparency=true, alpha=0.2)
		wireframe!(bx, xw, yw, zw, color=:black, linewidth=0.3)
	
	cx = Axis(fig2[2,1], aspect=1.1, xlabel=L"x_1", ylabel=L"x_2")
		limits!(cx, 2,8, 1, 8)
		contour!(cx, xs, ys, zs)
	
	fig2
end

# ╔═╡ 48e0fe25-9725-4e05-844e-62d7784ca793
md"""
### Example 1.12

In this example, we consider the feasible region $\Omega$ and the objective function of the NLP problem.
"""

# ╔═╡ f63b57b8-a773-4765-b3b7-f6327f9962a2
md"""
$\begin{align*}
	\operatorname*{minimize}_{x_1, x_2} & \qquad (x_1 + 2)^2 - x_2\\
	\text{subject to} & \qquad \frac{x_1^2}{4} + x_2 - 1 \leq 0 \\
                       & \qquad 2 + x_1 - 2x_2 \leq 0 
\end{align*}$
"""

# ╔═╡ ffd53430-b3df-4a31-bd3b-732a548d2f90
begin
	function example_qcp()
	    
		model = Model(Ipopt.Optimizer)
	    set_silent(model)
	    @variable(model, x1)
	    @variable(model, x2)
	
	    @objective(model, Min, (x1+2)^2 - x2)
	    @constraint(model, (x1/2)^2 + x2 - 1 <= 0)
	    @constraint(model, 2 + x1 - 2x2 <= 0)
		
	    optimize!(model)
	    print(model)
	    println("Objective value: ", objective_value(model))
	    println("x_1 = ", value(x1))
	    println("x_2 = ", value(x2))
	    #Test.@test termination_status(model) == LOCALLY_SOLVED
	    #Test.@test primal_status(model) == FEASIBLE_POINT
	    #Test.@test objective_value(model) ≈ 0.32699 atol = 1e-5
	    #Test.@test value(x1) ≈ 0.32699 atol = 1e-5
	    #Test.@test value(x2) ≈ 0.25707 atol = 1e-5
	    
		return value(x1), value(x2)
	end
	x1p, x2p = example_qcp()
	print(x1p)
end

# ╔═╡ 27c448f8-d434-4f77-b0e6-1e9b36c06de8
begin 
	f3(x1, x2) = (x1+2)^2 - x2
	g1(x1, x2) = x1^2/4 + x2 - 1
	g2(x1, x2) = 2 + x1 - 2x2
	
	fig3 = Figure(size = (600,400))
	ax3 = CairoMakie.Axis(fig3[1,1], xlabel = L"x_1", ylabel = L"x_2", 
				aspect = AxisAspect(1.5))
	limits!(ax3, -3.01, 1.01, -3.01, 3.01)

	xs3 = -3:0.01:1; ys3 = -3:0.01:3
	la = collect(-1:-0.1:-2.0);
	lb = [-0.2005, -0.5, -0.8]
	lc = [lb; la]
	

	contour!(ax3, xs3, ys3, f3, levels=lc, label=L"$$objective function contours")
	contour!(ax3, xs3, ys3, g1, levels= [0], linewidth=2, color=:red)
	contour!(ax3, xs3, ys3, g2, levels= [0], linewidth=2, color=:green)
	#limits!(ax3, -3,1.1,-2, 2)

	
	g1_up = x -> -x^2/4 + 1
	g2_low = x -> 1 + x/2
	xs3_bound = -2:0.1:0; 
	
	band!(ax3, xs3_bound, g1_up.(xs3_bound), g2_low.(xs3_bound), color=(:cyan, 0.7), 
		label = L"$$Feasible region")
	scatter!(ax3, [x1p], [x2p], marker=:cross)
	text!(x1p, x2p - 0.5, text="Optimal Point", color=:blue)
	
	axislegend(ax3, position=:rb)
	fig3
end

# ╔═╡ 5b803afe-0cc5-416e-9fbf-7b34b4de474a
begin
	f4(x1, x2) = (x1-2)^2 + (x2-1)^2
	g11(x1, x2) = x2 - 2x1
	g12(x1, x2) = x1^2 - x2
	g13(x1, x2) = x1 + x2 - 2
	
	fig4 = Figure(size = (600,400))
	ax4 = CairoMakie.Axis(fig4[1,1], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.2))
	limits!(ax4, -1.01,4.01,-1.01, 3.51)
    
	xs4 = -1:0.01:4; ys4 = -1:0.05:3.5
	la4 = collect(-0.5:0.1:5.0)
	contour!(ax4, xs4, ys4, f4, levels = la4, linewidth=1, label=L"$$ Objective function contours")
	contour!(ax4, xs4, ys4, g11, levels= [0], linewidth=2, color=:red)
	contour!(ax4, xs4, ys4, g12, levels= [0], linewidth=2, color=:blue)
	contour!(ax4, xs4, ys4, g13, levels= [0], linewidth=2, color=:green)

	
	g11_up = x -> 2x    # red
	g12_low = x -> x^2 # green
	g13_up = x -> -x + 2  # blue
	xs13_bound = 0:0.1:0.7
	xs14_bound = 0.7:0.1:1
	band!(ax4, xs13_bound, g11_up.(xs13_bound), g12_low.(xs13_bound), color=(:cyan, 0.7), label = "Feasible region")
	band!(ax4, xs14_bound, g12_low.(xs14_bound), g13_up.(xs14_bound), color=(:cyan, 0.7))
	
	axislegend(ax4, position=:rb)
	scatter!(ax4, [2], [1], marker=:cross, color=:red)
	text!(ax4, [2], [1], text="Unconstraint", color=:red)
	scatter!(ax4, [0.8], [1.6], marker=:cross, color=:green)
	text!(ax4, [0.5], [1.8], text="Only EC", color=:green)
	scatter!(ax4, [0.67], [1.33], marker=:cross, color=:blue)
	text!(ax4, [-0.8], [1.33], text="EC and IC", color=:blue)
	fig4
end

# ╔═╡ e4629ef4-845a-4cd7-a206-a5ffc24bf050
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

# ╔═╡ ea7320d4-3ce3-4c4c-a0ac-0e599e099817
md"""
## Reference

[1] R.L. Fox, Optimization Methods in Engineering Design, Addison Wesley, 1971


[2] A. D. Belegundu, and T. R. Chandrupatla , Optimization Concepts and Applications in Engineering, 3rd, 2019
"""

# ╔═╡ 98d333eb-4152-405f-8342-a05d2e4d26f6
TableOfContents(title="Local Descent", depth=4)

# ╔═╡ e6e91723-f6e5-416a-b36c-1247847707ae
html"""
<script>
var section = 0;
var subsection = 0;
var headers = document.querySelectorAll('h2, h3');
for (var i=0; i < headers.length; i++) {
    var header = headers[i];
    var text = header.innerText;
    var original = header.getAttribute("text-original");
    if (original === null) {
        // Save original header text
        header.setAttribute("text-original", text);
    } else {
        // Replace with original text before adding section number
        text = header.getAttribute("text-original");
    }
    var numbering = "";
    switch (header.tagName) {
        case 'H2':
            section += 1;
            numbering = section + ".";
            subsection = 0;
            break;
        case 'H3':
            subsection += 1;
            numbering = section + "." + subsection;
            break;
    }
    header.innerText = numbering + " " + text;
};
</script>
"""

# ╔═╡ Cell order:
# ╠═34cb4d35-319e-48ea-b532-2403d8950301
# ╠═806be89e-167a-11ed-1f7e-87dd862ce225
# ╟─41c4cc45-83f9-47ba-a708-19dca0286201
# ╟─28ca0c07-840f-459b-bfb1-ca06445877d6
# ╠═6c52e6a4-aa96-4c1b-8bed-9afeda5017ad
# ╟─02036da6-4704-48b5-847a-8ed5d45c86ad
# ╟─fa521ffd-f538-4115-ad94-d1b058315687
# ╠═3f7f91f4-1c8a-4276-a403-9620cf086096
# ╟─48e0fe25-9725-4e05-844e-62d7784ca793
# ╟─f63b57b8-a773-4765-b3b7-f6327f9962a2
# ╠═27c448f8-d434-4f77-b0e6-1e9b36c06de8
# ╠═2da61ffe-a205-48c7-b7df-3331882d0ec7
# ╠═ffd53430-b3df-4a31-bd3b-732a548d2f90
# ╠═5b803afe-0cc5-416e-9fbf-7b34b4de474a
# ╠═e4629ef4-845a-4cd7-a206-a5ffc24bf050
# ╟─ea7320d4-3ce3-4c4c-a0ac-0e599e099817
# ╟─98d333eb-4152-405f-8342-a05d2e4d26f6
# ╟─e6e91723-f6e5-416a-b36c-1247847707ae
