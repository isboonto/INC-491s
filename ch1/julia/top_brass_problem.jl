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

# ╔═╡ 867c6bac-75a5-11f0-2ab7-f3b10d441a44
using Pkg; Pkg.activate();

# ╔═╡ dc92b9f6-4723-44f7-9414-c0df7310cbb8
begin
	using Makie, GLMakie, CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using PlutoUI
end

# ╔═╡ 678708a0-6f5a-4c0e-9e27-1ff298b34eea
using JuMP, GLPK, Gurobi

# ╔═╡ 51e39bc8-7811-47dd-80af-100399dd8842
md"""
# Top Brass Problem
"""

# ╔═╡ b9db6c13-339e-493f-a521-2188138c6daf
md"""
Top Brass Trophy Company makes large championship trophies for youth athletic leagues. Currently, they are planning production for fall sports: US football and football.

Each football trophy has a wood base, an engraved plaque, a large brass football on top, and returns $12 in profit. Soccer trophies are similar, except that a brass soccer ball is on top, and the unit profit is only $9.

Since the football has an asymmetric shape, its base requires four board feet of wood;

The soccer base requires only two board feet. Currently, there are 1000 brass footballs in stock, 1500 soccer balls, 1750 plaques, and 4800 board feet of wood. What should trophies be produced from these supplies to maximize total profit, assuming that all that is made can be sold?
"""

# ╔═╡ 40a1b0ad-299e-48df-afb2-5fd5a62ec763
md"""
Recipe for building each trophy

|          |  wood | plaques | footballs | soccer balls |  profit |
|----------|:-----:|--------:|----------:|-------------:|--------:|
| **US football** |  4 ft |    1    |    1      |      0       |  \$12   |
| **football**   |  2 ft |    1    |    0      |      1       |  \$9    |

Quantity of each ingredient in stock

|          | wood | plaques | footballs |  soccer balls |
|----------|:-----:|-------:|----------:|--------------:|
|**in stock** | 4800 ft | 1750 | 1000 | 1500 |

"""

# ╔═╡ 346f10f4-320c-4a4d-b37a-7303b30e9638
md"""
# Top Brass model components
"""

# ╔═╡ 917db8e2-1115-4d18-9699-d741ed3b3da7
md"""
- **Decision variables**
	. f : number of US football trophies built
	. s : number of football trophies built

- **Constraints**
	. 4f + 2s <= 4800       (wood budget)
    . f + s <= 1750         (plaque budget)
    . 0 <= f <= 1000        (football budget)
    . 0 <= s <= 1500        (soccer ball budget)

- **Objective
	. Maximize 12f + 9s     (profit)
"""

# ╔═╡ b7c3c677-2114-4fa3-a87b-4557a2003c9d
md"""
$\begin{align*}
	&\operatorname*{maximize}\quad f(x) = 12f + 9s \\ 
	&\operatorname{subject to} \quad 4f + 2s <= 4800 \\
	&\hspace{2.95cm} f + s <= 1750 \\
	&\hspace{2.95cm} 0 <= f <= 1000 \\
	&\hspace{2.95cm} 0 <= s <= 1500 
\end{align*}$
"""

# ╔═╡ 4b7230d4-637c-49bc-b99d-07bdd0576636
let
	m = Model(GLPK.Optimizer)
	# Model(Gurobi.Optimizer)
	@variable(m, 0 <= f <= 1000)		# US football trophies
	@variable(m, 0 <= s <= 1500)		# football trophies

	@constraint(m, 4f + 2s <= 4800)		# total board feet of wood
	@constraint(m, f + s <= 1750)		# total number of plagues

	@objective(m, Max, 12f + 9s)		# maximize profit

	# Printing the prepared optimization model
	print(m)

	# Solving the optimization problem
	JuMP.optimize!(m)

	# Printing the optimal solutions obtained
	println("Build ", JuMP.value(f), " US football trophies.")
	println("Build ", JuMP.value(s), " footbal trophies.")
	println("Total profit will be \$", JuMP.objective_value(m))
end

# ╔═╡ c3c42e76-8a92-448b-a6b8-816b0602ccc6
md"""
Return Profit $p$ = $(@bind p PlutoUI.Slider(6000:100:30000; default=17700, show_value =true))
"""

# ╔═╡ 2773ac11-74c5-4521-8e9c-a9832b2833b3
begin
	ra = 27.3
	profit(fo,p) = -(12/9)*fo .+ p/9
	profit_t(fo,p) = (9/12)*fo .- (225/108)*p/ra .+ p/9
end

# ╔═╡ 8fe48c13-a63b-40f9-a46b-c3f9b5becfee
let 
	global fig1 = Figure(); #(size=(600, 400))
	ax1 = CairoMakie.Axis(fig1[1,1], title=L"$$ Geometry of Top Brass", 
						 ylabel=L"Football ($s$)", xlabel=L"US Football ($f$)", 
						 aspect = 1.5)
	ylims!(ax1, (0, 2500))
	xlims!(ax1, (0, 3500))
	

	f = 0:3500; 
	s1 = f-> -2f + 2400             		# number of s depending on f
	s2 = f-> -f + 1750
	s3 = f-> 0f + 1500
	s4 = f-> 0f + 1000
	
	
	# Total board feet of wood
	lines!(ax1, f, s1.(f), 
		   label=L"$$Total board feet of wood", color=:green, linewidth=2)
	band!(ax1, [0, 1200], [0, 0], [2400, 0], color=(:green, 0.1)); # band!(ax1, x, up, lw)

	# Total number of plague
	lines!(ax1,f, s2.(f),
		   label=L"$$Total number of plaques", color=:red, linewidth=2)
	band!(ax1, [0, 1750], [0, 0], [1750, 0], color=(:red, 0.1));

	# US football trophy
	lines!(ax1, f, s3.(f),
		   label=L"$$US Footbal trophies", linestyle = :dashdot, color = :blue, 
		   linewidth =2)

	# Football trophy
	lines!(ax1, 1000*ones(size(f,1)), f,
		  label=L"$$Footbal trophy", linestyle = :dash, color = :blue,
		  linewidth = 2)

	# feasible region
	xr = [0, 250, 650, 1000]
	yr = [1500, 1500, 1100, 400]
	band!(ax1, xr, fill(1.0, length(xr)), yr, color=(:darkblue, 0.3))

	# optimal point
	scatter!(ax1, [650], [1100], marker = :circle, strokecolor = :red, 
			 markersize = 10, strokewidth=1)

	lines!(ax1, f, profit(f,p), linewidth=2, linestyle=:dash, color=:black)
	
	# Define xx and compute y values
	xx = [p/ra, p/ra + 200]
	yy = [profit_t(xx[1], p), profit_t(xx[2], p)]

	# Plot the line with arrow
	lines!(ax1, xx, yy, linewidth = 1, color = :black)
	text!(ax1, xx[1], yy[1]-250, text="Profit = $(p)", rotation= (-51.1)*π/180, 
		  color=:red)


	# Define start and end points
	start = Point2f(p/ra, profit_t(p/ra, p))
	stop = Point2f(p/ra +200, profit_t(p/ra + 200, p))

	# Draw the arrow
	arrows!(ax1, [start], [stop-start], arrowsize=15,
		   arrowhead = :utriangle, linewidth=3, color=:black)
	scatter!(ax1, start[1], start[2], marker = :circle, strokecolor = :red, 
			 markersize = 10, strokewidth=1)

	axislegend(ax1, labelsize=14, rowgap=3)
	
end

# ╔═╡ 19f34b84-9e41-4666-8457-ce1bc41f78f1
fig1

# ╔═╡ 08b74bef-cebb-4f24-b2b7-233a8987f4be
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/top_brass1.pdf", fig1)

# ╔═╡ f2858957-b388-43a4-9a08-8d0986e05022
let 
	global fig2 = Figure(); #(size=(600, 400))
	ax1 = CairoMakie.Axis(fig2[1,1], title=L"$$ Geometry of Top Brass", 
						 ylabel=L"Football ($s$)", xlabel=L"US Football ($f$)", 
						 aspect = 1.5)
	ylims!(ax1, (0, 2500))
	xlims!(ax1, (0, 3500))
	

	f = 0:3500; 
	s1 = f-> -2f + 2400             		# number of s depending on f
	s2 = f-> -f + 1750
	s3 = f-> 0f + 1500
	s4 = f-> 0f + 1000
	
	
	# Total board feet of wood
	lines!(ax1, f, s1.(f), 
		   label=L"$$Total board feet of wood", color=:green, linewidth=2)
	band!(ax1, [0, 1200], [0, 0], [2400, 0], color=(:green, 0.1)); # band!(ax1, x, up, lw)

	# Total number of plague
	lines!(ax1,f, s2.(f),
		   label=L"$$Total number of plaques", color=:red, linewidth=2)
	band!(ax1, [0, 1750], [0, 0], [1750, 0], color=(:red, 0.1));

	# US football trophy
	lines!(ax1, f, s3.(f),
		   label=L"$$US Footbal trophies", linestyle = :dashdot, color = :blue, 
		   linewidth =2)

	# Football trophy
	lines!(ax1, 1000*ones(size(f,1)), f,
		  label=L"$$Footbal trophy", linestyle = :dash, color = :blue,
		  linewidth = 2)

	# feasible region
	xr = [0, 250, 650, 1000]
	yr = [1500, 1500, 1100, 400]
	band!(ax1, xr, fill(1.0, length(xr)), yr, color=(:darkblue, 0.3))

	# optimal point
	scatter!(ax1, [650], [1100], marker = :circle, strokecolor = :red, 
			 markersize = 10, strokewidth=1)

	#lines!(ax1, f, profit(f,p), linewidth=2, linestyle=:dash, color=:black)
	
	# Define xx and compute y values
	xx = [p/ra, p/ra + 200]
	yy = [profit_t(xx[1], p), profit_t(xx[2], p)]

	# Plot the line with arrow
	lines!(ax1, xx, yy, linewidth = 1, color = :black)


	# Define start and end points
	start = Point2f(p/ra, profit_t(p/ra, p))
	stop = Point2f(p/ra +200, profit_t(p/ra + 200, p))

	# Draw the arrow
	arrows!(ax1, [start], [stop-start], arrowsize=15,
		   arrowhead = :utriangle, linewidth=3, color=:black)
	scatter!(ax1, start[1], start[2], marker = :circle, strokecolor = :red, 
			 markersize = 10, strokewidth=1)

	axislegend(ax1, labelsize=14, rowgap=3)

	for p in [11000, 14000, 17700, 21000, 24000] 
		lines!(ax1, f, profit(f,p), linewidth=1, linestyle=:dash, color=:black, label="")
		text!(ax1, p/ra, profit(p/ra,p)-200, 
    		text = "Profit = $(p)", 
    		rotation = (-53.1)*π/180, 
    		color = :red, 
    		fontsize = 14)
	end
	
end

# ╔═╡ 82b58a00-86fa-4abc-bfdb-31a20fef6e60
fig2

# ╔═╡ 18a8f26a-c100-4366-a00e-8600efbf6373
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/top_brass2.pdf", fig2)

# ╔═╡ Cell order:
# ╠═867c6bac-75a5-11f0-2ab7-f3b10d441a44
# ╠═dc92b9f6-4723-44f7-9414-c0df7310cbb8
# ╟─51e39bc8-7811-47dd-80af-100399dd8842
# ╟─b9db6c13-339e-493f-a521-2188138c6daf
# ╟─40a1b0ad-299e-48df-afb2-5fd5a62ec763
# ╟─346f10f4-320c-4a4d-b37a-7303b30e9638
# ╟─917db8e2-1115-4d18-9699-d741ed3b3da7
# ╟─b7c3c677-2114-4fa3-a87b-4557a2003c9d
# ╠═678708a0-6f5a-4c0e-9e27-1ff298b34eea
# ╠═4b7230d4-637c-49bc-b99d-07bdd0576636
# ╠═8fe48c13-a63b-40f9-a46b-c3f9b5becfee
# ╟─c3c42e76-8a92-448b-a6b8-816b0602ccc6
# ╠═2773ac11-74c5-4521-8e9c-a9832b2833b3
# ╠═19f34b84-9e41-4666-8457-ce1bc41f78f1
# ╠═08b74bef-cebb-4f24-b2b7-233a8987f4be
# ╠═f2858957-b388-43a4-9a08-8d0986e05022
# ╠═82b58a00-86fa-4abc-bfdb-31a20fef6e60
# ╠═18a8f26a-c100-4366-a00e-8600efbf6373
