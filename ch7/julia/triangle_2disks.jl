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

# ╔═╡ ff74a7c6-a306-11ef-0c0e-5d9a25c427b3
using Pkg; Pkg.activate()

# ╔═╡ a82f9681-35f1-4259-be24-75ce6f84946d
begin
	using JuMP
	using Ipopt
	using CairoMakie;
		set_theme!(theme_latexfonts(), fontsize=18)
	using PlutoUI
end

# ╔═╡ e1bfe768-eaed-42bc-9388-2f79f8e9c9a9
md"""

k = $(@bind k PlutoUI.Slider(0:1:10, show_value=true, default=0))

"""

# ╔═╡ afb59030-2451-442f-8476-16b53809b12c
function create_circle(center, radius)
    θ = range(0, 2π, length=100)  # Generate 100 points from 0 to 2π
    Point2f.(center[1] .+ radius .* cos.(θ), center[2] .+ radius .* sin.(θ))
end

# ╔═╡ 5aafcba2-ec7b-4beb-8d8b-aef6a913718d
begin
	# Define the model
	model = Model(Ipopt.Optimizer)
	
	# Variables
	@variable(model, x[1:7])
	
	# Objective function: Minimize the area of the triangle
	@objective(model, Min, 0.5 * abs(x[1] * x[3]))
	
	# Constraints
	@constraint(model, x[1] >= 0)  # Disk 1 above x-axis
	@constraint(model, x[3] >= 0)  # Disk 2 above x-axis
	@constraint(model, (x[4] - x[6])^2 + (x[5] - x[7])^2  - 4 >= 0)  # Disk 1 left of x1
	@constraint(model, x[5] - 1 >= 0)  # Disk 2 left of x1
	@constraint(model, x[7] - 1 >= 0)  # Disk 1 right of (0,0)
	@constraint(model, (x[3]*x[4] - x[2]*x[5])/sqrt(x[3]^2 + x[2]^2) - 1 >= 0)
	@constraint(model, (x[3]*x[6] - x[2]*x[7])/sqrt(x[3]^2 + x[2]^2) - 1 >= 0)
	@constraint(model, (-x[3]*x[4] + (x[2]-x[1])*x[5] + x[1]*x[3])/(sqrt(x[3]^2 + (x[2]-x[1])^2)) -1 >= 0)  # Disk 2 right of (0,0)
	@constraint(model, (-x[3]*x[6] + (x[2]-x[1])*x[7] + x[1]*x[3])/(sqrt(x[3]^2 + (x[2]-x[1])^2)) -1 >= 0)  # Disk 1 right of (x2,x3)

	
	# Initial values
	set_start_value(x[1], 3)
	set_start_value(x[2], 0)
	set_start_value(x[3], 2)
	set_start_value(x[4], 1)
	set_start_value(x[5], 1)
	set_start_value(x[6], 4)
	set_start_value(x[7], 1)

	# Set the maximum number of iterations
	set_optimizer_attribute(model, "max_iter", k)
	
	# Solve the problem
	optimize!(model)
	
	# Get the results
	x_opt = value.(x)
	println("Optimal solution: ", x_opt)
	println("Objective function value at optimal solution: ", objective_value(model))
	
	# Plot the triangle and circles
	fig1 = Figure(size=(600,400), fontsize=18)
	ax1 = Axis(fig1[1,1], xlabel=L"x_1", ylabel=L"x_2", 
		title = L"$$ Optimized Triangle Containing Two Disjoint Disks in %$k Steps", aspect = DataAspect())
	limits!(ax1, -2, 8, -2, 6)
	
	# Plot triangle
	# Define the vertices of the triangle
	vertices = Point2f[(0, 0), (x_opt[1], 0), (x_opt[2], x_opt[3])]
	poly!(ax1, vertices, color = (:yellow, 0.2), strokewidth = 2, 
		strokecolor = :black, label=L"$$Triangle")
	# Define the centers and radius of the circles
	circle1_center = Point2f(x_opt[4], x_opt[5])
	circle2_center = Point2f(x_opt[6], x_opt[7])
	radius = 1.0
	
	
	# Plot the circles
	circle1 = create_circle(circle1_center, radius)
	circle2 = create_circle(circle2_center, radius)
	poly!(ax1, circle1, color = (:blue, 0.2), strokewidth = 2, strokecolor = :blue,
		label=L"$$Disk 1")
	poly!(ax1, circle2, color = (:red, 0.2) , strokewidth = 2, strokecolor = :red,
		label=L"$$Disk 2")

	# Plot points
	scatter!(ax1, [0, x_opt[1], x_opt[2]], [0, 0, x_opt[3]], color = :black, 
		markersize = 8, marker = :circle)
	scatter!(ax1, [x_opt[4]], [x_opt[5]], color = :black, markersize = 10, 
		marker = :circle)
	scatter!(ax1, [x_opt[6]], [x_opt[7]], color = :black, markersize = 10, 
		marker = :circle)
	
	axislegend(ax1, labelsize=12)
	fig1
	
	
end

# ╔═╡ Cell order:
# ╠═ff74a7c6-a306-11ef-0c0e-5d9a25c427b3
# ╠═a82f9681-35f1-4259-be24-75ce6f84946d
# ╟─e1bfe768-eaed-42bc-9388-2f79f8e9c9a9
# ╠═5aafcba2-ec7b-4beb-8d8b-aef6a913718d
# ╠═afb59030-2451-442f-8476-16b53809b12c
