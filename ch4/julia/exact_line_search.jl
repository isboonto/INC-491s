### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 9ff12418-7f2a-11f0-17c4-4501831aeb7e
using Pkg; Pkg.activate()

# ╔═╡ 31560044-2f92-4b69-8165-c17a68084ef0
using CairoMakie; set_theme!(theme_latexfonts(), fontsize=22)

# ╔═╡ d4c3579c-4c27-4200-9782-bff5dbceebc2
using Symbolics, NLsolve, LinearAlgebra, LaTeXStrings, ForwardDiff

# ╔═╡ b7e2d337-01cb-4029-959c-003a39a004a9
using Optim

# ╔═╡ da97e57b-fae6-44a1-8ec9-1452c93a79eb
f = x -> x^2 + 54/x

# ╔═╡ d82a03f7-e717-4945-b027-5ef6da7dc13e
function bracket_minimum(f, x=0, s=1e-2, k = 2.0)
	# Ensure we start by going downhill
	a, ya = x, f(x)
	b, yb = a + s, f(a + s)

	if yb > ya 			# swap the direction
		a, b = b, a
		ya, yb = yb, ya
		s = -s
	end
	while true
		c, yc = b + s, f(b + s)
		if yc > yb
			return a < c ? (a, c) : (c,a)
		end
		a, ya , b, yb = b, yb, c, yc
		s *= k
	end

end

# ╔═╡ 2ddbb3fe-b54c-4e08-8a84-331735678624
a, b = bracket_minimum(f)

# ╔═╡ 58fef3de-e88a-4dba-a37d-c9138bda81f7
md"""
## Exact Line search using NLsolve
"""

# ╔═╡ e6c9a8bc-15cf-4979-a096-e0d4e0d969b0
let
	f = α -> -cos(2 - α) - 2 * exp(5 - 2 * α) + 1
	
	sol = nlsolve(x -> [f(x[1])], [3.5])  # Initial guess between 3 and 4
	sol.zero
end

# ╔═╡ 7dcbf085-b38c-4f11-9952-69a085d0dc28
let
# using brent method from Optim.jl
	f = α -> sin(2-α) + exp(5 - 2α) + α - 3
	#@variables α

	#∇f = ForwardDiff.derivative(f, α)
	#println(∇f)
	
	res = Optim.optimize(f, 0, 4)   # minimize
	sol = Optim.minimizer(res)
end

# ╔═╡ 082eaeb9-34d3-4805-a586-c8504e7b1553
md"""
## Bisection Search
"""

# ╔═╡ f6c193fc-6c6f-4578-9ff5-b6a97cfc4ced
function plot_track(row, col, fig, f, x1, x2, track)
	ax1 = CairoMakie.Axis(fig[row, col], xlabel=L"x")
	
	xs = x1:0.1:x2
	∇f = ForwardDiff.derivative.(f, xs)
	lines!(ax1, xs, ∇f, linewidth=3) # gradient function
	f_max = maximum(f.(xs))
	f_min = minimum(f.(xs))

	band!(ax1, [track[1], track[2]], f_max, f_min, color=(:red, 0.5), 
		  label="Step $((row-1)*3 + col)")

	axislegend(ax1, labelsize=18)
	
	fig
end

# ╔═╡ 216d3cd6-e1c5-42c3-bde9-5664325e7655
function bisection(f, a₀, b₀, ϵ)

    function D(f, a)
        # Approximate the first derivative using central differences
        h = 0.001
        return (f(a + h) - f(a - h)) / (2 * h)
    end

    a = a₀
    b = b₀

    while (b - a) > ϵ
        c = (a + b) / 2.0

        if D(f, c) > 0
            b = c
        else
            a = c
        end
    end

    return (a + b) / 2 # this was changed
end

# ╔═╡ 37f66ca8-9fa7-4fd3-b38c-2a82bfbdb906
function bisection_track(f, a₀, b₀, ϵ)

    function D(f, a)
        # Approximate the first derivative using central differences
        h = 0.001
        return (f(a + h) - f(a - h)) / (2 * h)
    end
	track = [a₀, b₀]
    a = a₀; b = b₀; 
    while (b - a) > ϵ
        c = (a + b) / 2.0

        if D(f, c) > 0
            b = c
        else
            a = c
        end
		track = hcat(track, [a; b])
    end
	
    return track
end

# ╔═╡ 7e1e8285-ea44-48d7-9178-68f5b4e186e3
begin
	track = Any[2,]
	f1 = x -> (1/4)*(sin(x) + sin(x/2))
	∇f1 = x -> ForwardDiff(f1, x)
	x1 = -5; x2 = 9
	#sol = bisection(f1, x1, x2, 5e-2) 
	track = bisection_track(f1, x1, x2, 5e-2)
end

# ╔═╡ 58fc7267-4cab-4c2f-8b5e-51173618a840
begin
	global fig1 = Figure(size=(400*3, 200*size(track,2)/2)); nothing
end

# ╔═╡ 0ad76ceb-0906-4ee7-9b57-079952f32cda
# ╠═╡ disabled = true
#=╠═╡
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/bisection_ex.pdf", fig1)
  ╠═╡ =#

# ╔═╡ 1ea621ee-dc28-4cc8-8538-066167bad91d
let
	empty!(fig1)
	ncols = 4
	#∇f1 = x -> ForwardDiff.gradient(f1, x)
	
	for ii in 1:size(track,2)
		row = div(ii - 1, ncols) + 1
		col = mod(ii - 1, ncols) + 1
		println("Plot at: row = $row, col = $col")
		plot_track(row,col, fig1, f1, x1, x2, track[:,ii])
	end

	fig1
end

# ╔═╡ 4d4694a9-ac83-41dc-b325-0776e528fbe5
md"""
## Approximation Line Search
"""

# ╔═╡ 23caba9d-aa29-4c52-9b23-0f733b788bab
function backtracking_line_search(f, ∇f, x, d, α; p=0.5, β=1e-4)
	y, g = f(x), ∇f(x) 
	while f(x + α*d) > y + β*α*(g'*d)
		α *= p
	end
	
	return α
end

# ╔═╡ e46b6a95-d0ad-466e-822d-3fadd50cf225
begin
	@variables x₁, x₂, α
	
	f2 = x -> 0.1x[1]^6 - 1.5x[1]^4 + 5x[1]^2 + 0.1x[2]^4 + 3x[2]^2 - 9x[2] + 0.5x[1]*x[2]
	
	f2([x₁, x₂])
	latexstring("f(x) = $(f2([x₁, x₂]))")
	∇f = x -> ForwardDiff.gradient(f2, x)
end

# ╔═╡ cf279b71-c07d-4951-8d26-cd5a69f72f4d
ForwardDiff.derivative(f, α)

# ╔═╡ c951dd9a-e46e-43e8-a5c6-8dee3d58db81
fig2 = Figure(size=(800,350)); nothing

# ╔═╡ c12181cc-5fad-48c7-8703-e6d660fa7a0e
function plot_approx(fig, row, col, α)
	α1 = 0:0.001:1.2	
	x0 = [-1.25, 1.25]
	p = [4, 0.75]
	β = 1e-4
	z1 = [f2(x0 + α*p) for α in α1]
	z2 = [f2(x0) + α*1*∇f(x0)'*p for α in α1]
	z3 = [f2(x0) + α*β*∇f(x0)'*p for α in α1]
	
	title_text = latexstring(" α = $(α)")
	ax1 = CairoMakie.Axis(fig[row, col], title=title_text)
	
	lines!(ax1, α1, z1)
	lines!(ax1, α1, z2, color=:blue, linewidth=3)
	lines!(ax1, α1, z3, color=:red, linewidth=3)
	ρ = 0.7
	α2 = α
	while f2(x0 + α2*p) > f2(x0) + β*α2*∇f(x0)'*p 
		
		y1 = f2(x0); y2 = f2(x0 + α2*p);
		x1 = 0; x2 = α2
		formatted = string(round(α2, digits = 2))  
		scatterlines!(ax1, [x1, x2], [y1,y2], linewidth=2, linestyle=:dash, 
					  label="α = $(formatted)", markersize=15, strokewidth=2, strokecolor=:gray)
		α2 *= ρ
	end
	# plot last marker
	y1 = f2(x0); y2 = f2(x0 + α2*p);
	x1 = 0; x2 = α2
	formatted = string(round(α2, digits = 2))  
	scatterlines!(ax1, [x1, x2], [y1,y2], linewidth=2, linestyle=:dash, 
					  label="α = $(formatted)", markersize=15, strokewidth=2, strokecolor=:gray)
	
	axislegend(ax1, labelsize=14, position=:lt)

	fig
end

# ╔═╡ f9232a2d-e448-48da-b09f-a76f3293b23f
begin
	empty!(fig2)
	plot_approx(fig2, 1, 1, 1.2)
	plot_approx(fig2, 1, 2, 0.05)
	# The initial α is very important, too low, the algorithm will stop fast, but miss the lower point. 
end

# ╔═╡ c08249ee-f966-48c6-bb74-e93c069f9b91
# ╠═╡ disabled = true
#=╠═╡
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/approximate_linesearch.pdf", fig2)
  ╠═╡ =#

# ╔═╡ 45e1154e-682c-4658-8e1e-a5aa5c6b07a3
begin
	function check_strong_wolfe(f, ∇f, x, d, α; β = 1e-4, σ = 0.9)
		y, g = f(x), ∇f(x)
		yα = f(x + α*d); gα = ∇f(x + α*d)
		
		armijo = yα ≤ y + β*α*g'*d
		curvature = abs(gα'*d) ≤ σ*a*abs(g'*d)
	
		return armijo && curvature
	end
	
	function strong_wolfe_line_search(f, ∇f, x, d; α = 1.0, β = 1e-4, σ = 0.8, ρ = 0.5, α_min = 1e-8)
		α1 = α
		while α1 > α_min
			if check_strong_wolfe(f, ∇f, x, d, α1; β = β, σ = σ)
				return α
			end
			α *= ρ # shrink step size
		end
		return α_min # fallback
	end
	
end

# ╔═╡ 51e7004a-82fe-4446-97d1-bc5461981a49
begin
	f3 = x -> x[1]^2 + x[1]*x[2] + x[2]^2
	d = [4, 0.75]
	∇f3 = x -> ForwardDiff.gradient(f3, x)
end

# ╔═╡ 381bccc4-9777-4051-8dc3-d611d6c87201
fig3 = Figure(size=(500, 350)); nothing

# ╔═╡ 7184e7da-a665-435c-a7bd-476fe3def4b5
function plot_strong_wolfe(fig, row, col, α_init)
	α1 = 0:0.001:α_init	
	x0 = [1, 2]
	d = [-1, -1]
	β = 1e-4
	z1 = [f3(x0 + α*d) for α in α1]
	z2 = [f3(x0) + α*1*∇f3(x0)'*d for α in α1]
	z3 = [f3(x0) + α*β*∇f3(x0)'*d for α in α1]

	title_text = latexstring(" α = $(α_init)")
	ax1 = CairoMakie.Axis(fig[row, col], title=title_text)
	
	lines!(ax1, α1, z1)
	#lines!(ax1, α1, z2, color=:blue, linewidth=3)
	lines!(ax1, α1, z3, color=:red, linewidth=3)

	ρ = 0.5
	α2 = α_init; α_min = 1e-1
	
	
	while α2 > α_min 
		y1 = f3(x0); y2 = f3(x0 + α2*d);
		x1 = 0; x2 = α2
		formatted = string(round(α2, digits = 2))  
		scatterlines!(ax1, [x1, x2], [y1,y2], linewidth=2, linestyle=:dash, 
				label="α = $(formatted), x = $(x0+α2*d), $(f3(x0+α2*d))", markersize=15, strokewidth=2, strokecolor=:gray)	
		print(x0 + α2*d); println("  ", f3(x0 + α2*d))
		if check_strong_wolfe(f3, ∇f3, x0, d, α2; β = 1e-4, σ = 0.9)
			break
		end
		α2 *= ρ
	
	end
	
	axislegend(ax1, labelsize=14, position=:lt)
	fig
	
end

# ╔═╡ 0234897d-f756-4351-b7e0-d761631e40b5
begin
	empty!(fig3)
	plot_strong_wolfe(fig3, 1, 1, 10)
	
end

# ╔═╡ b17e6730-9564-4959-80ec-10390e538555
# ╠═╡ disabled = true
#=╠═╡
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/approximate_strong_wolfe.pdf", fig3)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═9ff12418-7f2a-11f0-17c4-4501831aeb7e
# ╠═31560044-2f92-4b69-8165-c17a68084ef0
# ╠═da97e57b-fae6-44a1-8ec9-1452c93a79eb
# ╠═2ddbb3fe-b54c-4e08-8a84-331735678624
# ╠═d82a03f7-e717-4945-b027-5ef6da7dc13e
# ╠═58fef3de-e88a-4dba-a37d-c9138bda81f7
# ╠═d4c3579c-4c27-4200-9782-bff5dbceebc2
# ╠═e6c9a8bc-15cf-4979-a096-e0d4e0d969b0
# ╠═cf279b71-c07d-4951-8d26-cd5a69f72f4d
# ╠═b7e2d337-01cb-4029-959c-003a39a004a9
# ╠═7dcbf085-b38c-4f11-9952-69a085d0dc28
# ╠═082eaeb9-34d3-4805-a586-c8504e7b1553
# ╠═7e1e8285-ea44-48d7-9178-68f5b4e186e3
# ╠═58fc7267-4cab-4c2f-8b5e-51173618a840
# ╠═f6c193fc-6c6f-4578-9ff5-b6a97cfc4ced
# ╠═1ea621ee-dc28-4cc8-8538-066167bad91d
# ╠═0ad76ceb-0906-4ee7-9b57-079952f32cda
# ╠═216d3cd6-e1c5-42c3-bde9-5664325e7655
# ╠═37f66ca8-9fa7-4fd3-b38c-2a82bfbdb906
# ╟─4d4694a9-ac83-41dc-b325-0776e528fbe5
# ╠═23caba9d-aa29-4c52-9b23-0f733b788bab
# ╠═e46b6a95-d0ad-466e-822d-3fadd50cf225
# ╠═c951dd9a-e46e-43e8-a5c6-8dee3d58db81
# ╠═c12181cc-5fad-48c7-8703-e6d660fa7a0e
# ╠═f9232a2d-e448-48da-b09f-a76f3293b23f
# ╠═c08249ee-f966-48c6-bb74-e93c069f9b91
# ╠═45e1154e-682c-4658-8e1e-a5aa5c6b07a3
# ╠═51e7004a-82fe-4446-97d1-bc5461981a49
# ╠═381bccc4-9777-4051-8dc3-d611d6c87201
# ╠═7184e7da-a665-435c-a7bd-476fe3def4b5
# ╠═0234897d-f756-4351-b7e0-d761631e40b5
# ╠═b17e6730-9564-4959-80ec-10390e538555
