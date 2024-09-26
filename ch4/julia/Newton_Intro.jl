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

# ╔═╡ c98bf95e-7ae6-11ef-1efb-ed52cfaa1bce
using Pkg; Pkg.activate()

# ╔═╡ 2317a938-c119-414d-8773-20acf9d4c5e8
begin
	using Printf, PlutoUI
	using ForwardDiff
	using Symbolics, LaTeXStrings
	using CairoMakie 
		set_theme!(theme_latexfonts(), fontsize=22)
end

# ╔═╡ e5d25072-0c4e-4715-bddb-62ecb3138279
begin
	@variables x
		
	f(x) = (x-2)^4 + 2x^2 - 4x + 4
	∇f(x) = ForwardDiff.derivative(f,x)
    ∇2f(x) = ForwardDiff.derivative(∇f,x)
end

# ╔═╡ 4cde65b4-075f-48bf-8c06-db750468e562
md"""
	
 $n$ = $(@bind n PlutoUI.Slider(0:10; show_value=true, default=0))

 $x₀$ = $(@bind x0 PlutoUI.Slider(-10:1:10; show_value=true, default=0))
"""

# ╔═╡ 73ff8404-fcd3-4d4c-8287-6879083f2b02
begin
	global x′ = 3
	@printf("----------------------------------------------\n")
	@printf(" i      xₖ            pₖ             fk\n")
	@printf("----------------------------------------------\n")
	@printf("%2.0d    % 1.4f      % 1.4f       % 1.4f\n", 0, round(x′, digits=4),
		-∇2f(x′)/∇f(x′), f(x′))
	for i = 1:n
		p = -∇2f(x′)\∇f(x′)
		global x′ = x′ + p
		@printf("%2.0d    % 1.4f      % 1.4f        % 1.4f\n", i, round(x′, digits=4),
		round(p, digits=4), f(x′))
	end
end

# ╔═╡ 107b2c28-92bd-4186-8e85-bd1ef266c3c0
straight(x0, y0, x, m) = y0 + m * (x - x0)

# ╔═╡ faf4cd15-7a7d-4ccc-94ba-42a0eda75d01
function standard_Newton1(f, n, x_range, x0, xmin=0, xmax=5, ymin=-10, ymax=10)
	
	f′ = x -> ForwardDiff.derivative(f,x)

	fig2 = Figure(size = (600,400))
	empty!(fig2)
	bx = Axis(fig2[1,1], xlabel = L"x_1", ylabel = L"\nabla f(x)",  
		aspect = AxisAspect(1.5))
	limits!(bx, xmin-0.2, xmax+0.2 ,ymin-0.2, ymax+0.2)

	lines!(bx, x_range, f, linewidth=3)
	lines!(bx, [xmin, xmax],[0,0], color = :magenta, linewidth=3, linestyle=:dash)
	scatter!(bx, [x0], [0], color = :green, markersize=15, strokecolor=:black, 
		strokewidth=1)
	text!(bx, x0, -5, text=L"x_0", fontsize=22)

	for i in 1:n
		lines!(bx, [x0, x0], [0, f(x0)], color=(:black, 0.8), linewidth=2)
		scatter!(bx, [x0], [f(x0)], color=:red, markersize=15, strokecolor=:black,
			strokewidth=1)
		m = f′(x0)

		lines!(bx, x_range, [straight(x0, f(x0), x, m) for x in x_range],
			color=(:blue,0.5), linestyle=:dash, linewidth=2)

		x1 = x0 - f(x0) / m     # update Newton 
		scatter!(bx, [x1], [0], color=:green, markersize=15, strokecolor=:black,
			strokewidth=1)
		text!(bx, x1, -5, text=L"x_%$(i)", fontsize=22)

		x0 = x1
	end
	
	fig2

end

# ╔═╡ 9a5808c1-68c0-407d-998f-ba44778603d2
begin
	xrange = -0:0.01:5
	
	standard_Newton1(∇f, n, xrange, x0, -0, 5, -20, 50)
end

# ╔═╡ Cell order:
# ╠═c98bf95e-7ae6-11ef-1efb-ed52cfaa1bce
# ╠═2317a938-c119-414d-8773-20acf9d4c5e8
# ╠═e5d25072-0c4e-4715-bddb-62ecb3138279
# ╟─4cde65b4-075f-48bf-8c06-db750468e562
# ╠═73ff8404-fcd3-4d4c-8287-6879083f2b02
# ╠═9a5808c1-68c0-407d-998f-ba44778603d2
# ╟─faf4cd15-7a7d-4ccc-94ba-42a0eda75d01
# ╠═107b2c28-92bd-4186-8e85-bd1ef266c3c0
