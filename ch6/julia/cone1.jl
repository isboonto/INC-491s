### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 60eeca4c-904a-11ef-1de9-c3667533212e
begin
	using Pkg; Pkg.activate()
end

# ╔═╡ 6373d11c-92a9-434a-95a6-f405b1731d2d
begin
	using ForwardDiff 
	using LaTeXStrings
	using CairoMakie
		set_theme!(theme_latexfonts(), fontsize=18)
	using LinearAlgebra
end

# ╔═╡ 6d9c1bb4-81fa-4883-93da-64b36ebffd57
begin
	f(x) = x[1]^2  + 2x[2]^2 
	∇f(x)= ForwardDiff.gradient(f, x)
	f1(x,y) = f([x,y])
	gf(x,y) = ∇f([x,y])
end

# ╔═╡ 44e8859d-6bee-450e-be2c-77e6570131fc
gf(-1,1)'*[1,1]

# ╔═╡ 146d5f22-62d4-4876-9fb5-92694997a554
begin
	x1s = LinRange(-2,4, 50); x2s = LinRange(-2,3, 50); 

	xv = [-1,1]
	
	z = [f1(x,y) for x in x1s, y in x2s];
	zl = [f1(-1,1) + gf(-1,1)'*[x,y] for x in x1s, y in x2s ]
	
	fig = Figure(size = (400,400))
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = AxisAspect(1.4))
	hidedecorations!(ax)  # hides ticks, grid and lables
	hidespines!(ax)  # hide the frame
	
	contour!(ax,x1s,x2s,z, levels= 0:1:20, color=(:blue, 0.3), linewidth =1)
	contour!(ax,x1s,x2s,zl, levels = [9], color=(:blue), linewidth=3, linestyle=:solid)
	limits!(ax, -2,0.5,0.5, 2)

	scatter!([-1],[1], color=:blue)
	yy = x -> (1.75-0.5)/(0.5-(-2.0))*x + 1.5  # offset
	
	band!(ax, x1s, yy.(x1s), 0.5*ones(length(x1s)), color=(:cyan, 0.2))
	fig[1,1] = ax
	

	
	pv = [1/√5, 2/√5] # steepest direction is -g 
	α = 1
	cosv = xv'*pv / (norm(xv)*norm(pv))
	alength = xv'*pv/norm(pv)
	g1 = gf(xv[1],xv[2])*0.1

	
	text!(ax, L"\nabla f^T\mathbf{p}  = 0", position=(-0.45,1.15), color=:blue)
	
	arrows!(ax, [xv[1]], [xv[2]], [g1[1]], [g1[2]], lengthscale = 1, 
		color=:black, linewidth=2, arrowsize=15)
	text!(ax, L"\nabla f", position=(-1.4,1.4), color=:black)
	
	arrows!(ax, [xv[1]], [xv[2]], 0.5*[pv[1]], [pv[2]], lengthscale = 0.5,
		color=:red, linewidth=2, arrowsize=15)
	text!(ax, L"\nabla f^T\mathbf{p} > 0", position=(-0.85,1.5), color=:red)

	arrows!(ax, [xv[1]], [xv[2]], [-1], [-1.5], lengthscale = 0.3, 
		color=:red, linewidth=2, arrowsize=15)
	text!(ax, L"\nabla f^T\mathbf{p} < 0", position=(-1.2,0.55), color=:red)

	
	
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
	
	save("first_con1b.pdf", fig)
	run(`pdfcrop  --clip first_con1b.pdf first_con1.pdf`)
	fig
end

# ╔═╡ Cell order:
# ╠═60eeca4c-904a-11ef-1de9-c3667533212e
# ╠═6373d11c-92a9-434a-95a6-f405b1731d2d
# ╠═6d9c1bb4-81fa-4883-93da-64b36ebffd57
# ╠═44e8859d-6bee-450e-be2c-77e6570131fc
# ╠═146d5f22-62d4-4876-9fb5-92694997a554
