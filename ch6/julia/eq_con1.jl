### A Pluto.jl notebook ###
# v0.20.1

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
	#using MathTeXEngine
	#to_font("New Computer Modern")
end

# ╔═╡ 6d9c1bb4-81fa-4883-93da-64b36ebffd57
begin
	f(x) = x[1] + 2x[2] 
	h(x) = (1/4)*x[1]^2 + x[2]^2 - 1;
	#hx2 = x -> sqrt(1 - (1/4)*x^2)
	
	∇f(x)= ForwardDiff.gradient(f, x)
	∇h(x)= ForwardDiff.gradient(h, x)
	f1(x,y) = f([x,y])
	h1(x,y) = h([x,y])
	gf(x,y) = ∇f([x,y])
	gh(x,y) = ∇h([x,y])
end

# ╔═╡ 44e8859d-6bee-450e-be2c-77e6570131fc
gf(-1,1)'*[1,1]

# ╔═╡ 146d5f22-62d4-4876-9fb5-92694997a554
begin
	x1s = LinRange(-3.1,4, 50); x2s = LinRange(-3.1,3.1, 50); 

	xA = [-√2,-√2/2]
	xB = -xA
	
	z = [f1(x,y) for x in x1s, y in x2s];
	
	
	fig = Figure(size = (600,400))
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.2))
	hidedecorations!(ax, ticklabels=false, ticks=false, grid=true, label=false) 
		
	contour!(ax,x1s,x2s,z, levels= -20:1:20, color=(:blue, 0.4), linewidth =1)
	contour!(ax,x1s,x2s,h1, levels=[0], color=(:red, 1.0), linewidth=1)
	limits!(ax, -3,3,-2.4, 2.4)

	
	gA = gf(xA[1],xA[2])*0.6
	hA = gh(xA[1],xA[2])*0.4
	gB = gf(xB[1],xB[2])*0.6
	hB = gh(xB[1],xB[2])*0.4

	
	
	
	arrows!(ax, [xA[1]], [xA[2]], [gA[1]], [gA[2]], lengthscale = 1, 
		color=:blue, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla f", position=(xA[1] + 0.6, xA[2] + 0.7), color=:blue)
	text!(ax, L"\mathbf{x}_A", position=(xA[1]-0.65, xA[2]-0.2), color=:blue)
	
	arrows!(ax, [xA[1]], [xA[2]], [hA[1]], [hA[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla h", position=(hA[1]-1.2, hA[2]-1), color=:red)

	arrows!(ax, [xB[1]], [xB[2]], [gB[1]], [gB[2]], lengthscale = 1, 
		color=:blue, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla f", position=(xB[1] + 0.6, xB[2] + 0.7), color=:blue)
	text!(ax, L"$\mathbf{x}_B$", position=(xB[1]-0.65, xB[2]-0.2), color=:blue)

	arrows!(ax, [xB[1]], [xB[2]], [hB[1]], [hB[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla h", position=(hB[1]+0.7, hB[2]+0.6), color=:red)

	text!(ax, L"$$ Maximum", position=(xB[1]+0.2,xB[2]), color=:black)
	text!(ax, L"$$ Minimum", position=(xA[1]+0.2,xA[2]), color=:black)
	
	scatter!(xA[1], xA[2], color=:blue)
	scatter!(xB[1], xB[2], color=:blue)
	
	
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
	
	#save("eq_con1b.pdf", fig)
	#run(`pdfcrop  --clip eq_con1b.pdf eq_con1.pdf`)
	fig
end

# ╔═╡ Cell order:
# ╠═60eeca4c-904a-11ef-1de9-c3667533212e
# ╠═6373d11c-92a9-434a-95a6-f405b1731d2d
# ╠═6d9c1bb4-81fa-4883-93da-64b36ebffd57
# ╠═44e8859d-6bee-450e-be2c-77e6570131fc
# ╠═146d5f22-62d4-4876-9fb5-92694997a554
