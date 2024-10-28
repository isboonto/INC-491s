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
		set_theme!(theme_latexfonts(), fontsize=18, font="amsmath")
	using LinearAlgebra
	#using MathTeXEngine
	#to_font("New Computer Modern")
end

# ╔═╡ 6d9c1bb4-81fa-4883-93da-64b36ebffd57
begin
	f(x) = x[1] + 2x[2] 
	g1(x) = (1/4)*x[1]^2 + x[2]^2 - 1;
	g2(x) = -x[2];
	#hx2 = x -> sqrt(1 - (1/4)*x^2)
	
	∇f(x)= ForwardDiff.gradient(f, x)
	∇g1(x)= ForwardDiff.gradient(g1, x)
	∇g2(x)= ForwardDiff.gradient(g2, x)
	f1(x,y) = f([x,y])
	g1(x,y) = g1([x,y])
	g2(x,y) = g2([x,y])
	gf(x,y) = ∇f([x,y])
	gg1(x,y) = ∇g1([x,y])
	gg2(x,y) = ∇g2([x,y])
end

# ╔═╡ 146d5f22-62d4-4876-9fb5-92694997a554
begin
	x1s = LinRange(-3.1,4, 50); x2s = LinRange(-3.1,3.1, 50); 

	# case 1
	xB = [-2 0; 2 0; √2 √2/2] # xA xC xB
	#xA = [-√2,-√2/2]
	#xB = -xA
	
	z = [f1(x,y) for x in x1s, y in x2s];
	
	
	fig = Figure(size = (600,400))
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", 
		aspect = AxisAspect(1.2))
	hidedecorations!(ax, ticklabels=false, ticks=false, grid=true, label=false) 

	g11 = x -> -sqrt(-(1/4)*x^2 + 1)  
	x11s = -2:0.05:2
	band!(ax, x11s,  -3*ones(length(x11s)), g11.(x11s), color=(:cyan, 0.4))
	band!(ax, x11s, g11.(x11s), 0*ones(length(x11s)) , color=(:cyan, 0.4))
	band!(ax, x11s, -g11.(x11s), 3*ones(length(x11s)), color=(:cyan, 0.4))
	x12s = -3:0.05:-2
	x13s = 2:0.05:3
	band!(ax, x12s, -3*ones(length(x12s)), 3*ones(length(x12s)), color=(:cyan, 0.4))
	band!(ax, x13s, -3*ones(length(x13s)), 3*ones(length(x13s)), color=(:cyan, 0.4))
		
	contour!(ax,x1s,x2s,z, levels= -20:1:20, color=(:blue, 0.4), linewidth =1)
	contour!(ax,x1s,x2s,g1, levels=[0], color=(:red, 1.0), linewidth=1)
	contour!(ax,x1s,x2s,g2, levels=[0], color=(:red, 1.0), linewidth=1)
	limits!(ax, -3,3,-2.4, 2.4)

	
	gfB1 = gf(xB[1,1],xB[1,2])*0.5
	ggB11 = gg1(xB[1,1],xB[1,2])*0.5
	ggB12 = gg2(xB[1,1],xB[1,2])*0.6
	
	gfB2 = gf(xB[2,1],xB[2,2])*0.5
	ggB21 = gg1(xB[2,1], xB[2,2])*0.5
	ggB22 = gg2(xB[2,1], xB[2,2])*0.6

	gfB3 = gf(xB[3,1],xB[3,1])*0.5
	ggB31 = gg1(xB[3,1], xB[3,2])*0.5
	ggB32 = gg2(xB[3,1], xB[3,2])*0.6
	

	# xA related plot
	arrows!(ax, [xB[1,1]], [xB[1,2]], [gfB1[1]], [gfB1[2]], lengthscale = 1, 
		color=:blue, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla f", position=(xB[1,1] + 0.1, xB[1,2] + 1.2), color=:blue)
	
	arrows!(ax, [xB[1,1]], [xB[1,2]], [ggB11[1]], [ggB11[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla g_1", position=(ggB11[1]-2.3, ggB11[2] + 0.1), color=:red)

	arrows!(ax, [xB[1,1]], [xB[1,2]], [ggB12[1]], [ggB12[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla g_2", position=(ggB12[1]-2.2, ggB12[2] - 0.5), color=:red)
	scatter!(xB[1,1], xB[1,2], color=:blue)
	text!(ax, L"\mathbf{x}^\ast", position=(xB[1,1] + 0.25
		, xB[1,2]-0.3), color=:blue)

	# xC related plot
	arrows!(ax, [xB[2,1]], [xB[2,2]], [gfB2[1]], [gfB2[2]], lengthscale = 1, 
		color=:blue, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla f", position=(xB[2,1] + 0.1, xB[2,2] + 1.2), color=:blue)
		
	arrows!(ax, [xB[2,1]], [xB[2,2]], [ggB21[1]], [ggB21[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla g_1", position=(ggB21[1]+1.9, ggB21[2] + 0.1), color=:red)

	arrows!(ax, [xB[2,1]], [xB[2,2]], [ggB22[1]], [ggB22[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla g_2", position=(ggB22[1] + 1.9, ggB22[2] - 0.5), color=:red)
	
	scatter!(xB[2,1], xB[2, 2], color=:blue)
	text!(ax, L"\mathbf{x}_C", position=(xB[2,1] + 0.15
		, xB[2,2]-0.3), color=:blue)
	
	# xB related plot
	arrows!(ax, [xB[3,1]], [xB[3,2]], [gfB3[1]], [gfB3[2]], lengthscale = 1, 
		color=:blue, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla f", position=(xB[3,1] - 0, xB[3,2] + 1.1), color=:blue)
	
	arrows!(ax, [xB[3,1]], [xB[3,2]], [ggB31[1]], [ggB31[2]], lengthscale = 1, 
		color=:red, linewidth=1, arrowsize=10)
	text!(ax, L"\nabla g_1", position=(ggB31[1]+0.8, ggB31[2]+0.6), color=:red)

	#arrows!(ax, [xB[3,1]], [xB[3,2]], [ggB32[1]], [ggB32[2]], lengthscale = 1, 
	#	color=:red, linewidth=1, arrowsize=10)
	#text!(ax, L"\nabla g_2", position=(ggB32[1]+0.8, ggB32[2]+0.6), color=:red)
	
	scatter!(xB[3,1], xB[3, 2], color=:blue)
	text!(ax, L"$\mathbf{x}_B$", position=(xB[3,1]-0.45, xB[3,2]-0.2), color=:blue)
	
	text!(ax, L"$$ Minimum", position=(xB[1,1]+0.2,xB[1,2]), color=:black)
	
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
	
	save("eq_con3b.pdf", fig)
	run(`pdfcrop  --clip eq_con3b.pdf eq_con3.pdf`)
	fig
end

# ╔═╡ Cell order:
# ╠═60eeca4c-904a-11ef-1de9-c3667533212e
# ╠═6373d11c-92a9-434a-95a6-f405b1731d2d
# ╠═6d9c1bb4-81fa-4883-93da-64b36ebffd57
# ╠═146d5f22-62d4-4876-9fb5-92694997a554
