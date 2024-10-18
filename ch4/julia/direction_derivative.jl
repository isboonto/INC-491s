### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ 933d22e8-f7b1-11ec-0597-b3d1bd741c42
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	set_theme!(theme_latexfonts())
	using Symbolics, ForwardDiff,LinearAlgebra
	using PlutoUI
end

# ╔═╡ 541e270e-8d76-4c80-a4cf-c40fd846ec5f
using Optim

# ╔═╡ 0e3e0db4-de9c-4e2a-b989-7e3dfed2f839
md"""
## Example 4.3 Directional derivative of a quadratic function

Consider the following function of two vriables:

$f(x_1,x_2) = x_1^2 + 2x_2^2 - x_1x_2$
"""

# ╔═╡ d6af1621-ffb8-4a72-86bb-7faa0b2292df
begin
	@variables x₁ x₂

	f(x) = x[1]^2 + 2x[2]^2 - x[1]*x[2]
end

# ╔═╡ 343017fd-4f2c-40bc-beea-9ee054832fe6
md"""
The gradient can be obtained using symbolic differentiation, yielding

$∇f(x_1, x_2) = \begin{bmatrix}2x_1 - x_2\\ 4x_2 - x_1\end{bmatrix}$
"""

# ╔═╡ 1d69b7a5-b7c1-456c-b6b4-9df161e61458
begin
	∇f(x) = ForwardDiff.gradient(f,x)
	∇f([x₁,x₂])   	# Check
end

# ╔═╡ 1868c80b-985d-401b-9f6e-25118474eba7
begin
	# definef functions for plot
	f1(x,y) = f([x,y])		# original function
	gf(x,y) = ∇f([x,y])		# gradient function 
end

# ╔═╡ c96c543e-51e7-4e42-81b7-eb251ebc72c7
let
	x1s = LinRange(-2,4, 50)
	x2s = LinRange(-2,3, 50)
	z = [f1(x,y) for x in x1s, y in x2s]

	fig = Figure(size = (600,400), font = "CMU Serif")
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = AxisAspect(1))
	contour!(ax,x1s,x2s,z, levels= 0:0.7:20, color=(:blue, 0.3), linewidth =1)
	limits!(ax, -2,0,0.5, 2.4)

	scatter!([-1],[1], color=:red)

	xv = [-1,1]
	pv = [1/√5, 2/√5] # steepest direction is -g 
	α = 1
	cosv = xv'*pv / (norm(xv)*norm(pv))
	alength = xv'*pv/norm(pv)
	g1 = gf(xv[1],xv[2])*0.1
    	
	arrows!(ax, [xv[1]], [xv[2]], [α*pv[1]], [α*pv[2]], color=:blue, linewidth=2)
	arrows!(ax, [xv[1]], [xv[2]], [g1[1]], [g1[2]], lengthscale = 1, 
		color=:black, linewidth=2)
	arrows!(ax, [xv[1]], [xv[2]], [pv[1]], [pv[2]], lengthscale = 1,
		color=:red, linewidth=2)
	lines!(ax, [xv[1] + g1[1], xv[1] + pv[1]*alength], 
		[xv[2]+g1[2], xv[2] + pv[2]*alength], linestyle=:dash, color=:black)
	scatter!(ax, xv[1]+pv[1]*alength, xv[2] + pv[2]*alength)
	
	text!(ax, L"\nabla f", position=(-1.6,1.3), color=:black)
	text!(ax, L"\mathbf{x}", position=(-0.98,0.85), color=:red)
	#text!(ax, L"x + \alpha p", position=(-0.3,0.73), color=:black)
	text!(ax, L"\mathbf{p}", position=(-0.5,1.89), color=:red)
	text!(ax, L"\nabla f^T\mathbf{p}", position=(-0.7, 1.3), color=:red)
	
	fig[1,1] = ax
	
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
	
	save("gradient3b.pdf", fig)
	run(`pdfcrop  --clip gradient3b.pdf gradient3.pdf`)
	fig
end

# ╔═╡ d9c830d6-ec91-4722-a554-6458cfdaa7ff
md"""
At point $x=[-1,1]$, the gradient is 

$\nabla f(-1,1) = \begin{bmatrix}-3 \\ 5\end{bmatrix}$
"""

# ╔═╡ 4c9f0680-e805-4ecb-94b4-3ebdc3c76ad3
gf(-1,1)

# ╔═╡ 22e69033-c1a2-4408-b4bf-e2117017d094
md"""
Taking the derivative in the normalized direction $p = [2/\sqrt{5}, -1/\sqrt{5}]$, we obtain

$\nabla f^Tp = \begin{bmatrix}-3 & 5\end{bmatrix}\begin{bmatrix}2/\sqrt{5}\\ -1/\sqrt{5}\end{bmatrix} = -\frac{11}{\sqrt{5}}$

"""

# ╔═╡ eabce287-97b4-4232-8745-3c0c2f81ea59
begin
	g1 = gf(-1,1)
	p = [2/√5, -1/√5]
	g1'*p          		# -11/sqrt{5}
end

# ╔═╡ a0ccee19-b92c-4bd9-a800-c313a4161726
let
	x1s = LinRange(-2,2, 50)
	x2s = LinRange(-2,2, 50)
	z = [f1(x,y) for x in x1s, y in x2s]

	fig = Figure(size = (600,400), font = "CMU Serif")
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = DataAspect())
	contour!(ax,x1s,x2s,z, levels= 0:0.3:10, color=(:blue, 0.3), linewidth =1)
	limits!(ax, -1.7,0.7,-0.5, 1.5)

	scatter!([-1],[1], color=:red)

	xv = [-1,1]
	pv = [2/√5, -1/√5] # steepest direction is -g 
	α = 0.8 
	
	arrows!(ax, [xv[1]], [xv[2]], [α*pv[1]], [α*pv[2]], color=:blue, linewidth=2)
	arrows!(ax, [xv[1]], [xv[2]], [g1[1]], [g1[2]], lengthscale = 0.05,
		color=:black, linewidth=2)
	arrows!(ax, [xv[1]], [xv[2]], [pv[1]], [pv[2]], lengthscale = 0.4,	
		color=:red, linewidth=2)
	
	text!(ax, L"\nabla f", position=(-1.3,1.3), color=:black)
	text!(ax, L"\mathbf{x}", position=(-0.98,1.05), color=:red)
	text!(ax, L"\mathbf{x} + \alpha \mathbf{p}", position=(-0.3,0.73), color=:black)
	text!(ax, L"\mathbf{p}", position=(-0.6,0.89), color=:red)
	
	fig[1,1] = ax
	
	cd("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/")
	
	save("gradient4b.pdf", fig)
	run(`pdfcrop  --clip gradient4b.pdf gradient4.pdf`)
	fig
end

# ╔═╡ b8778b21-7d94-4f98-9af1-8887b6574ffc
let
	p = normalize([-1/2, -√3/2])
	
	H = [2 -1; -1 4]
	Λ, X = eigen(H, sortby = x -> -abs(x)) # Sort eigenvalues from max to min 
										   # like [ 3-√2, 3+√2]
	X[:,1] = normalize(X[:,1]/X[2,1])      # Normalized eigenvector and fix x2 to 1
	X[:,2] = normalize(X[:,2]/X[2,2])	   # Normalized eigenvector and fix x2 to 1
	
	x1s = LinRange(-2,2, 50)
	x2s = LinRange(-2,2, 50)
	z = [f1(x,y) for x in x1s, y in x2s]

	fig = Figure(size = (600,400), font = "CMU Serif")
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = DataAspect())
	contour!(ax,x1s,x2s,z, levels= 0:0.5:10, color=:blue, linewidth =1)
	limits!(ax, -1.5,1.5,-1.5, 1.5)

	scatter!([0],[0], color=:red)

	arrows!(ax, [0], [0], [Λ[1]*X[1,1]], [Λ[1]*X[2,1]], lengthscale=0.2,
			color=:red, linewidth=2)
	arrows!(ax, [0], [0], [Λ[2]*X[1,2]], [Λ[2]*X[2,2]], lengthscale=0.2,
			color=:red, linewidth=2)
	arrows!(ax, [0], [0], [p[1]], [p[2]], lengthscale=1,
			color=:black, linewidth=2)

	text!(ax, L"\alpha_1 \hat{v}_1", position=(-0.5,0.85), color=:red)
	text!(ax, L"p", position=(-0.5,-1.05), color=:black)
	text!(ax, L"\alpha_2 \hat{v}_2", position=(0.3,0.1), color=:red)

	fig[1,1] = ax
	fig
end

# ╔═╡ 809073f1-5118-413e-aa04-0e59df4e2c7a
let
	H = [2 -1; -1 4]
		Λ = eigvals(H) #[3+√2, 3-√2]
		X = eigvecs(H)
		Xn = X/X[2,1]
		#X = [1-√2 1+√2; 1 1]
	V, E = eigen(H, sortby = x -> -abs(x))
	E[:,1] = E[:,1]/E[2,1]
	E[:,2] = E[:,2]/E[2,2]
	V, E
end

# ╔═╡ b8d17f4c-d746-4581-81aa-9471f7d94992
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

# ╔═╡ 8af8efc1-e2f3-4339-995e-1118526389f4
function line_search(f, x, d)
	objective = α -> f(x + α*d)
	a, b = bracket_minimum(objective)
	
	# using brent method from Optim.jl
	res = Optim.optimize(objective, a, b)     # minimize
	α = Optim.minimizer(res)
	return α
end

# ╔═╡ 6cb631d7-6f3f-462c-973b-65cd26351f7e
let
	@variables x₁ x₂ x₃ α₁

	f(x) = sin(x[1]*x[2]) + exp(x[2] + x[3]) - x[3]
	x = [1, 2, 3]
	d = [0, -1, -1]
	α = line_search(f, x, d)

	f1(α₁) = f(x + α₁*d)
 	α1 = 0:0.01:10
	fig1 = Figure(size = (600,400), font = "CMU Serif")
	ax1 = Axis(fig1[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = AxisAspect(1.5))
	lines!(ax1, α1, f1, linewidth=2)
	limits!(ax1, 0, 5, -10, 200)
	text!(ax1, 3.05,8, text=L"\alpha^\ast")
	scatter!(ax1, [3.127], [0], color=(:red, 1.0))
	
	fig1
end

# ╔═╡ 7742066a-e7a3-4582-9bd1-383b0a058034
TableOfContents(title="Local Descent", depth=4)

# ╔═╡ 01bb1467-1e97-4067-84f5-5c4147bd7833
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
# ╠═933d22e8-f7b1-11ec-0597-b3d1bd741c42
# ╠═541e270e-8d76-4c80-a4cf-c40fd846ec5f
# ╠═1868c80b-985d-401b-9f6e-25118474eba7
# ╠═c96c543e-51e7-4e42-81b7-eb251ebc72c7
# ╟─0e3e0db4-de9c-4e2a-b989-7e3dfed2f839
# ╠═d6af1621-ffb8-4a72-86bb-7faa0b2292df
# ╟─343017fd-4f2c-40bc-beea-9ee054832fe6
# ╠═1d69b7a5-b7c1-456c-b6b4-9df161e61458
# ╟─d9c830d6-ec91-4722-a554-6458cfdaa7ff
# ╠═4c9f0680-e805-4ecb-94b4-3ebdc3c76ad3
# ╟─22e69033-c1a2-4408-b4bf-e2117017d094
# ╠═eabce287-97b4-4232-8745-3c0c2f81ea59
# ╠═a0ccee19-b92c-4bd9-a800-c313a4161726
# ╠═b8778b21-7d94-4f98-9af1-8887b6574ffc
# ╠═809073f1-5118-413e-aa04-0e59df4e2c7a
# ╠═b8d17f4c-d746-4581-81aa-9471f7d94992
# ╠═8af8efc1-e2f3-4339-995e-1118526389f4
# ╠═6cb631d7-6f3f-462c-973b-65cd26351f7e
# ╠═7742066a-e7a3-4582-9bd1-383b0a058034
# ╠═01bb1467-1e97-4067-84f5-5c4147bd7833
