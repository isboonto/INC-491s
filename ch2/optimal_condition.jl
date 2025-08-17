### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ 0b06194d-d6aa-4e95-99a8-c1c02c703e37
using Pkg; Pkg.activate()

# ╔═╡ 58e070c0-1e1a-11ed-37cc-811aa6ef98f5
begin
	using PlutoUI, StaticArrays
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using LinearAlgebra, ForwardDiff
end

# ╔═╡ afaec5e7-5f1f-47e8-8660-5ccf096582f9
using Symbolics, Roots

# ╔═╡ c4912b42-b498-4e6e-9b5a-489c484f08c0
md"""
## Example 1
"""

# ╔═╡ 7e903893-4bbe-453e-8645-b0be042e4ebe
md"""
Consider a nonlinear function

$f(x_1,x_2) = 0.5x_1^4 + 2x_1^3 + 1.5x_1^2 + x_2^2 - 2x_1x_2$

"""

# ╔═╡ cdaf756d-7207-4c15-831b-70f0c7297157
begin
	f(x) = 0.5x[1]^4 + 2x[1]^3 + 1.5x[1]^2 + x[2]^2 - 2x[1]*x[2]
	∇f(x) = ForwardDiff.gradient(f,x)
	Hf(x) = ForwardDiff.hessian(f,x)
end

# ╔═╡ 820dd53c-3bf8-44b2-b890-c145d196bd19
md"""
 Let the gradient equal to zero.
"""

# ╔═╡ d0b4e7ee-5858-4687-abc5-4febb602053b
begin
	@variables x₁ x₂
	
	println("f(x) = $(f([x₁ x₂]))")
	println("∇f(x) = $(∇f([x₁ x₂]))")
	println("H(x) = $(Hf([x₁ x₂]))")
	
	# At the optimal ∇f(x) = 0
	# From the second row -2x1 + 2x2 = 0, we have x1 = x2
	∇f1(x₁) = substitute(∇f([x₁,x₂])[1],Dict(x₂=>x₁))
	result = find_zeros(∇f1, [-10, 1])
    println(result) 
	
	xA = [result[3], result[3]]
	xC = [result[2], result[2]]
	xB = [result[1], result[1]]

	HfxA = Hf([xA[1], 0])
	HfxB = Hf([xB[1], 0])
	HfxC = Hf([xC[1], 0])
end

# ╔═╡ c4397c5a-a869-412e-9055-95f19ebe7ad3
begin
	eigvals(HfxA)	# Positive definite, global minimum
	eigvals(HfxB)	# in-definite saddle point
	eigvals(HfxC)	# Positive definite, local minima
end

# ╔═╡ 48dcc503-210d-4f55-83e6-340aacb41a0e
begin
	x1s = LinRange(-4,2, 50)
	x2s = LinRange(-4,2, 50)
	z = [f([x,y]) for x in x1s, y in x2s]

	fig = Figure(size = (600,400))
	ax = Axis(fig[1,1], xlabel = L"x_1", ylabel = L"x_2", aspect = AxisAspect(1.5))
	contour!(ax,x1s,x2s,z, levels= -8:2:-2, linewidth = 2)
	contour!(ax,x1s,x2s,z, levels= 0.02:2:60, linewidth =2)
	limits!(ax, -4,2,-4, 2)

	scatter!(ax,result[1], result[1], color=(:red, 1.0), markersize=12)
	scatter!(ax,result[2:3], result[2:3], color=(:black, 1.0), markersize=12)
	text!(ax, result[1], result[1] + 0.25, text=L"$$global minimum", fontsize=14)
	text!(ax, result[2], result[2] - 0.5, text=L"$$saddle point", fontsize=14)
	text!(ax, result[3], result[3] + 0.25, text=L"$$local minimum", fontsize=14)

	fig
end

# ╔═╡ f62f5465-4a9b-42bd-bbc4-ba574e19cf10
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/condition_hessian_ex1.pdf", fig)

# ╔═╡ b0072bc0-0c20-49ec-8c62-351058fd4889
save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/condition_hessian_ex1.png", fig)

# ╔═╡ 00157155-b5eb-4b6c-a074-e225fc132374
TableOfContents(title="Local Descent", depth=4)

# ╔═╡ 72a9af8b-271a-4b22-9f48-cc0e57bce38e
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
# ╠═0b06194d-d6aa-4e95-99a8-c1c02c703e37
# ╠═58e070c0-1e1a-11ed-37cc-811aa6ef98f5
# ╠═afaec5e7-5f1f-47e8-8660-5ccf096582f9
# ╟─c4912b42-b498-4e6e-9b5a-489c484f08c0
# ╟─7e903893-4bbe-453e-8645-b0be042e4ebe
# ╠═cdaf756d-7207-4c15-831b-70f0c7297157
# ╟─820dd53c-3bf8-44b2-b890-c145d196bd19
# ╠═d0b4e7ee-5858-4687-abc5-4febb602053b
# ╠═c4397c5a-a869-412e-9055-95f19ebe7ad3
# ╠═48dcc503-210d-4f55-83e6-340aacb41a0e
# ╠═f62f5465-4a9b-42bd-bbc4-ba574e19cf10
# ╠═b0072bc0-0c20-49ec-8c62-351058fd4889
# ╠═00157155-b5eb-4b6c-a074-e225fc132374
# ╠═72a9af8b-271a-4b22-9f48-cc0e57bce38e
