### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ dbba4408-76c3-11f0-265d-cda3574f5dc7
using Pkg; Pkg.activate()

# ╔═╡ 80619445-119e-40eb-997e-71d77031fe7f
using Symbolics, ForwardDiff, LinearAlgebra

# ╔═╡ 3f890fb2-f1f9-4bda-9083-7aabd875b622
begin
	@variables x1 x2 x3 

	f(x) = 3x[1]^2 + 2x[1]x[2] + x[1]x[3] + 2.5x[2]^2 + 2x[2]x[3] + 2x[3]^2 - 8x[1] - 3x[2] - 3x[3]

	gf = (f, x) -> ForwardDiff.gradient(f, x)
	Hf = (f, x) -> ForwardDiff.hessian(f, x)

	x = [x1, x2, x3]
	f1 = f(x)
	g = gf(f, x)
	H = Hf(f, x)
	
	println("Function f(x) = $(f1)")
	println("Gradinet of f(x) = $(g)")
	println("Hessian of f(x) = $(H)")
end

# ╔═╡ 13f4cacd-f6e0-4d34-a71d-671855453896
begin
	f2(x) = 2x[1] + x[2]/x[1]

	
	x0 = [1, 0.5]
	flinear(x) = f2(x0) + gf(f2, x0)'*(x - x0)
	fquadratic(x) = f2(x0) + gf(f2, x0)'*(x - x0) + (1/2)*(x - x0)'*Hf(f2, x0)*(x - x0)

	println(expand(flinear([x1, x2])))
	println(expand(fquadratic([x1, x2])))
end

# ╔═╡ Cell order:
# ╠═dbba4408-76c3-11f0-265d-cda3574f5dc7
# ╠═80619445-119e-40eb-997e-71d77031fe7f
# ╠═3f890fb2-f1f9-4bda-9083-7aabd875b622
# ╠═13f4cacd-f6e0-4d34-a71d-671855453896
