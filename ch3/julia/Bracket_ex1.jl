### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ 3c294290-2691-11ed-08e0-bf8c43629af8
begin
	using Pkg; Pkg.activate()
	using CairoMakie
	set_theme!(theme_latexfonts())
	using Symbolics, ForwardDiff, LinearAlgebra
	using PlutoUI, LaTeXStrings
	using Images
end

# ╔═╡ 224f73a4-697c-4a8f-afc6-53b447150856
md"""
## Example 1

Consider the problem:

$\operatorname*{minimize}_x \qquad x^2 + \frac{54}{x}$

in the interval $(0,5)$. **Note:** The minimum point lies at $x^\ast = 3$.
"""

# ╔═╡ 7d20ac53-ad7f-49e0-8fd7-d1dabec8b8e3
begin
	f1 = x -> x^2 + 54/x
end

# ╔═╡ bbba2972-b3b1-4b65-90ba-b6f7faf451c0
function bracket_minimum(f, x=0; Δ=1e-2, γ=2.0)
	a, ya = x, f(x)
	b, yb = a + Δ, f(a + Δ)
	if yb > ya
		a, b = b, a
		ya, yb = yb, ya
		Δ = -Δ
	end
	Δ *= γ
	while true
		c, yc = b + Δ, f(b + Δ)
		if yc > yb
			println("$a  $c")
			return a < c ? (a, c) : (c, a)
		end
		a, ya, b, yb = b, yb, c, yc
		Δ *= γ
		println("$a  $c")
	end
end

# ╔═╡ 23531ff8-f3fe-4343-9a4a-e5b915b6e927
begin
	x1s = 0:0.1:11; 
	
	fig1 = Figure(size = (600,400))
	ax = CairoMakie.Axis(fig1[1,1], xlabel = L"x", ylabel = L"f(x)", 
		title=L"f(x)=x^2 + 54/x", aspect = AxisAspect(1.5), 
		xticks = 0:11, yticks = 0:25:200, 
		xlabelsize=20, ylabelsize=20, titlesize=20)
	lines!(ax,x1s,f1, color=(:blue, 1.0), linewidth =2)
	limits!(ax, 0,11,0, 200) 

	fig1
end

# ╔═╡ d2c7c947-0efc-422f-b663-bbb2573bd5ba
bound = bracket_minimum(f1, 0, Δ= 1e-2, γ=1.2)

# ╔═╡ e0122e2d-0074-4693-92b4-f05df3cf25b1
md"""
This bound $(bound) guarantees that the minimum point lies in the bracket. 
"""

# ╔═╡ de13d815-2088-433b-94dd-37c32a336954
function fibonacci_search(f, a, b, n; ϵ=0.01)
	φ = 1.61803
	s = (1-√5)/(1+√5)
	ρ = 1/ (φ*(1-s^(n+1))/(1-s^n))
	d = ρ*b + (1-ρ)*a
	yd = f(d)
	for i in 1:n-1
		if i == n-1
			c = ϵ*a + (1-ϵ)*d
		else
			c = ρ*a + (1-ρ)*b
		end
		yc = f(c)
		if yc < yd
			b, d, yd = d, c, yc
		else
			a, b = b, c
		end
		ρ = 1 / (φ*(1-s^(n-i+1))/(1-s^(n-i)))
	end
	return a < b ? (a, b) : (b, a)
end

# ╔═╡ 24a25e9b-0aef-4b42-a8d5-84347054dc06
let
	bound = fibonacci_search(f1, 0, 4.68, 19)
	#bound[2] - bound[1]
end

# ╔═╡ c42d8c9e-eafe-49b6-bd50-d0f97fef68ba
md"""
## Example: Projectile release
"""

# ╔═╡ a75c1f5a-c3ad-4f8f-bcb3-f874a56fbe99
let
	url = "Fibonacci3.png"
	house = load(url)
end

# ╔═╡ b2460d35-4e36-4072-b988-c3bc7ef71134
let
	h = 50; V = 90; g = 9.81 
	F = θ -> -((V*sin(θ*π/180)/g) + √((2h/g) + (V*sin(θ*π/180)/g)^2))*V*cos(θ*π/180)
	
	# 7 interval reductions -> n = 8
	# 19 interval reductions -> n = 20
	
	bound7 = fibonacci_search(F,0,80,8)
	println("Final interval for 7 reductions", bound7)
	bound19 = fibonacci_search(F,0,80,20)
	println("Final interval for 19  reductions", bound19)
	md"""
	7 steps $bound7 \

	19 steps $bound19
	"""
end

# ╔═╡ 813ef5cb-ca3e-4700-bfe7-b6504672c92a
md"""
## Example: Five function evaluations
"""

# ╔═╡ df435d58-33c8-49a7-924c-a7d217759ea4
md"""
Consider using Fibonacci search with five function evaluations to minimize $f(x) = e^{x-2} - x$ over the interval $[a, b] = [-2, 6]$. The first two function evaluations are made at $\frac{F_5}{F_6}$ and $1-\frac{F_5}{F_6}$, along the length of the initial bracketing interval:
"""

# ╔═╡ 313becbc-0ba1-4d09-b6ae-6c85055f5244
function golden_section_search(f, a, b, n)
	φ = 1.61803
	ρ = φ - 1
	d = ρ * b + (1 - ρ)*a
	yd = f(d)
	for i = 1:n-1
		c = ρ*a + (1 - ρ)*b
		yc = f(c)
		if yc < yd
			b, d, yd = d, c, yc
		else
			a, b = b, c
		end
		println( a < b ? (a, b) : (b, a))  # to test the nonunimodal function
	end
	return a < b ? (a, b) : (b, a)
end

# ╔═╡ 6cd732d9-7773-45cb-b9fa-d344020c27d4
let
	a1 = -2.0; b1 = 6.0; F56 = 8.0/13.0;
	ft = x -> exp(x-2) - x
	fx1 = ft(round(a1 + (b1 - a1)*(1 - F56)))
	fx2 = ft(round(a1 + (b1 - a1)*F56))
	bound4 = fibonacci_search(ft, -2,6, 15)
	boundg4 = golden_section_search(ft,-2,6,15)
	md"""
	Fibonacci = $bound4\
	Golden Ratio = $boundg4
	"""
end

# ╔═╡ f9fd736b-b6dd-403a-b06d-d4aa086db4c1
let
	f = x-> (sin(x) + sin(x/2))/4
	a₀, b₀ = -5.0, 9.0

	bound5 = golden_section_search(f, a₀, b₀, 6)
end

# ╔═╡ 129dfb6d-6da8-416d-953c-4e30921fb777
md"""
## Quadratic Fit Search
"""

# ╔═╡ 35bc8c74-a9a4-4d2a-a898-eee99a1996d7
function quadratic_fit_search(f, a, b, c, n)
	ya, yb, yc = f(a), f(b), f(c)
	for i in 1:n-3 # uncomment if you don't want to show plot.
		x = 0.5*(ya*(b^2 - c^2) + yb*(c^2 - a^2) + yc*(a^2 - b^2)) / 
			(ya*(b - c) + yb*(c -a) + yc*(a - b))
		yx = f(x)
		if x > b
			if yx > yb
				c , yc = x, yx
			else
				a, ya, b, yb = b, yb, x, yx
			end
		elseif x < b
			if yx > yb
				a, ya = x, yx
			else
				c, yc, b, yb = b, yb, x, yx
			end
		end
	end
	return (a, b, c)
end

# ╔═╡ bdc67781-3950-4c21-902e-67906c753d06
begin
	ftest = x->x - 2x^2 + 0.2x^3 
	
	new_p1 = quadratic_fit_search(ftest, 0.5, 4.0, 9.75, 4)
end

# ╔═╡ cde39806-03c9-4593-bb33-73e427e2dc5b
function quadratic_fit_search1(f, a, b, c, n)
	ya, yb, yc = f(a), f(b), f(c)
	#for i in 1:n-3 # uncomment if you don't want to show plot.
		x = 0.5*(ya*(b^2 - c^2) + yb*(c^2 - a^2) + yc*(a^2 - b^2)) / 
			(ya*(b - c) + yb*(c -a) + yc*(a - b))
		yx = f(x)
		if x > b
			if yx > yb
				c , yc = x, yx
			else
				a, ya, b, yb = b, yb, x, yx
			end
		elseif x < b
			if yx > yb
				a, ya = x, yx
			else
				c, yc, b, yb = b, yb, x, yx
			end
		end
	#end
	return (a, b, c)
end

# ╔═╡ 181a153a-8036-4062-b7b1-dea4b691d286
begin
	function get_q(a,b,c)
	    fa = f(a)
	    fb = f(b)
	    fc = f(c)
	    return x-> fa*(x-b)*(x-c)/((a-b)*(a-c)) +
	               fb*(x-a)*(x-c)/((b-a)*(b-c)) +
	               fc*(x-a)*(x-b)/((c-a)*(c-b))
	end

	
	f = x->x - 2x^2 + 0.2x^3 
	
	global a = 0.5 ; global b = 4.0 ; global c = 9.75
	global ya = f(a); global yb = f(b); global yc = f(c)
	
	 
	x1s1 = 0:0.1:11; 
	
	fig2 = Figure(size = (600,600))
	ax3 = CairoMakie.Axis(fig2[1,1], xlabel = L"x", ylabel = L"f(x)", title=L"$ $ Interation 1", aspect = AxisAspect(1.2), xticks = 0:3:11, yticks = -20:20:200, xlabelsize=20, ylabelsize=20)
	limits!(ax3,0,11,-25, 50)
	ax4 = CairoMakie.Axis(fig2[1,2], xlabel = L"x", ylabel = L"f(x)", title=L"$ $ Interation 2", aspect = AxisAspect(1.2), xticks = 0:3:11, yticks = -20:20:200,
	xlabelsize=20, ylabelsize=20)
	limits!(ax4,0,11,-25, 50)
	ax5 = CairoMakie.Axis(fig2[2,1], xlabel = L"x", ylabel = L"f(x)", title=L"$ $ Interation 3", aspect = AxisAspect(1.2), xticks = 0:3:11, yticks = -20:20:200,
	xlabelsize=20, ylabelsize=20)
	limits!(ax5,0,11,-25, 50)
	ax6 = CairoMakie.Axis(fig2[2,2], xlabel = L"x", ylabel = L"f(x)", title=L"$ $ Interation 4", aspect = AxisAspect(1.2), xticks = 0:3:11, yticks = -20:20:200,
	xlabelsize=20, ylabelsize=20)
	limits!(ax6,0,11,-25, 50)
	axx = [ax3, ax4, ax5, ax6]
	   

	for i in 1:4
		
		x = 0.5*(ya*(b^2-c^2)+yb*(c^2-a^2)+yc*(a^2-b^2)) /
	            (ya*(b-c)    +yb*(c-a)    +yc*(a-b))
		yx = f(x)
		scatter!(axx[i], [x], [yx], color = (:red, 1.0))
		scatter!(axx[i], [a, b, c], [ya, yb, yc], color=(:black, 1.0) )
		lines!(axx[i],x1s1,f, color=(:blue, 1.0), linewidth =2)
		lines!(axx[i], x1s1, get_q(a,b,c), color=(:red, 1.0), linewidth=2)
		
		new_p = quadratic_fit_search1(f, a, b, c, 1)
		global a = new_p[1]; global b = new_p[2]; global c = new_p[3]
		global ya = f(a); global yb = f(b); global yc = f(c)
		
	end
	
	#empty!(fig2)
	fig2
end

# ╔═╡ fb629c03-c58e-44e8-9a8c-dd7b4b98ea22
function bisection(f, a, b, ϵ)
	if a > b; a, b = b, a; end # ensure a < b

	ya, yb = f(a), f(b)
	if ya == 0; b = a; end
	if yb == 0; a = b; end

	while b - a > ϵ
		x = (a + b)/2
		y = f(x)
		if y == 0 
			a, b = x, x 
		elseif sign(y) == sign(ya)
			a = x
		else
			b = x
		end
	end

	return (a, b)
end

# ╔═╡ cce0e06c-02fb-4f10-882c-aee0529db953
f2 = x-> (sin(x) + sin(x/2))/4

# ╔═╡ 505eff72-a003-49f2-b9a5-cada1c6ddaf1
a7, b7 = bisection(f2, -1, 20, 1e-4)

# ╔═╡ 9e49fecd-6ad0-4e54-a2e3-ca3732519acd
TableOfContents(title="Local Descent", depth=4)

# ╔═╡ c344e5fc-3a24-4089-b3b6-d22843dd68c9
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
# ╠═3c294290-2691-11ed-08e0-bf8c43629af8
# ╟─224f73a4-697c-4a8f-afc6-53b447150856
# ╠═7d20ac53-ad7f-49e0-8fd7-d1dabec8b8e3
# ╠═bbba2972-b3b1-4b65-90ba-b6f7faf451c0
# ╠═23531ff8-f3fe-4343-9a4a-e5b915b6e927
# ╠═d2c7c947-0efc-422f-b663-bbb2573bd5ba
# ╠═e0122e2d-0074-4693-92b4-f05df3cf25b1
# ╠═de13d815-2088-433b-94dd-37c32a336954
# ╠═24a25e9b-0aef-4b42-a8d5-84347054dc06
# ╟─c42d8c9e-eafe-49b6-bd50-d0f97fef68ba
# ╟─a75c1f5a-c3ad-4f8f-bcb3-f874a56fbe99
# ╠═b2460d35-4e36-4072-b988-c3bc7ef71134
# ╠═813ef5cb-ca3e-4700-bfe7-b6504672c92a
# ╟─df435d58-33c8-49a7-924c-a7d217759ea4
# ╠═6cd732d9-7773-45cb-b9fa-d344020c27d4
# ╠═313becbc-0ba1-4d09-b6ae-6c85055f5244
# ╠═f9fd736b-b6dd-403a-b06d-d4aa086db4c1
# ╟─129dfb6d-6da8-416d-953c-4e30921fb777
# ╠═35bc8c74-a9a4-4d2a-a898-eee99a1996d7
# ╠═bdc67781-3950-4c21-902e-67906c753d06
# ╠═cde39806-03c9-4593-bb33-73e427e2dc5b
# ╠═181a153a-8036-4062-b7b1-dea4b691d286
# ╠═fb629c03-c58e-44e8-9a8c-dd7b4b98ea22
# ╠═cce0e06c-02fb-4f10-882c-aee0529db953
# ╠═505eff72-a003-49f2-b9a5-cada1c6ddaf1
# ╠═9e49fecd-6ad0-4e54-a2e3-ca3732519acd
# ╠═c344e5fc-3a24-4089-b3b6-d22843dd68c9
