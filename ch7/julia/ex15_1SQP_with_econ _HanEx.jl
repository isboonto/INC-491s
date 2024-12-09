### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 79e23548-a99a-11ef-232f-c7571aff27bb
using Pkg; Pkg.activate()

# ╔═╡ dd44d122-f8c3-49ab-825c-82e0a4c0e8a9
begin
	using LinearAlgebra, ForwardDiff, Symbolics
	using CairoMakie; 
		set_theme!(theme_latexfonts(), fontsize=18)
end

# ╔═╡ 30424885-32cb-49ca-a4bc-4e133218fdf2
begin
	f = x -> -x[1] - x[2]
	a2 = x -> x[1]^2 - x[2]
	a1 = x -> x[1]^2 + x[2]^2 - 1
	ax = x -> [a1(x), a2(x)]
	#ax = x -> [a1(x)]
end

# ╔═╡ b214d950-6d6a-46eb-9796-8676b5e4eac2
begin
	# For a scalar input x
	# For scalar input x
	function ∇f(f,x)
		return ForwardDiff.gradient(f,[x])
	end

	# For a vector input x
	function ∇f(f, x::Array)
		return ForwardDiff.gradient(f,x)
	end

	function ∇f(f::Array, x::Array)
		return ForwardDiff.jacobian(x) do f end
	end
	
	function ∇2f(f,x)
			# to specific that g is a gradient of f
		return ForwardDiff.hessian(f,[x])
	end
	
	# For a vector input x
	function ∇2f(f, x::Array)
		return ForwardDiff.hessian(f,x)
	end
		
end

# ╔═╡ 798666d4-ffe7-4e15-8585-9f9d564554ee
begin
	@variables x1, x2, x3, λ1, λ2
	
	xs0 = [x1, x2]; λs0 = [λ1, λ2]
	
	As = ForwardDiff.jacobian(ax, xs0)
	rd = -∇f(f, xs0) + As'*λs0
	Ws = ∇2f(f, xs0) + ∇2f(a1,xs0)*λs0[1] + ∇2f(a2, xs0)*λs0[2]
end

# ╔═╡ 73e9f289-fb29-4e4d-b3ce-32e227270069
md"""
Solving 

$\begin{align*}
&\text{minimize} \qquad  f(\textbf{x})\\
&\text{subject to} \quad a_i = 0 \quad \text{ where } i = 1,\ldots,p
\end{align*}$

$\begin{align*}
\begin{bmatrix}\mathbf{W}_k & \mathbf{A}_k^T \\ \mathbf{A}_k & 0 \end{bmatrix}\begin{bmatrix}\delta_x \\ \delta_\lambda\end{bmatrix} = \begin{bmatrix}-\mathbf{g}_k - \mathbf{A}_k^T\lambda_k \\ -\mathbf{a}_k\end{bmatrix}
\end{align*}$
If $\lambda_{k+1} = \lambda_k + \delta_k$, we can have

$\begin{align*}
\begin{bmatrix}\mathbf{W}_k & \mathbf{A}_k^T \\ \mathbf{A}_k & 0 \end{bmatrix}\begin{bmatrix}\delta_x \\ \lambda_{k+1}\end{bmatrix} = \begin{bmatrix}-\mathbf{g}_k \\ -\mathbf{a}_k\end{bmatrix}
\end{align*}$

where

$\begin{align*} \mathbf{W}_k &= \nabla_{xx}^2f(\textbf{x}_k) + \sum_{i=1}^p(\mathbf{\lambda}_k)_i\nabla_{xx}^2a_i(\textbf{x}_k)\\
\mathbf{A}_k &= \begin{bmatrix}\nabla_x^Ta_1(\textbf{x}_k) \\ \vdots \\ \nabla_x^Ta_p(\textbf{x}_k)\end{bmatrix}, \qquad \mathbf{g}_k = \nabla_xf(\textbf{x}) \\
\mathbf{a}_k &= \begin{bmatrix}a_1(x_k) & a_2(x_k) & \cdots & a_p(x_k)\end{bmatrix}^T\\
\lambda_k &\geq 0\end{align*}$
The problem can be considered as a quadratic problem:

$\begin{align*}
	&\text{minimize} \quad \frac{1}{2}\delta^T\mathbf{W}_k\delta + \delta^T\mathbf{g}_k \\
	&\text{subject to} \quad \mathbf{A}_k\delta = -\mathbf{a}_k
\end{align*}$

"""

# ╔═╡ 14eabf2d-ee66-43aa-8f38-23b371eb6856
begin
	function Qp(f, ax, x0, λ0::Array)
		# Should be ∇f(c,x) because we need to substitute the value of x after the jacobian
		Ak = ForwardDiff.jacobian(ax, x0)	
		
		rd = -∇f(f, x0) #- Ak'*λ0
		rp = -ax(x0)
		N = size(rp,1)
		
		aa = [ x->ax(x)[i]  for i = 1:N ] # make anonymous function to be array of fn
		Wk = ∇2f(f, x0)
	
		for i = 1:N
			Wk = Wk + ∇2f(aa[i],x0)*λ0[i]
		end
		
		δz = [Wk Ak'; Ak zeros(size(Ak,1), size(Ak',2))]\[rd; rp]
	    δx = δz[1:size(x0,1)]
	    λ = δz[size(x0,1)+1:end]
			  
		return x0 + δx, λ
	end
	
	# Analytical Solution	
	function Qp2(f, ax, x0, λ0::Array)
			gk = ∇f(f, x0)
			Ak = ForwardDiff.jacobian(ax,x0)
			
			rp = ax(x0)
			N = size(rp,1)

			aa = [x-> ax(x)[i] for i = 1:N]
			Wk = ∇2f(f,x0)
			for i = 1:N
				Wk = Wk + ∇2f(aa[i],x0)*λ0[i]
			end
			H = Ak*(Wk\(Ak'))
			λ = inv(H)*(rp - Ak*inv(Wk)*gk)
			dx = -inv(Wk)*(Ak'*λ + gk)
	
			return x0 + dx, λ
	end
end

# ╔═╡ b832474b-7bdd-41c9-a22d-2da760a4389d
begin
	kmax = 100
	x = zeros(2, kmax)
	λ = zeros(2, kmax)
	x0 = x[:,1] = [0.5, 1]; λ0 = λ[:,1] =  [1, 1]; 
	δx = 100*ones(3,1);
	k = 1
	while  (norm(δx) >= 1e-10) && (k <= kmax)
		global xnew, λnew = Qp(f, ax, x[:,k], λ[:,k])
		global δx = x[:,k] - xnew;
		global x[:,k+1] = xnew; 
		global λ[:,k+1] = λnew
		global k = k + 1
	end
	xp = x[:,1:k]
	λp = λ[:,1:k]; nothing
end

# ╔═╡ 0f5850ad-26ed-4625-8319-22ada4d65941
begin
	xp, λp, k, f(xp[:,end]), a1(xp[:,end]), a2(xp[:,end])
	@show xp[:,end]
	@show λp[:,end]
	@show k
	@show f(xp[:,end])
end

# ╔═╡ e851eb8f-0e3f-45de-934e-20138dd497d8
begin
	
	fig2 = Figure(size=(600,400), fontsize=20)
	ax2 = Axis(fig2[1,1], xlabel=L"$x_1$", ylabel=L"$x_2$", aspect=1.4)
	limits!(ax2, -1.5, 1.5, -0.20, 1.5)

	
	xs1 = -4:0.02:4; 
	fc = (x1, x2) -> f([x1, x2])
	a1c = (x1, x2) -> a1([x1, x2])
	a2c = (x1, x2) -> a2([x1, x2])
	contour!(ax2, xs1, xs1, fc, levels=-20:0.5:20)
	contour!(ax2, xs1, xs1, a1c, levels=[0], linewidth=2, color=:red)
	contour!(ax2, xs1, xs1, a2c, levels=[0], linewidth=2, color=:blue)

	xb = -0.8:0.1:0.8
	g1b = (x) -> x^2 
	g2b = (x) -> sqrt(1-x^2)
	band!(ax2, xb, g2b.(xb), g1b.(xb), color=(:orange, 0.4))
	
	scatterlines!(ax2, x[1,:], x[2,:], markersize=10, color=:blue, linestyle=:dash)
	scatter!(ax2, 1/sqrt(2), 1/sqrt(2), markersize=10, color=:red, 
		strokecolor=:black, strokewidth=1)
	text!(ax2, 0.52, 0.6, text=L"\textbf{x}^\ast")
	text!(ax2, x[1,1], x[2,1], text=L"\textbf{x}_0")
	fig2
end

# ╔═╡ 0e0d073d-6a28-4d7d-b46f-7ea7a0840b90
begin
	fig1 = Figure(size=(600,400), fontsize=18)
	ax1 = Axis(fig1[1,1], xlabel=L"$$ Step $k$", ylabel=L"$f(\textbf{x})$");
	N = size(xp,2)
	fp = [f(xp[:,i]) for i in 1:N]
	a1p = [a1(xp[:,i]) for i in 1:N]
	a2p = [a2(xp[:,i]) for i in 1:N]
	scatterlines!(ax1, fp, linewidth=2, markersize=15, strokecolor=:black, 
		strokewidth=1, label=L"objective $f(\textbf{x})$")
	scatterlines!(ax1, a1p, linewidth=2, markersize=15, strokecolor=:black,
		strokewidth=1, label=L"constraint $a_1(\textbf{x})$")
	scatterlines!(ax1, a2p, linewidth=2, markersize=15, strokecolor=:black,
		strokewidth=1, label=L"constraint $a_2(\textbf{x})$")
	lines!(ax1, 0*a2p, linestyle=:dash, linewidth=2, color=:black)
	limits!(ax1, 0, N+1, -3, 2)

	axislegend(ax1)
	fig1
end

# ╔═╡ Cell order:
# ╠═79e23548-a99a-11ef-232f-c7571aff27bb
# ╠═dd44d122-f8c3-49ab-825c-82e0a4c0e8a9
# ╠═30424885-32cb-49ca-a4bc-4e133218fdf2
# ╠═b214d950-6d6a-46eb-9796-8676b5e4eac2
# ╠═798666d4-ffe7-4e15-8585-9f9d564554ee
# ╟─73e9f289-fb29-4e4d-b3ce-32e227270069
# ╠═14eabf2d-ee66-43aa-8f38-23b371eb6856
# ╠═b832474b-7bdd-41c9-a22d-2da760a4389d
# ╠═0f5850ad-26ed-4625-8319-22ada4d65941
# ╠═e851eb8f-0e3f-45de-934e-20138dd497d8
# ╠═0e0d073d-6a28-4d7d-b46f-7ea7a0840b90
