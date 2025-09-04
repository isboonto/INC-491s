### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ fb4c1dbe-86f6-11f0-1de8-73cc39eb0e1c
using Pkg; Pkg.activate()

# ╔═╡ 3bc3d06a-eb9a-4b05-831b-2edd31f3fe19
begin
	using Random; Random.seed!(42)
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using LinearAlgebra
	using Statistics # Added for the mean function
	using Latexify, LaTeXStrings
end

# ╔═╡ e87b4097-83c0-4a0e-9cb6-9f998076481c
begin
	
	Random.seed!(42)  # for reproducibility
	m = 100           # number of instances

	X = vec(2 .* rand(m, 1))              # column vector of shape (100,)
	y = vec(4 .+ 3 .* X .+ randn(m, 1)); nothing  # add noise, keep shape consistent
end

# ╔═╡ ff3d6b6a-9718-4d2b-8231-e0106fafd381
begin
	function create_poly_features(x::AbstractVector; degree::Int=10)
		# return matrix with columns [x^0, x^1, ..., x^degree]
		return hcat([x.^i for i in 0:degree]...)  # we need ... !!!!!
	end

	X_b = create_poly_features(X; degree= 1)
	θ_best = (X_b'*X_b) \ X_b'*y        	# Normal Equation
	y_best = X_b * θ_best

	# add x0 = 1 to each instance 
	X_new_b = create_poly_features(X; degree=1)
	y_predict = X_new_b * θ_best
end

# ╔═╡ 3d934c3b-c157-4449-934a-2474be9a9f3f
θ_best     # Normal Equation

# ╔═╡ 29f88011-020d-472e-817a-d5426b1c19a7
let
	global fig1 = Figure(size=(600, 400))
	ax = CairoMakie.Axis(fig1[1,1], xlabel=L"$x$", ylabel=L"$y$")
	#limits!(ax, 0, 2, 0, 15)
	
	scatter!(ax, X,  y, color=:blue,
			 strokecolor=:black, strokewidth=1, label=L"$$ Real Data")
	lines!(ax, X, y_predict, color=:red, linewidth=2, label=L"$$ Predictions")
	#lines!(ax, y_fit, color=:red)

	axislegend(ax, labelsize=14, position=:rb)
	fig1
end

# ╔═╡ 667b1b9b-b276-4f90-8371-c57c6722bca8
md"""
## Batch Gradient Descent
"""

# ╔═╡ 69007dbc-9ade-442a-a7e5-87ed2df349ff
begin
	fig2 = Figure(size=(1000, 300)); nothing
end

# ╔═╡ 60b4c3b7-147c-48dc-9569-1701479b8e83
begin
	function learning_schedule(t)
		t0, t1 = 5, 50 					# learning schedule hyperparameters
		
		return t0 /(t + t1)
	end
	
	function stochastic_gradient(A, y; η = 0.1, n_epochs = 1000)
		N = size(A, 1)
		θ_history = []

		Random.seed!(42)				# to make same noise
		θ = vec(randn(size(A, 2), 1))	# randomly initialized model parameter
		push!(θ_history, θ)

		for epoch in 1:n_epochs
			indices = randperm(N) 		# Shuffle data
			for iteration in indices
				A_i = A[iteration]
				if A_i isa AbstractVector 
					A_i = A_i' 				# Change column to row vector
				end
				y_i = y[iteration] 					
				∇x = (2) * A_i' * (A_i * θ .- y_i)     	# note N = 1
				η = learning_schedule(epoch * N + iteration)
				θ -= η * ∇x
				push!(θ_history, θ)
			end
		end
			
		return θ, θ_history		
	end
end

# ╔═╡ b302d84b-48bf-445d-9972-0fb16fe4457f
θ4 , θ4_history = stochastic_gradient(X_b, y; η = 0.1)

# ╔═╡ 91cf50ec-e959-4b3a-be9d-4d0efac5c525
begin	
	
	function batch_gradient(A, y; batch_size = 50,  η = 0.01, n_epochs = 1000)
		N = size(A, 1)
		θ = zeros(size(A, 2), 1)
		θ_history = []
		
		
		Random.seed!(42)				# to make same noise
		θ = vec(randn(size(A, 2), 1)) 	# randomly initialized model parameters
		push!(θ_history, θ)
		
		for epoch in 1:n_epochs
			indices = randperm(N)
			for start in 1:batch_size:N 
				last = min(start + batch_size - 1, N) # avoid exceed array
				batch_idx = indices[start:last]
				A_batch = A[batch_idx, :]
				if A_batch isa AbstractVector 
					A_batch = A_batch'
				end
				y_batch = y[batch_idx]
				
				∇x = (2/length(A_batch)) * A_batch' * (A_batch * θ - y_batch)
				θ -= η*∇x
				push!(θ_history, θ)
			end
		end

		return θ, θ_history
	end
end

# ╔═╡ b7ba7b27-3319-49cc-b5fc-7989f4f519be
function plot_prediction(fig, ax,  y_predict)
	
	lines!(ax, X, y_predict, color=:red, linewidth=1, label=L"$$ Predictions")

end

# ╔═╡ a21b6619-b9e0-4178-a90f-de2c45bff86d
let
	empty!(fig2)

	η1 = 0.02; η2 = 0.1; η3 = 0.5
	b_size1 = 1; b_size2 = 50;
	global θ1, θ1_history = batch_gradient(X_b, y; batch_size = b_size1, η = η1)	# same like normal equation result
	# full batch
	global θ2, θ2_history = batch_gradient(X_b, y; batch_size = size(X_b,1), η = η1)
	global θ3, θ3_history = batch_gradient(X_b, y; batch_size = b_size2, η = η3)

	ax1 = CairoMakie.Axis(fig2[1,1], xlabel=L"$x$", ylabel=L"$y$", 
									 title = "η = $(η1), batch size = $(b_size2)")
	limits!(ax1, 0, 2, 0, 15)
	scatter!(ax1, X,  y, color=:blue, strokecolor=:black, strokewidth=1, 
				 label="Real Data")
	
	ax2 = CairoMakie.Axis(fig2[1,2], xlabel=L"$x$", ylabel=L"$y$", 
									 title = "η = $(η2), batch size = $(b_size2)")
	limits!(ax2, 0, 2, 0, 15)
	scatter!(ax2, X,  y, color=:blue, strokecolor=:black, strokewidth=1, 
				 label="Real Data")

	ax3 = CairoMakie.Axis(fig2[1,3], xlabel=L"$x$", ylabel=L"$y$", 
									 title = "η = $(η3), batch size = $(b_size2)")
	limits!(ax3, 0, 2, 0, 15)
	scatter!(ax3, X,  y, color=:blue, strokecolor=:black, strokewidth=1, 
				 label="Real Data")
	
	for ii in 1:20   			# first 20 epochs
		y_predict = X_new_b * θ1_history[ii]
		plot_prediction(fig2, ax1, y_predict)

		y_predict = X_new_b * θ2_history[ii]
		plot_prediction(fig2, ax2, y_predict)

		y_predict = X_new_b * θ3_history[ii]
		plot_prediction(fig2, ax3, y_predict)
	end

	fig2
end

# ╔═╡ c357bd46-85f5-4aad-aa56-238744af8780
begin
	M1 = stack(θ1_history, dims=1)
	M2 = stack(θ2_history, dims=1)
	M3 = stack(θ3_history, dims=1)
	M4 = stack(θ4_history, dims=1)

	fig3 = Figure(size=(600, 400))
	ax4 = CairoMakie.Axis(fig3[1,1], xlabel=L"$\theta_0$", ylabel=L"$\theta_1$")
	#limits!(ax4, 2.5, 3.55, 2.3, 4.4)
	
	scatterlines!(ax4, M2[:,1], M2[:,2], label=L"$$ Full Batch", markersize=15,
				 strokewidth=1, strokecolor=:black)
	scatterlines!(ax4, M1[:,1], M1[:,2], label=L"$$ Mini Batch", markersize=15,
				 strokewidth=1, strokecolor=:black)
	scatterlines!(ax4, M4[:,1], M4[:,2], label=L"$$ SGD", markersize=15, 
				 strokewidth=1, strokecolor=:black)
	
	axislegend(ax4, labelsize=14, position=:rb)

	fig3
end

# ╔═╡ 3ac6461d-e8c4-4387-870d-0f2bb8b25ac2
begin
	θ4_history[end]
	θ3_history[end]
end

# ╔═╡ 01003261-142b-46fc-b6df-4ecdd4c8b01c
begin
		fig4 = Figure(size=(600, 400))	
		ax5 = CairoMakie.Axis(fig4[1,1], xlabel=L"$x$", ylabel=L"$y$")
		limits!(ax5, 0, 2, 0, 15)
		scatter!(ax5, X,  y, color=:blue, strokecolor=:black, strokewidth=1, 
					 label="Real Data")
		for ii = 1:20
			y_predict = X_new_b * θ4_history[ii]
			plot_prediction(fig4, ax5, y_predict)
		end

		fig4
end

# ╔═╡ Cell order:
# ╠═fb4c1dbe-86f6-11f0-1de8-73cc39eb0e1c
# ╠═3bc3d06a-eb9a-4b05-831b-2edd31f3fe19
# ╠═e87b4097-83c0-4a0e-9cb6-9f998076481c
# ╠═ff3d6b6a-9718-4d2b-8231-e0106fafd381
# ╠═3d934c3b-c157-4449-934a-2474be9a9f3f
# ╠═29f88011-020d-472e-817a-d5426b1c19a7
# ╟─667b1b9b-b276-4f90-8371-c57c6722bca8
# ╠═69007dbc-9ade-442a-a7e5-87ed2df349ff
# ╠═a21b6619-b9e0-4178-a90f-de2c45bff86d
# ╠═60b4c3b7-147c-48dc-9569-1701479b8e83
# ╠═b302d84b-48bf-445d-9972-0fb16fe4457f
# ╠═91cf50ec-e959-4b3a-be9d-4d0efac5c525
# ╠═b7ba7b27-3319-49cc-b5fc-7989f4f519be
# ╠═c357bd46-85f5-4aad-aa56-238744af8780
# ╠═3ac6461d-e8c4-4387-870d-0f2bb8b25ac2
# ╠═01003261-142b-46fc-b6df-4ecdd4c8b01c
