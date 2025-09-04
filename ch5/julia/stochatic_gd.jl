### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ f0f4f39c-8689-11f0-169d-2fe53b79ca12
using Pkg; Pkg.activate()

# ╔═╡ def70522-1384-4325-8f91-9dfb09681b0f
begin
	using LinearAlgebra
	using CairoMakie; set_theme!(theme_latexfonts(), fontsize=18)
	using ForwardDiff, Random, LaTeXStrings, Latexify, Symbolics
end

# ╔═╡ ab0fcc0c-dd93-48e2-a844-28acb853564e
md"""
## SGD on a Linear Regression Problem
"""

# ╔═╡ 73ffb1e3-e328-4925-838d-504e78748cee
begin
	N = 2000 				# number of samples
	A = randn(N, 2)			# design matrix
	true_x = [2.0, -3.0] 	# true slope
	true_b = 5.0

	# True y = Ax + b plus some noise
	noise_std = 0.5
	noise = noise_std * randn(N)
	y = A*true_x .+ true_b .+ noise
end

# ╔═╡ 7460106e-9309-451b-8672-a01a0a9655f1
begin
	function predict(A, x, b)
		# A: shape (batch_size, 2)
		# x: shape (2, 0)
		# b: scalar
			
		return A*x .+ b       # shape (batch_size, )
	end

	

	function mean_squared_error(y_true, y_pred)
		
		return sum((y_true - y_pred).^2)/length(y_true)
	end

	function batch_gradients(A_batch, y_batch, x, b)
		# A_batch, x, b: same shape
		# y_batch: shape (batch_size,)

		batch_size = size(A_batch, 1)
		y_pred = predict(A_batch, x, b)

		# ensure residuals is 1D
		residuals = y_pred  .- y_batch

		∇x = (2.0 / batch_size) * (A_batch' * residuals)
		∇b = (2.0 / batch_size) * sum(residuals)

		return ∇x, ∇b
	end

	function full_batch_gd(A, y; α=0.01, n_epochs=50)
		x = zeros(size(A, 2))  			# unknown parameter vector
		b = 0.0
		mse_history = Float64[]

		for epoch in 1:n_epochs
			∇x, ∇b = batch_gradients(A, y, x, b)
			x = x - α * ∇x
			b = b - α * ∇b 
		
			y_pred = predict(A, x, b)
			mse = mean_squared_error(y, y_pred)
			push!(mse_history, mse)
		end

		return x, b, mse_history
	end

end

# ╔═╡ fd7872af-3085-49ba-bf7f-97fc3761075f
begin
	function mini_batch_sgd(A, y; batch_size = 50, α = 0.01, n_epochs=50)
		x = zeros(size(A, 2)) 
		b = 0.0 
		N = size(A, 1)
		mse_history = Float64[]

		for epoch in 1:n_epochs
			indices = randperm(N)			# Shuffle data
			for start in 1:batch_size:N
				last = min(start + batch_size - 1, N)  # avoid exceed array size
				batch_idx = indices[start:last]
				A_batch = A[batch_idx, :]
				if A_batch isa AbstractVector 
					A_batch = A_batch'
				end
				y_batch = y[batch_idx]

				∇x, ∇b = batch_gradients(A_batch, y_batch, x, b)
				x -= α*∇x
				b -= α*∇b
			end

			# Check MSE after each epoch
			y_pred = predict(A, x, b)
			mse = mean_squared_error(y, y_pred)
			push!(mse_history, mse)
		end

		return x, b, mse_history
	end

	function sgd_batch1(A, y; α = 0.001, n_epochs=50)
		x = zeros(size(A, 2))
		b = 0.0
		N = size(A, 1)
		mse_history = Float64[]

		for epoch in 1:n_epochs
			indices = randperm(N) 			# Shuffle data
			for i in indices
				A_i = A[i, :]  				# shape (1,2) Note A[1, : ] is a col vec!
				if A_i isa AbstractVector
					A_i = A_i'
				end
				y_i = y[i] 					# shape (1, )
				∇x, ∇b =  batch_gradients(A_i, y_i, x, b)
				x -= α * ∇x                # element-wise
				b -= α * ∇b
			end

			# MSE after epoch
			y_pred = predict(A, x, b)
			mse = mean_squared_error(y, y_pred)
			push!(mse_history, mse)
		end

		return x, b, mse_history
	end
end

# ╔═╡ 1fb81a1d-7cbb-4ac2-b64b-69b5659e0527
begin
	xfull, bfull, mse_history_full = full_batch_gd(A, y)
	xmini, bmini, mse_history_mini = mini_batch_sgd(A, y)
	x_sdg, b_sdg, mse_history_sdg = sgd_batch1(A, y)
end

# ╔═╡ 9f4cb262-9c50-42d7-873a-9782ff99ade8
xfull, xmini, x_sdg

# ╔═╡ 371ae476-99ea-44c3-b115-eadd77557e65
let
	global fig1 = Figure(size=(600, 500)); nothing
end

# ╔═╡ 3e2a96ea-8988-4703-9726-c48d7a3c80e7
let 
	empty!(fig1)

	ax2 = CairoMakie.Axis(fig1[1,1])
	scatterlines!(ax2, mse_history_full, label = L"$$ GD (full-batch)")
	scatterlines!(ax2, mse_history_mini, label = L"$$ Mini-Batch GD")
	scatterlines!(ax2, mse_history_sdg, label = L"$$ SGD (batch=1)")
	axislegend(ax2, labelsize=14)
	
	ax1 = CairoMakie.Axis(fig1[2,1], yscale=log10,
						  yminorticksvisible = true, yminorgridvisible = true, 
						 yminorticks = IntervalsBetween(10))
	scatterlines!(ax1, mse_history_full, label = L"$$ GD (full-batch)")
	scatterlines!(ax1, mse_history_mini, label = L"$$ Mini-Batch GD")
	scatterlines!(ax1, mse_history_sdg, label = L"$$ SGD (batch=1)")
	axislegend(ax1, labelsize=14)

	fig1
end

# ╔═╡ bce65c1b-1786-453e-a28f-bf902f7f29d6
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient.pdf", fig1)
	input = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient.pdf"
	output = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient.pdf"
	run(`pdfcrop $input $output`)
end
  ╠═╡ =#

# ╔═╡ f2c59982-1027-49b9-ae16-5479d2bf129d
begin
	@time xmini_1, bmini_1, mse_history_mini_1 = mini_batch_sgd(A, y; batch_size = 1)
	@time xmini_10, bmini_10, mse_history_mini_10 = mini_batch_sgd(A, y; batch_size = 10)
	@time xmini_50, bmini_50, mse_history_mini_50 = mini_batch_sgd(A, y; batch_size = 50)
	@time xmini_200, bmini_200, mse_history_mini_200 = mini_batch_sgd(A, y; batch_size = 200)
	@time xmini_2000, bmini_2000, mse_history_mini_2000 = mini_batch_sgd(A, y; batch_size = 2000)
end

# ╔═╡ 9defc413-49b8-41bb-864a-4c351734fc35
let
	global fig2 = Figure(size=(600, 300))
	ax1 = CairoMakie.Axis(fig2[1,1])
	scatterlines!(ax1, mse_history_mini_1, label = L"$$ Batch size = 1")
	scatterlines!(ax1, mse_history_mini_10, label = L"$$ Batch size = 10")
	scatterlines!(ax1, mse_history_mini_50, label = L"$$ Batch size = 50")
	scatterlines!(ax1, mse_history_mini_200, label = L"$$ Batch size = 200")
	scatterlines!(ax1, mse_history_mini_2000, label = L"$$ Batch size = 2000")
	axislegend(ax1, labelsize=14)

	fig2
end

# ╔═╡ 5d91442f-81bd-41b3-b3a6-a9a1cff02a2a
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient2.pdf", fig2)
	input = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient2.pdf"
	output = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient2.pdf"
	run(`pdfcrop $input $output`)
end
  ╠═╡ =#

# ╔═╡ 77bd4fc3-2079-4c1b-a3e3-34e03e78d818
md"""
## Example: SGD on a Nonlinear Regression Problem
"""

# ╔═╡ 22b98f00-d849-48e4-bdfd-88fdbba4f930
let
	@variables x
	
	latexstring("f(x) = " * latexify(-10x^2 + 2x^4))
end

# ╔═╡ 3b402d95-cfc8-440b-ac5c-3855aa1495dd
begin
	function true_poly(x)
		return (-10.0x.^2 + 2x.^4)
	end

	num_train = 50 				# Small training set makes the problem challenging
	num_test = 50

	x_train = collect(LinRange(-2, 2, num_train)) 	# data in range [-2, 2]
	y_train_true = true_poly(x_train)

	#noise_std = 2.0 			# Add relatively large noise
	noise_train = noise_std * randn(num_train)
	y_train = y_train_true .+ noise_train

	Random.seed!(2023) 					# Test data (different seed)
	x_test = collect(LinRange(-2, 2, num_test))
	y_test_true = true_poly(x_test)
	noise_test = noise_std .+ randn(num_test)
	y_test = y_test_true .+ noise_test

end

# ╔═╡ 4ccc09af-2455-48f7-9d59-02be3a291820
begin
	function create_poly_features(x::AbstractVector; degree::Int=5)
		# return matrix with columns [x^0, x^1, ..., x^degree]
		return hcat([x.^i for i in 0:degree]...)  # we need ... !!!!!
	end
	
	X_train = create_poly_features(x_train)
	X_test = create_poly_features(x_test)	

	function batch_gradient(X, y, w)
		# Compute gradient of MSE w.r.t. w: -2/N *X'(y - Xw)
		N = size(X,1)
		residuals = y .- X * w
		
		∇w = -(2.0 / N) * X' * residuals
	
		return ∇w
	end
end

# ╔═╡ dad3194a-f716-4143-bae2-a705ab988381
begin
	function mini_batch_sgd_only_X(A, y; batch_size = 1, α = 0.01, n_epochs=1000)
		x = zeros(size(A, 2)) 
		N = size(A, 1)
		mse_history = Float64[]

		for epoch in 1:n_epochs
			indices = randperm(N)			# Shuffle data
			for start in 1:batch_size:N
				last = min(start + batch_size - 1, N)  # avoid exceed array size
				batch_idx = indices[start:last]
				A_batch = A[batch_idx, :]
				if size(A_batch, 2) == 1
					A_batch = A_batch'
				end
				y_batch = y[batch_idx]
				∇w = batch_gradient(A_batch, y_batch, x)
				x -= α*∇w
			end

			# Check MSE after each epoch
			y_pred = predict(A, x, zeros(size(A, 1)))
			mse = mean_squared_error(y, y_pred)
			push!(mse_history, mse)
		end

		return x, mse_history
	end

	function learning_schedule(t, t0, t1)
		return t0 /(t + t1)
	end
	
	function stochastic_gradient(A, y; n_epochs = 50)
		N = size(A, 1)
		mse_history = Float64[]


		#Random.seed!(42)				# to make same noise
		x = vec(randn(size(A, 2), 1))	# randomly initialized model parameter
		for epoch in 1:n_epochs
			indices = randperm(N) 		# Shuffle data
			for iteration in indices
				A_i = A[iteration, :]
				if A_i isa AbstractVector 
					A_i = A_i' 				# Change column to row vector
				end
				y_i = y[iteration] 		
				∇x = (2) * A_i' * (A_i * x .- y_i)     	# note N = 1
				# n_epochs is a tuning parameter
				η = learning_schedule(n_epochs + iteration, 1, n_epochs)
				x -= η * ∇x
			end
			# Check MSE after each epoch
			y_pred = predict(A, x, zeros(size(A, 1)))
			mse = mean_squared_error(y, y_pred)
			push!(mse_history, mse)
		end
			
		return x, mse_history		
	end
end

# ╔═╡ 105a146a-1435-4169-8630-4d678c6794be
function adam_gradient(A, y; n_epochs = 50, α = 0.01, β1 = 0.9, β2 = 0.999, ϵ = 1e-8)
    N = size(A, 1)
    d = size(A, 2)
    mse_history = Float64[]

    x = vec(randn(d, 1))         # Initial model parameters 
    m = zeros(d)                 # First moment estimate 
    v = zeros(d)                 # Second moment estimate 

    t = 0                        # Time step 

    for epoch in 1:n_epochs
        indices = randperm(N)    # Shuffle data 
        for iteration in indices
            t += 1
            A_i = A[iteration, :]
            if A_i isa AbstractVector
                A_i = A_i'
            end
            y_i = y[iteration]

            ∇x = 2 * A_i' * (A_i * x .- y_i)  # Gradient of squared error 

            # Update biased first and second moment estimates
            m = β1 * m + (1 - β1) * ∇x
            v = β2 * v + (1 - β2) * (∇x .^ 2)

            # Bias correction
            m̂ = m ./ (1 - β1^t)
            v̂ = v ./ (1 - β2^t)

            # Parameter update
            x -= α * m̂ ./ (sqrt.(v̂) .+ ϵ)
        end

        # Evaluate MSE after each epoch
        y_pred = predict(A, x, zeros(N))
        mse = mean_squared_error(y, y_pred)
        push!(mse_history, mse)
    end

    return x, mse_history
end

# ╔═╡ 272ade81-b8f1-428d-a744-2ccd8c756260
begin
	x_mini_nlp, mse_history_nlp = mini_batch_sgd_only_X(X_train, y_train, 
											batch_size = 50, α = 0.008)
	x_sto, mse_history_sto = stochastic_gradient(X_train, y_train;
												 n_epochs=1000)
	x_adam, mse_history_adam = adam_gradient(X_train, y_train; 
											 n_epochs = 1000, α = 0.01)
	
	y_mini_nlp_test = X_test * x_mini_nlp; 
	y_mini_nlp_train = X_train * x_mini_nlp; 
	println(x_mini_nlp)
	
	y_sto_test = X_test * x_sto;
	y_sto_train = X_train * x_sto;
	println(x_sto)

	y_adam_test = X_test * x_adam
	y_adam_train = X_train * x_adam
	println(x_adam)
end

# ╔═╡ c7cbf19c-efdf-49b4-86c9-cd58615322df
let
	global fig3 = Figure(size=(600,600))

	# train data
	ax2 = CairoMakie.Axis(fig3[1,1], xlabel=L"$x$", ylabel=L"$y$")
	scatter!(ax2, y_train, strokewidth=1, strokecolor=:black, 
			 label=L"$$ Train (noise)")
	lines!(ax2, y_train_true, linestyle=:dash, color=:black, linewidth=3,
		  label=L"$$ True function")
	lines!(ax2, y_mini_nlp_train, linewidth = 3, label=L"$$ Full batch GD")
	lines!(ax2, y_sto_train, linewidth = 3,  label=L"$$ SGD")
	lines!(ax2, y_adam_train, linewidth = 3, linestyle=:dash, label=L"$$ ADAM", color=:red)
	
	axislegend(ax2, labelsize=14)
	
	# test data
	ax1 = CairoMakie.Axis(fig3[2,1], xlabel=L"$x$", ylabel=L"$y$")
	scatter!(ax1, y_test, strokewidth=1, strokecolor=:black, 
			 label=L"$$ Test (noise)")
	lines!(ax1, y_test_true, linestyle=:dash, color=:black, linewidth=3,
		  label=L"$$ True function")
	lines!(ax1, y_mini_nlp_test, linewidth = 3, label=L"$$ Full batch GD")
	lines!(ax1, y_sto_test, linewidth = 3,  label=L"$$ SGD")
	lines!(ax1, y_adam_test, linewidth = 3, linestyle=:dash, 
		   label=L"$$ ADAM", color=:red)
	
	axislegend(ax1, labelsize=14)
	fig3
end

# ╔═╡ a7c60771-57e3-4b33-8181-544129f86f00
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient3.pdf", fig3)
	input = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient3.pdf"
	output = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient3.pdf"
	run(`pdfcrop $input $output`)
end
  ╠═╡ =#

# ╔═╡ fa8fce40-2427-43c1-87e8-955adb1c713d
let
	global fig4 = Figure(size=(600, 300))
	
	ax1 = CairoMakie.Axis(fig4[1,1], yscale=log10,
					yminorticksvisible = true, yminorgridvisible = true, 
					yminorticks = IntervalsBetween(10), 
					xlabel=L"$$ Iteration", ylabel=L"$$ MSE")
	#limits!(ax1, -20, 1000, -10, 50)
	lines!(ax1, mse_history_adam, label = L"$$ ADAM", color=:red)
	lines!(ax1, mse_history_nlp, label = L"$$ Batch size = 50")
	lines!(ax1, mse_history_sto, label = L"$$ SGD")
	

	axislegend(ax1, labelsize=14)

	fig4
end

# ╔═╡ c51cf85c-eaa4-4837-951b-fd4363a86db6
# ╠═╡ disabled = true
#=╠═╡
begin
	save("/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient4.pdf", fig4)
	input = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient4.pdf"
	output = "/mnt/e/OneDrive/Public/workKMUTT/INC Selection Optimization/Lecture2022/images/batch_gradient4.pdf"
	run(`pdfcrop $input $output`)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═f0f4f39c-8689-11f0-169d-2fe53b79ca12
# ╠═def70522-1384-4325-8f91-9dfb09681b0f
# ╟─ab0fcc0c-dd93-48e2-a844-28acb853564e
# ╠═73ffb1e3-e328-4925-838d-504e78748cee
# ╠═7460106e-9309-451b-8672-a01a0a9655f1
# ╠═fd7872af-3085-49ba-bf7f-97fc3761075f
# ╠═1fb81a1d-7cbb-4ac2-b64b-69b5659e0527
# ╠═9f4cb262-9c50-42d7-873a-9782ff99ade8
# ╠═371ae476-99ea-44c3-b115-eadd77557e65
# ╠═3e2a96ea-8988-4703-9726-c48d7a3c80e7
# ╠═bce65c1b-1786-453e-a28f-bf902f7f29d6
# ╠═f2c59982-1027-49b9-ae16-5479d2bf129d
# ╠═9defc413-49b8-41bb-864a-4c351734fc35
# ╠═5d91442f-81bd-41b3-b3a6-a9a1cff02a2a
# ╟─77bd4fc3-2079-4c1b-a3e3-34e03e78d818
# ╟─22b98f00-d849-48e4-bdfd-88fdbba4f930
# ╠═3b402d95-cfc8-440b-ac5c-3855aa1495dd
# ╠═4ccc09af-2455-48f7-9d59-02be3a291820
# ╠═dad3194a-f716-4143-bae2-a705ab988381
# ╠═105a146a-1435-4169-8630-4d678c6794be
# ╠═272ade81-b8f1-428d-a744-2ccd8c756260
# ╠═c7cbf19c-efdf-49b4-86c9-cd58615322df
# ╠═a7c60771-57e3-4b33-8181-544129f86f00
# ╠═fa8fce40-2427-43c1-87e8-955adb1c713d
# ╠═c51cf85c-eaa4-4837-951b-fd4363a86db6
