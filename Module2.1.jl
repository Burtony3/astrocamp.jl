### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ c8900b10-ca5b-11eb-06e8-ad64264938f4
begin
	using LinearAlgebra, DifferentialEquations, Plots
	md"""
	## Astrocamp Module 2.1 Coding Assignment
	### Problem 1
	"""
end

# ╔═╡ 703f5fbe-ed47-48ef-9d9c-6104b0bcbb77
begin
	using Random
	# CONSTANTS
	μ = 0.0121505856   # Nondimm System μ
	r⃗₀ = [0.8369 0 0]  # Initial Position (L₁)
	v⃗₀ = [-0.025 0.05 0]   # Initial Velocity
	σ = 1e-2           # Velocity Perturbation Stnd. Dev.
	iters = 100
	
	# PROBLEM SETUP
	genIC() = [r⃗₀ v⃗₀+σ*[rand([-1.0 1.0])*randn() rand([-1.0 1.0])*randn() 0]]
	tspan = (0, 25.0)
end

# ╔═╡ da6c88c9-0fb4-4669-8753-1ed0fcc3adcf
function CR3BP!(du, u, p, t)
	# DECODING INPUTS
	r⃗ = u[1:3]
	v⃗ = u[4:6]
	μ = p
	
	# CALCULATING MIDDLE VALUES
	r⃗13 = r⃗ - [-μ, 0, 0]
	r13 = norm(r⃗13)
	r⃗23 = r⃗ - [1-μ, 0, 0]
	r23 = norm(r⃗23)
	
	# EXPORTING
	du[1:3] = v⃗
	du[4:6] = -((1-μ)/r13^3)*r⃗13 - (μ/r23^3)*r⃗23 + [2*v⃗[2]+r⃗[1], -2*v⃗[1]+r⃗[2], 0]
end

# ╔═╡ 1ac9ecf1-5f5d-48a8-b7fb-f77e100f2a02
begin
	# GENERATING STORAGE
	out = Any[]          # State Arrays
	ICs = Any[]		 	 # Initial Conditions
	
	# RUNNING LOOP
	for i = 1:iters
		# INTEGRATING
		IC = genIC()
		prob = ODEProblem(CR3BP!, IC, tspan, μ)
		sol = solve(prob)
		
		# STORING RESULTS
		push!(out, vcat(sol.u...))    # "..." turns each indici into an input
		push!(ICs, IC)
	end
end

# ╔═╡ 60f16607-12e7-408d-bb3b-3372f561b946
begin
	gr()
	plot([r⃗₀[1]], [r⃗₀[2]], marker=:circle, aspect_ratio=:equal, label="L₁")
	plot!(out[1][:, 1], out[1][:, 2], lw=0.5, alpha=0.1, color=:grey, label="Trajectories")
	for i = 2:iters
		plot!(out[i][:, 1], out[i][:, 2], lw=0.5, alpha=0.1, color=:grey, label="")
	end
	plot!([1-μ], [0], lw=0, marker=:circle, color=:black, label="Bodies")
	plot!([-μ], [0], lw=0, marker=:circle, color=:black, label="")
end

# ╔═╡ Cell order:
# ╠═c8900b10-ca5b-11eb-06e8-ad64264938f4
# ╠═da6c88c9-0fb4-4669-8753-1ed0fcc3adcf
# ╠═703f5fbe-ed47-48ef-9d9c-6104b0bcbb77
# ╠═1ac9ecf1-5f5d-48a8-b7fb-f77e100f2a02
# ╠═60f16607-12e7-408d-bb3b-3372f561b946
