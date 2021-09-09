### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 80b5ea60-ca67-11eb-278f-7d8e19a1b996
begin
	using LinearAlgebra, DifferentialEquations, Plots
	md"""
	## Astrocamp Module 2.2 Coding Assignment
	### Problem 1
	"""
end

# ╔═╡ eaedda48-5520-4592-9018-3186be1f8455
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

# ╔═╡ 2bbf72d3-c0a3-439d-9fb3-954b3bed8cea
function pseudoPotential(X, μ)
	# FINDING AUXILIARY DISTANCES
	r⃗13 = X[1:3] - [-μ, 0, 0]
	r13 = norm(r⃗13)
	r⃗23 = X[1:3] - [1-μ, 0, 0]
	r23 = norm(r⃗23)
	
	U = (1-μ)/r13 + μ/r23 + norm(X[1:3])^2/2
end

# ╔═╡ 2dafe87f-8bdb-4be7-8aa3-dcb3e02a37b2
begin
	# CONSTANTS
	μ = 0.0121505856   # Nondimm System μ
	X = [1.1745 0 0 0 -0.111503288786392 0]
	J₀ = 2*pseudoPotential(X, μ) - norm(X[4:6])^2
	tspan = (0, 3.393216168933766)
	
	# ODE PROBLEM
	prob = ODEProblem(CR3BP!, X, tspan, μ)
	sol = solve(prob)#, reltol=1e-8, abstol=1e-8)
	md"Initial Jacobi Constant: $J₀"
end

# ╔═╡ 3cd63cae-79fd-4b89-9a0f-805b974f30ef
begin
	gr()
	Lₓ = [0.836915 1.15568 -1.00506] # QUICK AND DIRTY L-POINT X POS
	plot([1-μ], [0], lw=0, marker=:circle, color=:black, label="Moon")
	plot!(sol, vars=(1,2), lw=2, label="Spacecraft", aspect_ratio=:equal, xlim=(0.95, 1.25))
	plot!([Lₓ[2]], [0], marker=:circle, color=:red, label="L₂", title="Plot of Orbit")
end

# ╔═╡ 942c041d-52a2-424b-9d46-5b98c46b960e
begin
	J = Any[]
	for i = 1:length(sol.u)
		push!(J, 2*pseudoPotential(sol.u[i], μ) - norm(sol.u[i][4:6])^2)
	end
	plot(sol.t, J.-J₀, lw=0.25, label="", marker=:dot, 
		xlabel="Time", ylabel="Error", title="Jacobi Constant Integration Error")
end

# ╔═╡ Cell order:
# ╠═80b5ea60-ca67-11eb-278f-7d8e19a1b996
# ╠═eaedda48-5520-4592-9018-3186be1f8455
# ╠═2bbf72d3-c0a3-439d-9fb3-954b3bed8cea
# ╠═2dafe87f-8bdb-4be7-8aa3-dcb3e02a37b2
# ╠═3cd63cae-79fd-4b89-9a0f-805b974f30ef
# ╠═942c041d-52a2-424b-9d46-5b98c46b960e
