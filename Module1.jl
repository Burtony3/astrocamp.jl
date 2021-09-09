### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 95d9150a-31d6-4d7f-b287-016bba8ba14a
begin
	using LinearAlgebra, DifferentialEquations, Plots
	md"## Astrocamp Module 1 Coding Assignment"
end

# ╔═╡ 43da0930-c717-11eb-26e2-b1ab3489dde7
md"### Problem 1"

# ╔═╡ a22c79fa-5307-4912-9d82-ddba8cfa1476
function twobody!(du, u, p, t)
	# REDEFINING EQUATION CAUSES DIFFEQ.jl to SLOW CONSIDERABLY
	μ, = p
	r⃗ = u[1:3]
	v⃗ = u[4:6]
	du[1:3] = v⃗
	du[4:6] = -(μ/norm(r⃗)^3)*r⃗
end

# ╔═╡ 2b65bd0e-7a86-4283-aea6-7e3ce9eb2732
begin	
	# INPUTS
	R⃗ = [-7327.031 -813.869 0.]
	V⃗ = [1.137 -10.237 0.]
	tspan = (0., 40. * 86400)
	μₑ = 398600
	
	# SETTING UP ODE PROBLEM
	params1 = (μₑ)
	prob1 = ODEProblem(twobody!, [R⃗ V⃗], tspan, params1)
	
end

# ╔═╡ a293b549-96ac-4690-b823-5a9f3a86a9c8
begin
	sol1 = solve(prob1,reltol=1e-8,abstol=1e-8)
	sol1[end]	
end

# ╔═╡ 96f8dd53-4fe7-4ebc-8251-8325a6b01d08
begin
	gr()
	plot(sol1,vars=(1,2), lw=2, label="Spacecraft")
	θMoon2bd = LinRange(0, 2*π, 100)
	xMoon2bd = 385000*cos.(θMoon2bd) # cos. is required to broadcast across range
	yMoon2bd = 385000*sin.(θMoon2bd)
	plot!(xMoon2bd, yMoon2bd, aspect_ratio=:equal, label="Moon", 
		xlim=385000*[-1.01 1.01])
end

# ╔═╡ 79a32517-dc62-41ae-9276-4053794bf5fd
md"### Problem 2"

# ╔═╡ e47ea334-493d-4fbe-bb80-9fe8e1b89d96
function threebody_intertial!(du, u, p, t)
	# DECODING INPUTS
	r⃗ = u[1:3]
	v⃗ = u[4:6]
	r⃗1 = u[7:9]
	v⃗1 = u[10:12]
	r⃗2 = u[13:15]
	v⃗2 = u[16:18]
	
	# CALCULATING IMPORTANT VALUES
	μ1, μ2 = p
	r⃗13 = r⃗ - r⃗1
	r⃗23 = r⃗ - r⃗2
	θ1 = atan(r⃗1[2], r⃗1[1])
	a1 = norm(v⃗1)^2/norm(r⃗1)
	θ2 = atan(r⃗2[2], r⃗2[1])
	a2 = norm(v⃗2)^2/norm(r⃗2)
	
	# EXPORTING
	du[1:3]   = v⃗
	du[4:6]   = -(μ1/norm(r⃗13)^3)*r⃗13 - (μ2/norm(r⃗23)^3)*r⃗23
	du[7:9]   = v⃗1
	du[10:12] = a1*[-cos(θ1) -sin(θ1) 0]
	du[13:15] = v⃗2
	du[16:18] = a2*[-cos(θ2) -sin(θ2) 0]
end

# ╔═╡ bde40216-aa82-4365-a65f-1c9d9c34149b
begin
	# USING SAME AS ABOVE INPUTS
	
	# SETTING UP ODE
	Rₑ = 4677.975
	R⃗ₑ = Rₑ*[-1 0 0]
	V⃗ₑ = [0 -0.013 0]
	μₘ = 4.903e3
	ωₑ = sqrt((μₑ+μₘ)/Rₑ^3)
	Rₘ = 380322.025
	R⃗ₘ = Rₘ*[1 0 0]
	V⃗ₘ = [0 1.012 0]
	ωₘ = sqrt((μₑ+μₘ)/Rₘ^3)
	params2 = (μₑ, μₘ)
	IC =[R⃗+R⃗ₑ V⃗+V⃗ₑ R⃗ₑ V⃗ₑ R⃗ₘ V⃗ₘ]
	# IC = [R⃗ V⃗]
	prob2 = ODEProblem(threebody_intertial!, IC, tspan, params2)
end

# ╔═╡ 99f2a869-5959-436e-9ed2-db6b091b61f9
begin
	sol2 = solve(prob2,reltol=1e-8,abstol=1e-8)
	sol2[end][1:6]
end

# ╔═╡ 9621b40d-51d0-489c-b0d2-6ceea25dd2fd
begin
	plot(sol2, vars=(1, 2), lw=2, label="Spacecraft")
	plot!(sol2,vars=(7, 8), label="Earth")
	plot!(sol2,vars=(13, 14), label="Moon", aspect_ratio=:equal)
end

# ╔═╡ a3dc528e-36b2-40f3-86c6-afb074440291
md"### Problem 3"

# ╔═╡ 51d5a30f-2cc9-426d-b9e1-d1e335815549
begin
	# GETTING NON-DIMENSIONALIZED SCALARS & TRANSFORMATION MATRICIES
	ls = Rₑ + Rₘ
	ts = sqrt(ls^3 / (μₑ + μₘ))
	vs = ls/ts
	
	# CALCULATING NEW VALUES
	r⃗ = (R⃗+R⃗ₑ)/ls
	v⃗ = (V⃗)/vs
	[r⃗ v⃗]
end

# ╔═╡ 7d48dcf5-7c99-4303-beb5-75a31b2c8d31
md"### Problem 4"

# ╔═╡ 4aecfe26-fcce-428a-8470-a3b368d029d5
function threebody_rotational!(du, u, p, t)
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
	du[4] = 2*v⃗[2] + r⃗[1] - ((1-μ)/r13^3)*(r⃗[1] + μ) - (μ/r23^3)*(r⃗[1] + μ - 1)
	du[5] = -2*v⃗[1] + r⃗[2] - ((1-μ)/r13^3)*(r⃗[2]) - (μ/r23^3)*(r⃗[2])
	du[6] = ((1-μ)/r13^3)*(r⃗[3]) - (μ/r23^3)*(r⃗[3])
	# du[4:6] = -((1-μ)/r13^3)*r⃗13 - (μ/r23^3)*r⃗23
end

# ╔═╡ 1d4dfca9-19d8-4854-bd51-7f920acd8246
begin	
	# INPUTS
	tspan4 = (0, 40*86400/ts)
	
	# SETTING UP ODE PROBLEM
	params4 = (μₘ/(μₑ + μₘ))
	prob4 = ODEProblem(threebody_rotational!, [-0.03118 -0.002114 0 1.108 -9.974 0], tspan4, params4)
	
end

# ╔═╡ 0351d545-a00b-4df3-bc02-3cba95d1b678
begin
	sol4 = solve(prob4,reltol=1e-8,abstol=1e-8)
	sol4[end]
end

# ╔═╡ 3e42b9e1-392e-4d94-94a8-83964a9a7dba
begin
	plot(sol4, vars=(1, 2), lw=2, label="Spacecraft")
	plot!([1-params4[1]], [0], lw=0, marker=:hexagon, label="Moon")
	plot!([-params4[1]], [0], lw=0, marker=:hexagon, label="Earth", 
		aspect_ratio=:equal, xlim=[-1.5 1.5])
end

# ╔═╡ 785dc94a-0a4b-4773-80e9-2037854e019c
md"### Problem 5"

# ╔═╡ 14ebcdbe-dfcb-4e0a-a96c-27b079964333
begin
	# CREATING ANONYMOUS FUNCTIONS
	A(τ) = [cos(τ) -sin(τ) 0;
			sin(τ) cos(τ) 0;
			0 0 1]
	Ȧ(τ) = [-sin(τ) -cos(τ) 0;
			cos(τ) -sin(τ) 0;
			0 0 0]
	r2R(r⃗, τ) = ls*A(τ)*r⃗
	v2V(v⃗, r⃗, τ) = vs*Ȧ(τ)*r⃗ + vs*A(τ)*v⃗
	
	# PRE-ALLOCATION
	R⃗new = zeros(3, length(sol4))
	# V⃗new = R⃗new # This will NOT work, as this is a referential relation any changes 
	# 				to R⃗new will affect V⃗new in later lines
	V⃗new = zeros(3, length(sol4))
	
	# CONVERTING
	for i = 1:length(sol4)
		R⃗new[:, i] = r2R(sol4.u[i][1:3], sol4.t[i])
		V⃗new[:, i] = v2V(sol4.u[i][4:6], sol4.u[i][1:3], sol4.t[i])
	end
	[R⃗new; V⃗new]
	# r2R(sol4.u[end][1:3], sol4.t[end])
	# v2V(sol4.u[end][4:6], sol4.u[end][1:3], sol4.t[end])
end

# ╔═╡ c8886267-beff-4c66-ab66-eb6f532b13af
begin
	plot(R⃗new[1, :], R⃗new[2, :], lw=2, label="Spacecraft", aspect_ratio=:equal)
end

# ╔═╡ Cell order:
# ╠═95d9150a-31d6-4d7f-b287-016bba8ba14a
# ╟─43da0930-c717-11eb-26e2-b1ab3489dde7
# ╠═a22c79fa-5307-4912-9d82-ddba8cfa1476
# ╠═2b65bd0e-7a86-4283-aea6-7e3ce9eb2732
# ╠═a293b549-96ac-4690-b823-5a9f3a86a9c8
# ╠═96f8dd53-4fe7-4ebc-8251-8325a6b01d08
# ╟─79a32517-dc62-41ae-9276-4053794bf5fd
# ╠═e47ea334-493d-4fbe-bb80-9fe8e1b89d96
# ╠═bde40216-aa82-4365-a65f-1c9d9c34149b
# ╠═99f2a869-5959-436e-9ed2-db6b091b61f9
# ╠═9621b40d-51d0-489c-b0d2-6ceea25dd2fd
# ╟─a3dc528e-36b2-40f3-86c6-afb074440291
# ╠═51d5a30f-2cc9-426d-b9e1-d1e335815549
# ╟─7d48dcf5-7c99-4303-beb5-75a31b2c8d31
# ╠═4aecfe26-fcce-428a-8470-a3b368d029d5
# ╠═1d4dfca9-19d8-4854-bd51-7f920acd8246
# ╠═0351d545-a00b-4df3-bc02-3cba95d1b678
# ╠═3e42b9e1-392e-4d94-94a8-83964a9a7dba
# ╟─785dc94a-0a4b-4773-80e9-2037854e019c
# ╠═14ebcdbe-dfcb-4e0a-a96c-27b079964333
# ╠═c8886267-beff-4c66-ab66-eb6f532b13af
