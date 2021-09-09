### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 233ec700-ccd1-11eb-2983-914a848c6aab
begin
	using LinearAlgebra, DifferentialEquations, Plots
	md"""
	## Astrocamp Module 2 Coding Assignment
	### Problem 1
	"""
end

# ╔═╡ 366664b7-6b31-4ae7-b479-58a53a717ffb
md"""
#### Deriving Jacobian Matrix

Function for Psuedo-Potential Energy:

$U = \dfrac{1-\mu}{\sqrt{(x+\mu)^2 + y^2 + z^2}} + \dfrac{\mu}{\sqrt{(x+\mu-1)^2 + y^2 + z^2}} + \dfrac{1}{2}(x^2 + y^2 + z^2)$

Derivatives w.r.t. each independent variable

$\dfrac{\partial U}{\partial x} = U_x = -\dfrac{\left(1-{\mu}\right)\left(x+{\mu}\right)}{\left(\left(x+{\mu}\right)^2+z^2+y^2\right)^\frac{3}{2}}-\dfrac{{\mu}\left(x+{\mu}-1\right)}{\left(\left(x+{\mu}-1\right)^2+z^2+y^2\right)^\frac{3}{2}}+x$

$\dfrac{\partial U}{\partial y} = U_y = -\dfrac{\left(1-{\mu}\right)y}{\left(\left(x+{\mu}\right)^2 + y^2+z^2\right)^\frac{3}{2}}-\dfrac{{\mu}y}{\left(\left(x+{\mu}-1\right)^2 + y^2+z^2\right)^\frac{3}{2}}+y$

$\dfrac{\partial U}{\partial z} = U_z = -\dfrac{\left(1-{\mu}\right)z}{\left(\left(x+{\mu}\right)^2 + y^2+z^2\right)^\frac{3}{2}}-\dfrac{{\mu}z}{\left(\left(x+{\mu}-1\right)^2 + y^2+z^2\right)^\frac{3}{2}}+z$
"""

# ╔═╡ cacfe7ff-5939-4be7-9a28-08829d9e7e11
function jacobian(u, μ)
	# UNWRAPPING INPUTS
	x, y, z = u[1:3]
	r⃗13 = [x+μ,   y, z]; r13 = norm(r⃗13)
	r⃗23 = [x+μ-1, y, z]; r23 = norm(r⃗23)
	#=
	# Derivatives of Ux
	Uxx = -(1-μ)/((x+μ)^2+z^2+y^2)^(3/2)+(3*(1-μ)*(x+μ)^2)/((x+μ)^2+z^2+y^2)^(5/2)-μ/((x+μ-1)^2+z^2+y^2)^(3/2)+(3*μ*(x+μ-1)^2)/((x+μ-1)^2+z^2+y^2)^(5/2)+1
	Uxy = (3*(1-μ)*(x+μ)*y)/(y^2+z^2+(x+μ)^2)^(5/2)+(3*μ*(x+μ-1)*y)/(y^2+z^2+(x+μ-1)^2)^(5/2)
	Uxz = (3*(1-μ)*(x+μ)*z)/(z^2+y^2+(x+μ)^2)^(5/2)+(3*μ*(x+μ-1)*z)/(z^2+y^2+(x+μ-1)^2)^(5/2)
	
	# Derivatives of Uy
	Uyx = (3*(1-μ)*y*(x+μ))/((x+μ)^2+z^2+y^2)^(5/2)+(3*μ*y*(x+μ-1))/((x+μ-1)^2+z^2+y^2)^(5/2)
	Uyy = -(1-μ)/(y^2+z^2+(x+μ)^2)^(3/2)+(3*(1-μ)*y^2)/(y^2+z^2+(x+μ)^2)^(5/2)-μ/(y^2+z^2+(x+μ-1)^2)^(3/2)+(3*μ*y^2)/(y^2+z^2+(x+μ-1)^2)^(5/2)+1
	Uyz = (3*(1-μ)*y*z)/(z^2+y^2+(x+μ)^2)^(5/2)+(3*μ*y*z)/(z^2+y^2+(x+μ-1)^2)^(5/2)
	
	# Derivatives of Uz
	Uzx = (3*(1-μ)*z*(x+μ))/((x+μ)^2+z^2+y^2)^(5/2)+(3*μ*z*(x+μ-1))/((x+μ-1)^2+z^2+y^2)^(5/2)
	Uzy = (3*(1-μ)*z*y)/(y^2+z^2+(x+μ)^2)^(5/2)+(3*μ*z*y)/(y^2+z^2+(x+μ-1)^2)^(5/2)
	Uzz = -(1-μ)/(z^2+y^2+(x+μ)^2)^(3/2)+(3*(1-μ)*z^2)/(z^2+y^2+(x+μ)^2)^(5/2)-μ/(z^2+y^2+(x+μ-1)^2)^(3/2)+(3*μ*z^2)/(z^2+y^2+(x+μ-1)^2)^(5/2)+1
	=#

	# Derivatives of Ux
	Uxx = -(1-μ)/r13^3+(3*(1-μ)*(x+μ)^2)/r13^5 - μ/r23^3+(3*μ*(x+μ-1)^2)/r23^5 + 1
	Uxy = (3*(1-μ)*(x+μ)*y)/r13^5 + (3*μ*(x+μ-1)*y)/r23^5
	Uxz = (3*(1-μ)*(x+μ)*z)/r13^5 + (3*μ*(x+μ-1)*z)/r23^5
	
	# Derivatives of Uy
	Uyx = (3*(1-μ)*y*(x+μ))/r13^5 + (3*μ*y*(x+μ-1))/r23^5
	Uyy = -(1-μ)/r13^3 + (3*(1-μ)*y^2)/r13^5 - μ/r23^3 + (3*μ*y^2)/r23^5 + 1
	Uyz = (3*(1-μ)*y*z)/r13^5 + (3*μ*y*z)/r23^5
	
	# Derivatives of Uz
	Uzx = (3*(1-μ)*z*(x+μ))/r13^5 + (3*μ*z*(x+μ-1))/r23^5
	Uzy = (3*(1-μ)*z*y)/r13^5 + (3*μ*z*y)/r23^5
	Uzz = -(1-μ)/r13^3 + (3*(1-μ)*z^2)/r13^5 - μ/r23^3 + (3*μ*z^2)/r23^5 + 1

	#CONSTRUCTING JACOBIAN
	A = [ 0   0   0   1   0  0;
		  0   0   0   0   1  0;
		  0   0   0   0   0  1;
		 Uxx Uxy Uxz  0   2  0;
		 Uyx Uyy Uyz  -2  0  0;
		 Uzx Uzy Uzz  0   0  0 ];
	
end

# ╔═╡ 6b72128e-e2bb-47a5-8a80-315d70348923
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
	
	# EXPORTING TRAJECTORY EOM
	du[1:3] = v⃗
	du[4:6] = -((1-μ)/r13^3)*r⃗13 - (μ/r23^3)*r⃗23 + [2*v⃗[2]+r⃗[1], -2*v⃗[1]+r⃗[2], 0]
	
	# CALCULATING STM
	A = jacobian(u, μ)
	ϕ = reshape(u[7:42], (6, 6))
	dϕ = A*ϕ
	du[7:42] = reshape(dϕ, (:, 1))
end

# ╔═╡ fa642306-c5d0-461a-8d14-0b085a8147e2
begin
	u = [1.17, 0, 0, 0, -0.489780292125578, 0]
	μ = 0.0121505856
	tspan = (0.0, 25.0)
	append!(u, [1.0*Matrix(I, 6, 6)...])
	prob = ODEProblem(CR3BP!, u, tspan, μ)
	
	
# 	# FUNCTION TO SAVE EVENTS TO OUTSIDE VARIABLE
# 	func!(integrator) = push!(x, [integrator.t, integrator.u[1], integrator.u[4]])
	
# 	# CREATING CALLBACK
# 	cb = ContinuousCallback((u,t,integrator) -> u[2], func!)
	
	sol = solve(prob, reltol=1e-8, abstol=1e-8)
	# plot(sol, vars=(1,2))
	# Φ = Any[]
	# FTLE = Any[]
	# for i = 1:length(sol)
	# 	push!(Φ, reshape(sol.u[i][7:42], (6, 6)))
	# 	Δ = Φ[i]*transpose(Φ[i])
	# 	λ = eigvals(Δ)
	# 	push!(FTLE, (1/abs(sol.t[i]))* log(max(λ...)))
	# end
	# plot(sol.t, FTLE)
end

# ╔═╡ Cell order:
# ╟─233ec700-ccd1-11eb-2983-914a848c6aab
# ╠═6b72128e-e2bb-47a5-8a80-315d70348923
# ╟─366664b7-6b31-4ae7-b479-58a53a717ffb
# ╠═cacfe7ff-5939-4be7-9a28-08829d9e7e11
# ╠═fa642306-c5d0-461a-8d14-0b085a8147e2
