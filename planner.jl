using Distributed

using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD

include("type_def.jl")
include("reporting_routines.jl")

function value(θv, x, ystar, γ, β, κ, itp_v)
	πv, yv, θp = x
	if θp < 0 || θp > 200
		itp_v = extrapolate(itp_v, Interpolations.Flat())
	end

	return (yv-ystar)^2 + γ*πv^2 + θp * (πv - κ*yv) - θv*πv + β*itp_v(θp)
end

function FOCs(θv, x, ystar, γ, β, κ, itp_gπ)
	πv, yv, θp = x

	F = zeros(3)
	if θp < 0 || θp > 200
		itp_gπ = extrapolate(itp_gπ, Interpolations.Flat())
	end

	F[1] = yv - ystar - κ * θp
	F[2] = γ * πv + θp - θv
	F[3] = πv - κ*yv - β * itp_gπ(θp)

	return F'F
end

get_π(θp, θv, γ) = (θv - θp) / γ
get_y(θp, κ, ystar) = ystar + κ*θp
get_vars(θp, θv, κ, γ, ystar) = get_π(θp, θv, γ), get_y(θp, κ, ystar)

function allFOCs(θv, θp, ystar, γ, β, κ, itp_gπ)

	πv, yv = get_vars(θp, θv, κ, γ, ystar)

	return πv - κ*yv - β * itp_gπ(θp)
end

function optim_step(rp::Ramsey, itp_v, itp_gπ; choose::Bool=true)
	g, v = zeros(size(rp.g)), zeros(size(rp.v))

	ystar, γ, β, κ = rp.ystar, rp.γ, rp.β, rp.κ
	minπ = -0.1 * Nash(rp)
	miny = 1/κ * (minπ - β*Nash(rp))
	maxπ = 1.1*Nash(rp)
	maxy = 1/κ * (maxπ - β*minπ)

	minθ = minimum(rp.θgrid)
	maxθ = maximum(rp.θgrid)

	# maxθ = maxθ - 0.01 * (maxθ-minθ)
	# minθ = 0.0
	
	minx = [minπ, miny, minθ]
	maxx = [maxπ, maxy, maxθ]
	# Threads.@threads for jθ in 1:length(rp.θgrid)
	for jθ in 1:length(rp.θgrid)
		θv = rp.θgrid[jθ]
		xguess = rp.g[jθ,:]
		# xguess[3] = (maxθ + minθ)/2


		res = Optim.optimize(θp -> allFOCs(θv, θp, ystar, γ, β, κ, itp_gπ)^2, minθ, maxθ, GoldenSection())

		θp = res.minimizer
		if abs(res.minimum) > 1e-8
			println(res.minimum)
		end
		πv, yv = get_vars(θp, θv, κ, γ, ystar)

		G = [πv, yv, θp]
		# res = Optim.optimize(x -> FOCs(θv, x, ystar, γ, β, κ, itp_gπ), minx, maxx, xguess, Fminbox(NelderMead()))

		# obj_f(x; choose::Bool=true, get_others=false) = value(θv, x, ystar, γ, β, κ, itp_v)
		# res = Optim.optimize(x -> obj_f(x), minx, maxx, xguess, Fminbox(NelderMead()))

		# if Optim.converged(res)
		# else
		# 	println("caution")
		# end

		# G = res.minimizer
		# vv = res.minimum
		# if vv > 1e-8
		# 	println(vv)
		# end

		g[jθ,:] = G
		v[jθ] = value(θv, G, ystar, γ, β, κ, itp_v)

	end

	return g, v
end

function vfi_iter(rp::Ramsey)
	knots = (rp.θgrid,)
	itp_v = interpolate(knots, rp.v, Gridded(Linear()))
	itp_gπ = interpolate(knots, rp.g[:,1], Gridded(Linear()))

	new_g, new_v = optim_step(rp, itp_v, itp_gπ)

	return new_g, new_v
end

function vfi!(rp::Ramsey; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, upd_η = 0.75)

	iter = 0
	dist = 1e8

	while dist > tol && iter < maxiter
		iter += 1

		old_g, old_v = copy(rp.g), copy(rp.v)

		new_g, new_v = vfi_iter(rp)

		dist_g = sqrt.(sum( (new_g  - old_g ).^2 )) / sqrt.(sum(old_g .^2))
		dist_v = sqrt.(sum( (new_v  - old_v ).^2 )) / sqrt.(sum(old_v .^2))

		dist = max(dist_g, dist_v)

		rp.v = old_v + upd_η * (new_v - old_v)
		rp.g = new_g

		if verbose# && iter % 50 == 0
			print("\nAfter $iter iterations, d(v) = $(@sprintf("%0.3g",dist_v))")
		end

	end

	nothing
end

function simul_plan(rp::Ramsey, T = 4*10)

	θv = zeros(T)
	πv = zeros(T)
	knots = (rp.θgrid,)
	itp_gπ = interpolate(knots, rp.g[:,1], Gridded(Linear()))
	itp_gθ = interpolate(knots, rp.g[:,3], Gridded(Linear()))

	θt = 0
	for jt in 1:T
		πt = itp_gπ(θt)

		θv[jt] = θt
		πv[jt] = πt

		θt = itp_gθ(θt)
	end

	return πv, θv
end