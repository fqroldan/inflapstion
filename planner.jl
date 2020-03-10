using Distributed

using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD

include("type_def.jl")
include("reporting_routines.jl")

function value(θv, x, ystar, γ, β, κ, itp_v)
	πv, yv, θp = x

	return (yv-ystar)^2 + γ*πv^2 - θp * (πv - κ*yv) + θv*πv + β*itp_v(θp)
end

function optim_step(rp::Ramsey, itp_v; choose::Bool=true)
	g, v = zeros(size(rp.g)), zeros(size(rp.v))

	ystar, γ, β, κ = rp.ystar, rp.γ, rp.β, rp.κ
	minπ = -0.1 * Nash(rp)
	miny = 1/κ * (minπ - β*Nash(rp))
	maxπ = 1.1*Nash(rp)
	maxy = 1/κ * (maxπ - β*minπ)

	minθ = 0.0
	maxθ = maximum(rp.θgrid)

	maxθ = maxθ - 0.01 * (maxθ-minθ)
	minθ = 1e-5
	
	minx = [minπ, miny, minθ]
	maxx = [maxπ, maxy, maxθ]

	for (jθ, θv) in enumerate(rp.θgrid)

		xguess = rp.g[jθ,:]
		xguess[3] = (maxθ + minθ)/2

		obj_f(x; choose::Bool=true, get_others=false) = value(θv, x, ystar, γ, β, κ, itp_v)
		if choose
			res = Optim.optimize(x -> obj_f(x), minx, maxx, xguess, Fminbox(NelderMead()))

			G = res.minimizer
			vv = res.minimum
		else
			G = xguess
			vv = obj_f(xguess, choose=false)
		end

		g[jθ,:] = G
		v[jθ] = vv

	end

	return g, v
end

function vfi_iter(rp::Ramsey; choose::Bool=true)
	knots = (rp.θgrid,)
	itp_v = interpolate(knots, rp.v, Gridded(Linear()))

	new_g, new_v = optim_step(rp, itp_v; choose=choose)

	return new_g, new_v
end

function vfi!(rp::Ramsey; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, upd_η = 0.75)

	iter = 0
	dist = 1e8

	while dist > tol && iter < maxiter
		iter += 1

		old_g, old_v = copy(rp.g), copy(rp.v)

		new_g, new_v = vfi_iter(rp, choose=true)

		dist_g = sqrt.(sum( (new_g  - old_g ).^2 )) / sqrt.(sum(old_g .^2))
		dist_v = sqrt.(sum( (new_v  - old_v ).^2 )) / sqrt.(sum(old_v .^2))

		dist = max(dist_g, dist_v)

		rp.v = old_v + upd_η * (new_v - old_v)
		rp.g = new_g

		if verbose && iter % 50 == 0
			print("\nAfter $iter iterations, d(v) = $(@sprintf("%0.3g",dist_v))")
		end

	end

	nothing
end




