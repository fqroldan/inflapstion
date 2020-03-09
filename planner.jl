using Distributed

using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD

include("type_def.jl")
include("reporting_routines.jl")

function value(x, θv, θp, ystar, γ, β, κ, itp_vf)
	πv, yv = x
	return (yv-ystar)^2 + γ*πv^2 + θp * (πv - κ*yv) - θv*πv + β*itp_vf(θp)
end

function opt_value(θv, θp, minx, maxx, xguess, ystar, γ, β, κ, itp_vf; choose::Bool=true, get_others::Bool=false)

	if choose
		println(xguess)
		println(θp)
		res = Optim.optimize(
			x -> value(x, θv, θp, ystar, γ, β, κ, itp_vf),
			minx, maxx, xguess, Fminbox(LBFGS()))

		vf = res.minimum
		if get_others
			πv, yv = res.minimizer
			return πv, yv
		end
	else
		vf = value(xguess, θv, θp, ystar, γ, β, κ, itp_vf)
	end

	return vf
end

function optim_step(rp::Ramsey, itp_vf; choose::Bool=true)
	gπ, gy, vf = zeros(size(rp.gπ)), zeros(size(rp.gy)), zeros(size(rp.vf))

	ystar, γ, β, κ = rp.ystar, rp.γ, rp.β, rp.κ
	minπ = -0.1 * Nash(rp)
	miny = 1/κ * (minπ - β*Nash(rp))
	maxπ = 1.1*Nash(rp)
	maxy = 1/κ * (maxπ - β*minπ)

	minx = [minπ, miny]
	maxx = [maxπ, maxy]
	
	minθ = 0.0
	maxθ = maximum(rp.θgrid)

	maxθ = maxθ - 0.01 * (maxθ-minθ)

	for (jθ, θv) in enumerate(rp.θgrid)
		πg = rp.gπ[jθ]
		yg = rp.gy[jθ]
		xguess = [πg, yg]

		θg = (maxθ + minθ)/2

		obj_f(x; choose::Bool=true, get_others=false) = opt_value(θv, x, minx, maxx, xguess, ystar, γ, β, κ, itp_vf, choose=choose, get_others=get_others)
		if choose
			res = Optim.optimize(x -> -obj_f(x), minθ, maxθ, GoldenSection())

			θp = res.minimizer
			vv = -res.minimum

			πv, yv = obj_f(θp, get_others=true)
		else
			πv, yv = xguess
			vv = obj_f(guess, choose=false)
		end

		gπ[jθ] = πv
		gy[jθ] = yv
		vf[jθ] = vv

	end

	return gπ, gy, vf
end

function vfi_iter(rp::Ramsey; choose::Bool=true)
	knots = (rp.θgrid,)
	itp_vf = interpolate(knots, rp.vf, Gridded(Linear()))
	# itp_vf = extrapolate(itp_vf, Interpolations.Line())

	new_gπ, new_gy, new_vf = optim_step(rp, itp_vf; choose=choose)

	return new_gπ, new_gy, new_vf
end

function vfi!(rp::Ramsey; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, upd_η = 1)

	iter = 0
	dist = 1e8

	while dist > tol && iter < maxiter
		iter += 1

		for jj in 1:4
			# _, _, _ = vfi_iter(rp; choose=false)
		end

		old_gπ, old_gy, old_vf = copy(rp.gπ), copy(rp.gy), copy(rp.vf)

		new_gπ, new_gy, new_vf = vfi_iter(rp, choose=true)

		dist_gπ = sqrt.(sum( (new_gπ  - old_gπ ).^2 )) / sqrt.(sum(old_gπ .^2))
		dist_gy = sqrt.(sum( (new_gy  - old_gy ).^2 )) / sqrt.(sum(old_gy .^2))
		dist_vf = sqrt.(sum( (new_vf  - old_vf ).^2 )) / sqrt.(sum(old_vf .^2))

		dist = max(dist_gπ, dist_gy, dist_vf)

		rp.vf = old_vf + upd_η * (new_vf - old_vf)
		rp.gπ = new_gπ
		rp.gy = new_gy

		if verbose
			print("\nAfter $iter iterations, d(v) = $(@sprintf("%0.3g",dist_vf))")
		end

	end

	nothing
end




