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

function optim_step(rp::Ramsey, itp_v, itp_gπ)
	g, v = zeros(size(rp.g)), zeros(size(rp.v))

	ystar, γ, β, κ = rp.ystar, rp.γ, rp.β, rp.κ
	minπ = -0.1 * Nash(rp)
	miny = 1/κ * (minπ - β*Nash(rp))
	maxπ = 1.1*Nash(rp)
	maxy = 1/κ * (maxπ - β*minπ)

	minθ = minimum(rp.θgrid)
	maxθ = maximum(rp.θgrid)
	
	minx = [minπ, miny, minθ]
	maxx = [maxπ, maxy, maxθ]
	for jθ in 1:length(rp.θgrid)
		θv = rp.θgrid[jθ]
		xguess = rp.g[jθ,:]

		res = Optim.optimize(θp -> allFOCs(θv, θp, ystar, γ, β, κ, itp_gπ)^2, minθ, maxθ, GoldenSection())

		θp = res.minimizer
		if abs(res.minimum) > 1e-8
			println(res.minimum)
		end
		πv, yv = get_vars(θp, θv, κ, γ, ystar)

		G = [πv, yv, θp]

		g[jθ,:] = G
		v[jθ] = value(θv, G, ystar, γ, β, κ, itp_v)

	end

	return g, v
end

function value_dev(πv, sp::Sustainable)
	ystar, γ, β, κ, b, ξ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ
	return value_dev(πv, ystar, γ, β, κ, b, ξ)
end

function value_dev(πv, ystar, γ, β, κ, b, ξ)
	yv = (πv - β * ξ) / κ
	return (yv-ystar)^2 + γ*πv^2 + β*b
end

function value_comply(sp::Sustainable{T}, ap, av, ystar, γ, β, κ, itp_v, itp_gπ) where T <: PhillipsCurve
	πv = av
	if T == Fwd_strategy
		Eπ = itp_gπ(ap)
	elseif T == Fwd_literal
		Eπ = ap
	end
	yv = (πv - β * Eπ) / κ

	return (yv-ystar)^2 + γ*πv^2 + β * itp_v(ap)
end

function optim_step(sp::Sustainable, itp_v, itp_gπ)
	g, v = zeros(size(sp.g)), zeros(size(sp.v))

	ystar, γ, β, κ, b, ξ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ
	minπ = -0.1 * Nash(sp)
	maxπ = 1.1*Nash(sp)

	mina = minimum(sp.agrid)
	maxa = maximum(sp.agrid)

	# First choose best deviation
	obj_f(πv) = value_dev(πv, sp)
	res = Optim.optimize(obj_f, minπ, maxπ, GoldenSection())

	πd = res.minimizer
	vd = res.minimum

	# Now optimize given a
	for ja in 1:length(sp.agrid)
		av = sp.agrid[ja]

		# If complying
		obj_f2(ap) = value_comply(sp, ap, av, ystar, γ, β, κ, itp_v, itp_gπ)
		res = Optim.optimize(obj_f2, mina, maxa, GoldenSection())

		ap = res.minimizer
		vc = res.minimum

		if vc > vd
			yd = (πd - β * ξ) / κ
			g[ja, :] .= πd, yd, av
			v[ja] = vd
		else
			yc = (av - β * ap) / κ
			g[ja, :] .= av, yc, ap
			v[ja] = vc
		end
	end

	return g, v
end

function vfi_iter(pp::Union{Ramsey, Sustainable})
	itp_v = make_itp(pp, pp.v)
	itp_gπ = make_itp(pp, pp.g[:,1])

	new_g, new_v = optim_step(pp, itp_v, itp_gπ)

	return new_g, new_v
end

function vfi!(pp::Union{Ramsey, Sustainable}; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, upd_η = 0.75)

	iter = 0
	dist = 1e8

	while dist > tol && iter < maxiter
		iter += 1

		old_g, old_v = copy(pp.g), copy(pp.v)

		new_g, new_v = vfi_iter(pp)

		dist_g = sqrt.(sum( (new_g  - old_g ).^2 )) / sqrt.(sum(old_g .^2))
		dist_v = sqrt.(sum( (new_v  - old_v ).^2 )) / sqrt.(sum(old_v .^2))

		dist = max(dist_g, dist_v)

		pp.v = old_v + upd_η * (new_v - old_v)
		pp.g = new_g

		if verbose# && iter % 50 == 0
			print("\nAfter $iter iterations, d(v) = $(@sprintf("%0.3g",dist_v))")
		end

	end

	nothing
end

initial_state(rp::Ramsey) = 0

function initial_state(sp::Sustainable)

	Lmin, ja = findmin(sp.v)

	return sp.agrid[ja]
end

function best_response_ξ(sp::Sustainable)
	# yξ = (pp.ystar - pp.β*pp.κ*pp.γ*pp.ξ) / (1+pp.κ^2*pp.γ)
	# πξ = pp.κ*yξ + pp.β*pp.ξ
	minπ = -0.1 * Nash(sp)
	maxπ = 1.1*Nash(sp)
	obj_f(πv) = value_dev(πv, sp)
	res = Optim.optimize(obj_f, minπ, maxπ, GoldenSection())

	πd = res.minimizer
	vd = res.minimum
	return πd
end

policy(pp::Ramsey, θ, itp_gπ) = itp_gπ(θ)
new_state(pp::Ramsey, θ, π, itp_gθ) = itp_gθ(θ)

policy(pp::Sustainable, θ, itp_gπ) = ifelse(θ == -1.0, best_response_ξ(pp), itp_gπ(θ))
new_state(pp::Sustainable, θ, π, itp_gθ) = ifelse(θ != -1 || isapprox(π,θ), itp_gθ(θ), -1.0)

function simul_plan(pp::Union{Ramsey, Sustainable}, T = 4*10)

	θv = zeros(T)
	πv = zeros(T)
	# knots = (rp.θgrid,)
	# itp_gπ = interpolate(knots, rp.g[:,1], Gridded(Linear()))
	# itp_gθ = interpolate(knots, rp.g[:,3], Gridded(Linear()))
	itp_gπ = make_itp(pp, pp.g[:,1])
	itp_gθ = make_itp(pp, pp.g[:,3])

	θt = initial_state(pp)
	for jt in 1:T
		πt = policy(pp, θt, itp_gπ)

		θv[jt] = θt
		πv[jt] = πt

		θt = new_state(pp, θt, πt, itp_gθ)
	end

	return πv, θv
end


