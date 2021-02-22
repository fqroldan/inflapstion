using Distributed

using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD2, FileIO

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

function πbounds(sp::Sustainable)
	minπ = -0.1 * Nash(sp)
	maxπ =  1.1 * Nash(sp)
	return minπ, maxπ
end

function abounds(sp::Sustainable)
	mina = minimum(sp.agrid)
	maxa = maximum(sp.agrid)
	return mina, maxa
end

function value_dev(πv, sp::Sustainable, itp_v)
	ystar, γ, β, κ, b, ξ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ
	return value_dev(πv, ystar, γ, β, κ, b, ξ)
end

function value_dev(πv, ystar, γ, β, κ, b, ξ)
	yv = (πv - β * ξ) / κ
	return (yv-ystar)^2 + γ*πv^2 + β*b
end

function value_dev(πv, ap, sp::Sustainable{Fwd_GP}, itp_v)
	ystar, γ, β, κ, b, ξ, θ, D, σ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ, sp.θ, sp.D, sp.σ
	mina, maxa = abounds(sp)
	if ap < mina || ap > maxa
		itp_v  = extrapolate(itp_v,  Interpolations.Flat())
	end

	yv = (πv - β * ξ) / κ

	vB = (yv - ystar)^2 + γ * πv^2 + β * (θ*itp_v(ap) + (1-θ)*b) + σ^2 * (γ + 1/κ^2)
	return vB
end

function find_best_dev(sp::Sustainable, itp_v)
	minπ, maxπ = πbounds(sp)

	obj_f(πv) = value_dev(πv, sp, itp_v)
	res = Optim.optimize(obj_f, minπ, maxπ, GoldenSection())

	πd = res.minimizer
	vd = res.minimum

	return vd, πd
end

function find_best_dev(sp::Sustainable{Fwd_GP}, itp_v)
	minπ, maxπ = πbounds(sp)
	mina, maxa = abounds(sp)
	xguess = [(minπ + maxπ)/2, (mina+maxa)/2]

	obj_f(x) = value_dev(x[1], x[2], sp, itp_v)
	res = Optim.optimize(obj_f, [minπ, mina], [maxπ, maxa], xguess)

	πd, ap = res.minimizer
	vd = res.minimum
	sp.b = vd
	return vd, πd
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

function find_best_comp(sp::Sustainable, av, itp_v, itp_gπ)
	mina, maxa = abounds(sp)
	ystar, γ, β, κ, b, ξ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ

	obj_f2(ap) = value_comply(sp, ap, av, ystar, γ, β, κ, itp_v, itp_gπ)
	res = Optim.optimize(obj_f2, mina, maxa, GoldenSection())

	ap = res.minimizer
	vc = res.minimum
	return ap, av, vc
end

function cond_value(πv, ap, sp::Sustainable{Fwd_GP}, av, itp_v, itp_gπ)
	ystar, γ, β, κ, b, ξ, θ, D, σ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ, sp.θ, sp.D, sp.σ
	if abs(πv-av) > av*D
	# if abs(πv-av) > D
		pp = 0.0
	else
		pp = 1.0
	end
	mina, maxa = abounds(sp)
	if ap < mina || ap > maxa
		itp_gπ = extrapolate(itp_gπ, Interpolations.Flat())
		itp_v  = extrapolate(itp_v,  Interpolations.Flat())
	end
	yv = πv - β * ( pp * itp_gπ(ap) + (1-pp) * ξ )

	vG = (yv-ystar)^2 + γ * πv^2 + β * ( pp * itp_v(ap) + (1-pp) * b)
	return vG
end

function exp_value(πv, ap, sp::Sustainable{Fwd_GP}, av, itp_v, itp_gπ)
	pdf_ϵ(sp, ϵv) = pdf(Normal(0, sp.σ), ϵv)
	f(ϵv) = cond_value(πv + ϵv, ap, sp, av, itp_v, itp_gπ) * pdf_ϵ(sp, ϵv)
	(val, err) = hquadrature(f, -3.09*sp.σ, 3.09*sp.σ, rtol=1e-10, atol=0, maxevals=0)

	sum_prob = cdf_ϵ(sp, 3.09*sp.σ) - cdf_ϵ(sp, -3.09*sp.σ)
	val = val / sum_prob
	return val
end

function find_best_comp(sp::Sustainable{Fwd_GP}, av, itp_v, itp_gπ)
	minπ, maxπ = πbounds(sp)
	ttπ = 0.01 * (maxπ - minπ)
	mina, maxa = abounds(sp)
	tta = 0.01 * (maxa - mina)
	xguess = [min(max(itp_gπ(av), minπ+ttπ), maxπ-ttπ), min(max(av, mina+tta), maxa-tta)]

	obj_f(x) = exp_value(x[1], x[2], sp, av, itp_v, itp_gπ)
	res = Optim.optimize(obj_f, [minπ, mina], [maxπ, maxa], xguess, Fminbox(NelderMead()))

	πv, ap = res.minimizer
	vv = res.minimum
	return ap, πv, vv
end


function update_comp_dev(sp::Sustainable, vc, vd, πd, πc, av, ap)
	ystar, γ, β, κ, b, ξ = sp.ystar, sp.γ, sp.β, sp.κ, sp.b, sp.ξ
	if vc > vd
		yd = (πd - β * ξ) / κ
		gg = [πd, yd, av]
		vv = vd
	else
		yc = (πc - β * ap) / κ
		gg = [πc, yc, ap]
		vv = vc
	end

	return vv, gg
end

function optim_step(sp::Sustainable, itp_v, itp_gπ)
	g, v = zeros(size(sp.g)), zeros(size(sp.v))

	# First choose best deviation
	vd, πd = find_best_dev(sp, itp_v)

	# Now optimize given a

	# for ja in 1:length(sp.agrid)
	Threads.@threads for ja in 1:length(sp.agrid)
		av = sp.agrid[ja]

		ap, πc, vc = find_best_comp(sp, av, itp_v, itp_gπ)

		vv, gg = update_comp_dev(sp, vc, vd, πd, πc, av, ap)

		v[ja] = vv
		g[ja,:] = gg

	end

	return g, v
end

function vfi_iter(pp::Union{Ramsey, Sustainable})
	itp_v = make_itp(pp, pp.v)
	itp_gπ = make_itp(pp, pp.g[:,1])

	new_g, new_v = optim_step(pp, itp_v, itp_gπ)

	return new_g, new_v
end

function vfi!(pp::Union{Ramsey, Sustainable}; tol::Float64=15e-4, maxiter::Int64=2500, verbose::Bool=true, upd_η = 0.75)

	cens_v = 1e-4
	if typeof(pp) == Sustainable{Fwd_GP}
		upd_η = 0.5
		cens_v = 1e-2
	end

	iter = 0
	dist = 10+tol

	while dist > tol && iter < maxiter
		iter += 1

		old_g, old_v = copy(pp.g), copy(pp.v);

		new_g, new_v = vfi_iter(pp);

		norm_g = max(sqrt.(sum(old_g .^2))/length(old_g), cens_v);
		norm_v = max(sqrt.(sum(old_v .^2))/length(old_v), cens_v)

		dist_g = sqrt.(sum( (new_g  - old_g ).^2 ))/length(old_g) / norm_g;
		dist_v = sqrt.(sum( (new_v  - old_v ).^2 ))/length(old_v) / norm_v

		dist = max(dist_g, dist_v)

		pp.v = old_v + upd_η * (new_v - old_v);
		pp.g = old_g + upd_η * (new_g - old_g);

		if verbose && iter % 50 == 0
			print("After $iter iterations, d(v,g) = $(@sprintf("%0.3g",dist_v)) $(@sprintf("%0.3g",dist_g)) at $(Dates.format(now(),"HH:MM"))\n")
		end
		if typeof(pp)<:Sustainable && iter > 200
			upd_η = 0.25
			if iter == 400
				upd_η = 0.1
			elseif iter > 600
				upd_η = max(1e-5, upd_η * 0.975)
			end
		end
	end
	if verbose
		final_report(pp, iter, maxiter, dist)
	end
	return dist <= tol
end
final_report(pp::Sustainable, iter, maxiter, dist) = final_report(iter,maxiter,dist)

function show_value(rp::Ramsey)
	knots = (rp.θgrid,)
	itp = interpolate(knots, rp.v, Gridded(Linear()))
	return itp(0.0)
end
function final_report(rp::Ramsey, iter, maxiter, dist)
	vR = show_value(rp)
	print("value attained = $(@sprintf("%0.3g",vR))\n")
	final_report(iter, maxiter, dist)
end

function final_report(iter, maxiter, dist)
	if iter < maxiter
		print("Converged in $iter iterations at $(Dates.format(now(),"HH:MM"))\n")
	else
		print("Failed to converge\n")
	end
end

initial_state(rp::Ramsey) = 0

function initial_state(sp::Sustainable)

	Lmin, ja = findmin(sp.v)

	return sp.agrid[ja]
end

function best_response_ξ(sp::Sustainable, itp_v)
	# yξ = (pp.ystar - pp.β*pp.κ*pp.γ*pp.ξ) / (1+pp.κ^2*pp.γ)
	# πξ = pp.κ*yξ + pp.β*pp.ξ
	minπ = -0.1 * Nash(sp)
	maxπ = 1.1*Nash(sp)
	obj_f(πv) = value_dev(πv, sp, itp_v)
	res = Optim.optimize(obj_f, minπ, maxπ, GoldenSection())

	πd = res.minimizer
	vd = res.minimum
	return πd
end

policy(pp::Ramsey, θ, itp_gπ, itp_v) = itp_gπ(θ)
new_state(pp::Ramsey, θ, π, itp_gθ) = itp_gθ(θ)

policy(pp::Sustainable, θ, itp_gπ, itp_v) = ifelse(θ == -1.0, best_response_ξ(pp, itp_v), itp_gπ(θ))
new_state(pp::Sustainable, θ, π, itp_gθ) = ifelse(θ != -1 || isapprox(π,θ), itp_gθ(θ), -1.0)

function simul_plan(pp::Union{Ramsey, Sustainable}, T = 4*10)

	θv = zeros(T)
	πv = zeros(T)
	# knots = (rp.θgrid,)
	# itp_gπ = interpolate(knots, rp.g[:,1], Gridded(Linear()))
	# itp_gθ = interpolate(knots, rp.g[:,3], Gridded(Linear()))
	itp_v = make_itp(pp, pp.v)
	itp_gπ = make_itp(pp, pp.g[:,1])
	itp_gθ = make_itp(pp, pp.g[:,3])

	θt = initial_state(pp)
	for jt in 1:T
		πt = policy(pp, θt, itp_gπ, itp_v)

		θv[jt] = θt
		πv[jt] = πt

		θt = new_state(pp, θt, πt, itp_gθ)
	end

	return πv, θv
end


