function value(θv, x, ystar, γ, β, κ, itp_v)
	πv, yv, θp = x
	if θp < 0 || θp > 200
		itp_v = extrapolate(itp_v, Interpolations.Flat())
	end

	return (yv-ystar)^2 + γ*πv^2 + θp * (πv - κ*yv) - θv*πv + β*itp_v(θp)
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

	ystar, γ, β, κ = (rp.pars[k] for k in (:ystar, :γ, :β, :κ))
	minπ = -0.1 * Nash(rp)
	miny = 1/κ * (minπ - β*Nash(rp))
	maxπ = 1.1*Nash(rp)
	maxy = 1/κ * (maxπ - β*minπ)

	minθ = minimum(rp.θgrid)
	maxθ = maximum(rp.θgrid)
	
	# minx = [minπ, miny, minθ]
	# maxx = [maxπ, maxy, maxθ]
	for jθ in 1:length(rp.θgrid)
		θv = rp.θgrid[jθ]
		# xguess = rp.g[jθ,:]

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

function vfi_iter(pp::Ramsey)
	itp_v = make_itp(pp, pp.v)
	itp_gπ = make_itp(pp, pp.g[:,1])

	new_g, new_v = optim_step(pp, itp_v, itp_gπ)

	return new_g, new_v
end

function vfi!(pp::Ramsey; tol::Float64=1e-8, maxiter::Int64=2500, verbose::Bool=true, upd_η=1)
    iter = 0
    dist = 10 + tol

    while dist > tol && iter < maxiter
        iter += 1

        old_g, old_v = copy(pp.g), copy(pp.v)

        new_g, new_v = vfi_iter(pp)

        norm_g = max(norm(pp.g), 1)
        norm_v = max(norm(pp.v), 1)

        dist_g = norm(new_g - pp.g) / norm_g
        dist_v = norm(new_v - pp.v) / norm_v

        dist = max(dist_g, dist_v)

        pp.v = pp.v + upd_η * (new_v - pp.v)
        pp.g = pp.g + upd_η * (new_g - pp.g)

        if verbose
            print("After $iter iterations, d(v,g) = $(@sprintf("%0.3g",dist_v)) $(@sprintf("%0.3g",dist_g)) at $(Dates.format(now(),"HH:MM"))\n")
        end
    end
    if verbose
        final_report(pp, iter, maxiter, dist)
    end
    return dist <= tol
end
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

policy(::Ramsey, θ, itp_gπ, _) = itp_gπ(θ)
new_state(::Ramsey, θ, _, itp_gθ) = itp_gθ(θ)

function make_itp(rp::Ramsey, y)
    knots = (rp.θgrid,)
    itp = interpolate(knots, y, Gridded(Linear()))
    return itp
end

function simul_plan(pp::Ramsey, T=4 * 10)

    θv = zeros(T)
    πv = zeros(T)

    itp_v = make_itp(pp, pp.v)
    itp_gπ = make_itp(pp, pp.g[:, 1])
    itp_gθ = make_itp(pp, pp.g[:, 3])

    θt = initial_state(pp)
    for jt in 1:T
        πt = policy(pp, θt, itp_gπ, itp_v)

        θv[jt] = θt
        πv[jt] = πt

        θt = new_state(pp, θt, πt, itp_gθ)
    end

    return πv, θv
end

function checkRamsey(mt::MultiType)

    β, γ, κ = (mt.ct.pars[k] for k in (:β, :γ, :κ))
    λ = γ * κ^2

    return checkRamsey(β, λ)
end

checkRamsey(β, λ) = (1 + β + λ - sqrt((1+β+λ)^2 - 4*β)) / (2*β)
