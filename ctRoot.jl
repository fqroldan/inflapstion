using Distributions, Interpolations, Optim, HCubature, QuantEcon, Printf, PlotlyJS, Dates, JLD2, LinearAlgebra

include("type_def.jl")
include("reporting_routines.jl")
include("simul.jl")
include("plotsct.jl")
include("prequel.jl")
include("planner.jl")

function Bayes(ct::Plan, obs_π, exp_π, pv, av)
	
	if isapprox(pv, 0.0)
		p′ = 0.0
	elseif isapprox(pv, 1.0)
		p′ = 1.0
	else
		numer = pv * pdf_ϵ(ct, obs_π - av)
		denomin = numer + (1.0-pv) * pdf_ϵ(ct, obs_π - exp_π)
		if isapprox(denomin, 0.0)
			p′ = 0.0
		else
			p′ = numer / denomin
		end
	end

	return p′
end

function next_a(ct::Plan, av, apv, π)
	if haskey(ct.pars, :ψ)
		return ϕ(ct, av) + ct.pars[:ψ] * (π - av)
	else
		return ϕ(ct, av)
	end
end

next_a(ct::DovisKirpalani, av, apv, π) = apv

function cond_L(ct::Plan{T}, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime, ge, πe′; use_ϕ = true) where {T<:PhillipsCurve}
	ystar, γ, β = (ct.pars[k] for k in (:ystar, :γ, :β))

    pprime = Bayes(ct, obs_π, ge, pv, av)

    πe = pv * av + (1 - pv) * ge

    a_min, a_max = extrema(ct.gr[:a])

    π_today = max(a_min, min(a_max, obs_π))
	if use_ϕ
    	aprime = next_a(ct, av, aprime, π_today)
    	# aprime = next_a(ct, av, aprime, πe)
	end
    aprime = max(a_min, min(a_max, aprime))

    # if aprime <= minimum(ct.agrid) || aprime >= maximum(ct.agrid)
    itp_L = extrapolate(itp_L, Interpolations.Line())
    itp_gπ = extrapolate(itp_gπ, Interpolations.Flat())
    itp_C = extrapolate(itp_C, Interpolations.Line())
    # end

    L′ = itp_L(pprime, aprime)
    exp_π′ = 0.0
    if T == Forward
        exp_π′ = pprime * aprime + (1.0 - pprime) * itp_gπ(pprime, aprime)
    end

    y = PC(ct, obs_π, πe, exp_π′, πe′) # Automatically uses method for forward or backward
    L = (ystar - y)^2 + γ * obs_π^2 + β * L′
    C′ = itp_C(pprime, aprime)

    return L, pprime, y, C′
end

get_sumprob(ct::Plan) = cdf_ϵ(ct, 3.09*ct.pars[:σ]) - cdf_ϵ(ct, -3.09*ct.pars[:σ])

function exp_L_y(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, aprime, ge, πe′; use_ϕ = true)

	sum_prob = get_sumprob(ct)

	f_C(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime, ge, πe′, use_ϕ = use_ϕ)[4] * pdf_ϵ(ct, ϵv)
	Ec, err = hquadrature(f_C, -3.09*ct.pars[:σ], 3.09*ct.pars[:σ], rtol=1e-10, atol=0, maxevals=0)

	Ec = Ec / sum_prob

	return Ec
end

function exp_L(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, aprime, ge, πe′; use_ϕ = true)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime, ge, πe′, use_ϕ = use_ϕ)[1] * pdf_ϵ(ct, ϵv)
	val, err = hquadrature(f, -3.09*ct.pars[:σ], 3.09*ct.pars[:σ], rtol=1e-10, atol=0, maxevals=0)
	sum_prob = get_sumprob(ct)

	return val/sum_prob
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′; use_ϕ = true)
    π_guess = xguess[1]
	aprime = xguess[2]

    gπ, L = π_guess, itp_L(pv, av)

    obj_f(x) = exp_L(ct, itp_gπ, itp_L, itp_C, first(x), pv, av, aprime, ge, πe′, use_ϕ = use_ϕ)
	
    πN = Nash(ct)
    h = 0.05 * πN
    minπ, maxπ = -h, πN + h

    res = Optim.optimize(obj_f, minπ, maxπ, GoldenSection())
    gπ, L = res.minimizer, res.minimum

    return gπ, L, aprime
end

function exp_π_prime(ct::Plan{T}, pv, av, itp_gπ, ge, aprime) where T<:PhillipsCurve
	return 0.0
end

function π_prime(ct, ϵv, pv, av, itp_gπ, ge, aprime)
	exp_π = pv*av + (1-pv)*ge
	pprime = Bayes(ct, ge+ϵv, exp_π, pv, av)

	ge′ = itp_gπ(pprime, aprime)

	return pv*aprime + (1-pv)*ge′
end

function exp_π_prime(ct::Plan{SemiForward}, pv, av, itp_gπ, ge, aprime)
	f(ϵv) = π_prime(ct, ϵv, pv, av, itp_gπ, ge, aprime) * pdf_ϵ(ct, ϵv)
	val::Float64, err::Float64 = hquadrature(f, -3.09*ct.pars[:σ], 3.09*ct.pars[:σ], rtol=1e-10, atol=0, maxevals=0)
	sum_prob = get_sumprob(ct)

	return val/sum_prob
end

function optim_step(ct::Plan, itp_gπ, itp_L, itp_C, optim::Bool = true)
	gπ, ga = zeros(size(ct.gπ)), zeros(size(ct.ga))
	L  	   = zeros(size(ct.L))
	πN 	   = Nash(ct)
	
	h = 0.05 * πN
	Gmin, Gmax = -h, πN + h
	
	apgrid = gridmake(1:N(ct, :p), 1:N(ct, :a))
	Threads.@threads for js in axes(apgrid,1)
		# for js in axes(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.gr[:p][jp], ct.gr[:a][ja]
		
		aprime = ct.ga[jp, ja]
		ge = ct.gπ[jp, ja]

		xguess = [ge, aprime]

		πe′ = exp_π_prime(ct, pv, av, itp_gπ, ge, aprime)

		if optim
			obj(G) = (G - opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, G, πe′)[1])^2

			res = Optim.optimize(obj, Gmin, Gmax, GoldenSection())

			if res.minimum < 1e-6
				Gc = res.minimizer
			else
				Gc = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′)[1]
			end
			
			_, Lc, ap = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, Gc, πe′)
		else
			Gc = ge
			Lc = exp_L(ct, itp_gπ, itp_L, itp_C, ge, pv, av, aprime, ge, πe′)
			ap = aprime
		end

		gπ[jp, ja] = Gc
		L[jp, ja] = Lc
		ga[jp, ja] = ap

	end
	return gπ, L
end

function pf_iter(ct::Plan; optim=true)
	knts = (ct.gr[:p], ct.gr[:a]);
	itp_gπ = interpolate(knts, ct.gπ, Gridded(Linear()));
	itp_L  = interpolate(knts, ct.L, Gridded(Linear()));
	itp_C  = interpolate(knts, ct.C, Gridded(Linear()));

	new_gπ, new_L = optim_step(ct, itp_gπ, itp_L, itp_C, optim)

	return new_gπ, new_L
end

# function update_others!(ct::Plan, new_others, upd_η2)
# 	new_y, new_π, new_p, new_C, new_a = new_others[:]
# 	ct.Ey = new_y
# 	ct.Eπ = new_π
# 	ct.Ep = new_p
# 	ct.C  = new_C
# 	ct.ga = ct.ga + upd_η2 * (new_a - ct.ga)
# 	nothing
# end

function iter_cred!(new_C, ct::Plan, itp_C, itp_gπ, itp_L)
	β = ct.pars[:β]
	πN = Nash(ct)

	apgrid = gridmake(1:N(ct,:p), 1:N(ct,:a))
	Threads.@threads for js in axes(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.gr[:p][jp], ct.gr[:a][ja]

		aprime = ct.ga[jp, ja]
		ge = ct.gπ[jp, ja]
		
		πe′ = exp_π_prime(ct, pv, av, itp_gπ, ge, aprime)
		
		EC′ = exp_L_y(ct, itp_gπ, itp_L, itp_C, ge, pv, av, aprime, ge, πe′)

		Eπ = pv * av + (1-pv) * ge

		if πN > av
			new_C[jp, ja] = (1-β)*(πN - Eπ)/(πN-av) + β * EC′
		else
			new_C[jp, ja] = (1-β) * 1 + β * EC′
		end
	end
end
		

function cred_vfi!(ct::Plan; tol = 1e-5, maxiter = 2000, verbose = false)
	dist, iter = 1+tol, 0

	knts = (ct.gr[:p], ct.gr[:a])
	itp_gπ = interpolate(knts, ct.gπ, Gridded(Linear()));
	itp_L  = interpolate(knts, ct.L, Gridded(Linear()));

	new_C = similar(ct.C)
	while dist > tol && iter < maxiter
		iter += 1
		verbose && print("iter $iter: ")

		itp_C = interpolate(knts, ct.C, Gridded(Linear()));

		iter_cred!(new_C, ct, itp_C, itp_gπ, itp_L)

		norm_C = max(1, norm(ct.C))
		dist = norm(new_C - ct.C) / norm_C

		verbose && print("d, n = $dist, $norm_C\n")
		ct.C .= new_C
	end
	verbose && print("Converged in $iter iterations.\n")
	norm(ct.C)
end

function pfi!(ct::Plan; miniter::Int = 2, tol::Float64=1e-5, maxiter::Int64=2000, verbose::Bool=true, accelerate = true)
	dist = 10.
	iter = 0
	upd_η = 1

	rep = "Starting PFI (tol = $(@sprintf("%0.3g",tol)))\n"
	verbose && print(rep)

	while iter < miniter || (dist > tol && iter < maxiter)
		iter += 1
		if accelerate && dist < min(1e-2, tol * 100)
			verbose && print("Acceleration step\n")
			for _ in 1:10
			_, new_L = pf_iter(ct, optim = false)
			ct.L  = upd_η * new_L  + (1.0-upd_η) * ct.L
			end
		end
		
		new_gπ, new_L = pf_iter(ct, optim = true)
		# update_others!(ct, new_others, upd_η)

		norm_L = max(1, norm(ct.L))
		norm_g = max(1, norm(ct.gπ))

		dist_L = norm(new_L  - ct.L ) / norm_L
		dist_g = norm(new_gπ - ct.gπ) / norm_g

		dist = max(dist_L, dist_g)

		ct.L  = upd_η * new_L  + (1.0-upd_η) * ct.L
		ct.gπ = upd_η * new_gπ + (1.0-upd_η) * ct.gπ

		verbose && print("Iter $iter: d(L, g) = ($(@sprintf("%0.3g",dist_L)), $(@sprintf("%0.3g",dist_g))) at |L| = $(@sprintf("%0.3g",norm_L))\n")
	end

	if verbose && dist <= tol
		print("Converged in $iter iterations.\n")
	elseif verbose
		print("After $iter iterations, d(L) = $(@sprintf("%0.3g",dist))\n")
	end

	cred_vfi!(ct)

	return (dist <= tol)
end

function solve_all!(mt::MultiType; verbose = true, check = false)
	verbose && print("Going over all plans at $(Dates.format(now(), "HH:MM"))\n")
	iter = 0
	tot  = length(mt.ωgrid) * length(mt.χgrid)
	for (jω, ωv) in enumerate(mt.ωgrid), (jχ, χv) in enumerate(mt.χgrid)
		iter += 1
		
		show_ω = @sprintf("%.3g", perc_rate(ωv))
		show_χ = @sprintf("%.3g", annualized(χv))
		
		verbose && print("Plan with (ω, χ) = ($show_ω%, $show_χ%)")
		ct = mt.ct

		if jω > 1 && jχ == 1
			ct.L .= mt.L_mat[jω - 1, jχ, :, :]
		end

		ct.pars[:ω] = ωv
		ct.pars[:χ] = χv
		update_ga!(ct)
		if check
			ct.L .= mt.L_mat[jω, jχ, :, :]
			ct.gπ .= mt.g_mat[jω, jχ, :, :]
		end

		flag = pfi!(ct, verbose = false)
		verbose && flag && print(": ✓")
		verbose && !flag && print(": no convergence.")
		
		mt.L_mat[jω, jχ, :, :] .= ct.L
		mt.C_mat[jω, jχ, :, :] .= ct.C
		mt.g_mat[jω, jχ, :, :] .= ct.gπ
		
		perc = 100 * iter / tot
		verbose && print(" $(@sprintf("%.3g",perc))% completed.\n")
	end
end

function solve_all!(mt::Multiψ; verbose = true, check = false, save_progress = false)
	verbose && print("Going over all plans at $(Dates.format(now(), "HH:MM"))\n")
	iter = 0
	tot  = length(mt.ωgrid) * length(mt.χgrid) * length(mt.ψgrid)
	for (jω, ωv) in enumerate(mt.ωgrid), (jχ, χv) in enumerate(mt.χgrid), (jψ, ψv) in enumerate(mt.ψgrid)
		iter += 1
		
		show_ω = @sprintf("%.3g", ωv)
		show_χ = @sprintf("%.3g", annualized(χv))
		show_ψ = @sprintf("%.3g", 100*ψv)
		
		verbose && print("Plan with (ω, χ, ψ) = ($show_ω, $show_χ%, $show_ψ%)")
		ct = mt.ct

		if jω > 1 && jχ == 1
			ct.L .= mt.L[jω - 1, jχ, jψ, :, :]
		end

		ct.pars[:ω] = ωv
		ct.pars[:χ] = χv
		ct.pars[:ψ] = ψv
		update_ga!(ct)
		if check
			ct.L .= mt.L[jω, jχ, jψ, :, :]
			ct.gπ .= mt.g[jω, jχ, jψ, :, :]
		end

		flag = pfi!(ct, verbose = false)
		verbose && flag && print(": ✓")
		verbose && !flag && print(": no convergence.")
		
		mt.L[jω, jχ, jψ, :, :] .= ct.L
		mt.C[jω, jχ, jψ, :, :] .= ct.C
		mt.g[jω, jχ, jψ, :, :] .= ct.gπ
		
		perc = 100 * iter / tot
		verbose && print(" $(@sprintf("%.3g",perc))% completed.\n")
		
        save_progress && (jψ == length(mt.ψgrid)) && save("temp.jld2", "mt", mt)
	end
end


Bayes_plan(ν, z, μ) = z*ν / (z*ν + (1-z)*μ)

function eval_k_to_mu(mt::MultiType, k, itp_L; get_mu::Bool=false, verbose::Bool=false)

	ωgrid, χgrid, L_mat = mt.ωgrid, mt.χgrid, mt.L_mat
	pgrid, agrid = mt.ct.gr[:p], mt.ct.gr[:a]

	μ, p0 = [zeros(length(ωgrid), length(χgrid), length(agrid)) for jj in 1:2]

	for (ja, av) in enumerate(agrid), (jχ, χv) in enumerate(χgrid), (jω, ωv) in enumerate(ωgrid)
		if L_mat[jω, jχ, end, ja] > k
			pv = 0.0
			μ[jω, jχ, ja] = 0.0
		else
			res = Optim.optimize(
				p -> (itp_L(ωv, χv, p, av) - k)^2, 0, 1, GoldenSection())

			disp = sqrt(res.minimum)
			if verbose && disp > 1e-4
				print("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)\n")
			end
			pv = res.minimizer # p0 to get value k in the current type

			νv = mt.ν[jω, jχ, ja] # Get density of the current type
			μv = (1-pv)/pv * mt.z/(1-mt.z) * νv # probability of announcement to start with p0
			if (Bayes_plan(νv, mt.z, μv) - pv)^2 > 1e-4
				print("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)\n")
			end
			# res = Optim.optimize(
			# 	μ -> (Bayes_plan(νv, mt.z, μ) - pv)^2, 0, 1, GoldenSection())
			# disp = res.minimum
			# if verbose && disp > 1e-4
			# 	print("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)\n")
			# end
			# μv = res.minimizer # probability of announcement to start with p0

			μ[jω, jχ, ja] = μv
		end
		p0[jω, jχ, ja] = pv # Save p0
	end
	if get_mu
		knts = (ωgrid[end:-1:1], χgrid, pgrid, agrid)
		itp_C = interpolate(knts, mt.C_mat[end:-1:1,:,:,:], Gridded(Linear()))
		C_eqm = [ itp_C(ωv, χv, p0[jω, jχ, ja], av) for (jω, ωv) in enumerate(ωgrid), (jχ, χv) in enumerate(χgrid), (ja, av) in enumerate(agrid) ]
		C_mat = [ sum( [itp_C(ωv, χv, p0[jω, jχ, ja], av) * sum(μ[jω, jχ, ja]) for (ja, av) in enumerate(agrid)] ) for (jω, ωv) in enumerate(ωgrid), (jχ, χv) in enumerate(χgrid) ]
		return μ, C_eqm, C_mat
	else
		return sum(μ)
	end
end

function find_plan_μ(mt::MultiType; annualize::Bool=false, decay::Bool=false)
	pgrid, agrid = mt.ct.gr[:p], mt.ct.gr[:a]
	ωgrid, χgrid = mt.ωgrid, mt.χgrid

	if decay
		ωgrid = perc_rate.(ωgrid)
	end
	if annualize
		agrid = annualized.(agrid)
		χgrid = annualized.(χgrid)
	end

	mean_ω, mean_χ, mean_a = zeros(3)
	m2_ω, m2_χ, m2_a = zeros(3)

	sum_prob = sum(mt.μ)
	for (ja, av) in enumerate(agrid), (jχ, χv) in enumerate(χgrid), (jω, ωv) in enumerate(ωgrid)

		mean_ω += mt.μ[jω, jχ, ja] * ωv
		mean_a += mt.μ[jω, jχ, ja] * av
		mean_χ += mt.μ[jω, jχ, ja] * χv

		m2_ω += mt.μ[jω, jχ, ja] * ωv^2
		m2_a += mt.μ[jω, jχ, ja] * av^2
		m2_χ += mt.μ[jω, jχ, ja] * χv^2
	end

	mean_ω *= 1/sum_prob
	mean_a *= 1/sum_prob
	mean_χ *= 1/sum_prob
	m2_ω *= 1/sum_prob
	m2_a *= 1/sum_prob
	m2_χ *= 1/sum_prob

	sd_ω = sqrt(m2_ω - mean_ω^2)
	sd_a = sqrt(m2_a - mean_a^2)
	sd_χ = sqrt(m2_χ - mean_χ^2)

	return mean_ω, mean_a, mean_χ, sd_ω, sd_a, sd_χ
end

function find_equil!(mt::MultiType, z0=mt.ct.gr[:p][3])
	mt.z = z0
	pgrid, agrid = mt.ct.gr[:p], mt.ct.gr[:a]
	ωgrid, χgrid = mt.ωgrid, mt.χgrid
	L_mat = mt.L_mat

	k_max = mean(L_mat[:,:,1,:]) # mean L when p = 0 (should be constant across plans)
	k_min = minimum(L_mat[:,:,end,:]) # lower loss 
	V = var(L_mat[:,:,1,:])
	if V > 1e-4
		print("WARNING: variance of Lᴺ = $(@sprintf("%0.3g",V)), should be 0\n")
	end

	knots = (ωgrid[end:-1:1], χgrid, pgrid, agrid)
	itp_L = interpolate(knots, L_mat[end:-1:1,:,:,:], Gridded(Linear()))

	res = Optim.optimize(
		k -> (eval_k_to_mu(mt, k, itp_L)-1)^2,
		k_min, k_max, GoldenSection())

	if res.minimum > 1e-4
		print("WARNING: Couldn't find μ at z = $z0\n")
	end

	k_star = res.minimizer

	mt.μ, C_eqm, C_mat = eval_k_to_mu(mt, k_star, itp_L; get_mu = true, verbose=true)

	return k_star
end

function mimic_z(mt::MultiType, N=100; decay::Bool=false, annualize::Bool=false)

	zgrid = cdf.(Beta(4,1), range(0,1,length=N))
	# move_grids!(zgrid, xmax=0.9, xmin=mt.ct.pgrid[3])
	move_grids!(zgrid, xmax=0.95, xmin=1e-9)

	data = zeros(N,6)
	datanames = ["ω", "a", "χ", "s_ω", "s_a", "s_χ"]

	for (jz, zv) in enumerate(zgrid)
		find_equil!(mt, zv)
		data[jz,:] .= find_plan_μ(mt, decay=decay, annualize=annualize)
	end

	return data, datanames, zgrid
end

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


function extend_ψ(ct::CrazyType; ψmin = 0.0, ψmax = 0.35, Nψ = 100)

	Lψ = Float64[]
	ψvec = Float64[]

	print("Preparing run.\n")

	ct.pars[:ψ] = 0.
	pfi!(ct, verbose = false)

	ψgrid = range(ψmin, ψmax, length = Nψ)

	for (jψ, ψv) in enumerate(ψgrid)

		frac = @sprintf("%.3g", 100 * jψ / Nψ)

		print("Starting with ψ = $(@sprintf("%.3g", ψv))")

		ct.pars[:ψ] = ψv

		flag = pfi!(ct, verbose = false)

		if flag
			push!(Lψ, minimum(ct.L[2,:]))
			push!(ψvec, ψv)

			print(": convergence reached.")
		else
			print(": failed.")
		end
		print(" Done with $frac%\n")
	end

	return Lψ, ψvec
end