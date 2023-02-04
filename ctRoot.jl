using Distributions, Interpolations, Optim, HCubature, QuantEcon, Printf, PlotlyJS, Dates, JLD2, LinearAlgebra

include("type_def.jl")
include("reporting_routines.jl")
include("simul.jl")
include("plotsct.jl")
# include("planner.jl")

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

next_a(ct::Plan, av, apv, π) = ifelse(ct.use_a, ϕ(ct, av), ϕ(ct, av) + ct.ψ * (π-av))
next_a(ct::DovisKirpalani, av, apv, π) = apv

function cond_L(ct::Plan{T}, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime, ge, πe′) where {T<:PhillipsCurve}
    # ge = itp_gπ(pv, av)
    pprime = Bayes(ct, obs_π, ge, pv, av)

    πe = pv * av + (1 - pv) * ge

    a_min, a_max = extrema(ct.agrid)

    π_today = max(a_min, min(a_max, obs_π))
    aprime = next_a(ct, av, aprime, π_today)
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
    L = (ct.ystar - y)^2 + ct.γ * obs_π^2 + ct.β * L′
    C′ = itp_C(pprime, aprime)

    return L, pprime, y, C′
end

get_sumprob(ct::Plan) = cdf_ϵ(ct, 3.09*ct.σ) - cdf_ϵ(ct, -3.09*ct.σ)

function exp_L_y(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, aprime, ge, πe′)

	sum_prob = get_sumprob(ct)

	f_y(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime, ge, πe′)[3] * pdf_ϵ(ct, ϵv)
	Ey, err = hquadrature(f_y, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	f_p(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime, ge, πe′)[2] * pdf_ϵ(ct, ϵv)
	Ep, err = hquadrature(f_p, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	f_C(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime, ge, πe′)[4] * pdf_ϵ(ct, ϵv)
	Ec, err = hquadrature(f_C, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)

	Ey = Ey / sum_prob
	Ep = Ep / sum_prob
	Ec = Ec / sum_prob

	return Ey, Ep, Ec
end

function exp_L(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, aprime, ge, πe′)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime, ge, πe′)[1] * pdf_ϵ(ct, ϵv)
	val, err = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	sum_prob = get_sumprob(ct)

	return val/sum_prob
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′)
    π_guess = xguess[1]
	aprime = xguess[2]

    gπ, L = π_guess, itp_L(pv, av)

    obj_f(x) = exp_L(ct, itp_gπ, itp_L, itp_C, first(x), pv, av, aprime, ge, πe′)
	
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
	val::Float64, err::Float64 = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	sum_prob = get_sumprob(ct)

	return val/sum_prob
end

function optim_step(ct::Plan, itp_gπ, itp_L, itp_C, optim = true)
	gπ, ga = zeros(size(ct.gπ)), zeros(size(ct.ga))
	L  	   = zeros(size(ct.L))
	Ey, Eπ = zeros(size(ct.Ey)), zeros(size(ct.Eπ))
	Ep, C  = zeros(size(ct.Ep)), zeros(size(ct.C))
	πN 	   = Nash(ct)
	
	h = 0.05 * πN
	Gmin, Gmax = -h, πN + h
	
	apgrid = gridmake(1:ct.Np, 1:ct.Na)
	Threads.@threads for js in axes(apgrid,1)
		# for js in axes(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.pgrid[jp], ct.agrid[ja]
		
		aprime = ct.ga[jp, ja]
		π_guess = ct.gπ[jp, ja]

		xguess = [π_guess, aprime]

		ge = π_guess
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

		Ey[jp, ja], Ep[jp, ja], EC′ = exp_L_y(ct, itp_gπ, itp_L, itp_C, Gc, pv, av, aprime, ge, πe′)
		Eπ[jp, ja] = pv * av + (1.0 - pv) * Gc

		if av >= πN || isapprox(av, πN)
			C[jp, ja] = (1-ct.β)*1 + ct.β * EC′
		else
			C[jp, ja] = (1-ct.β)*(πN - Eπ[jp,ja])/(πN-av) + ct.β * EC′
		end
	end
	return gπ, L, Ey, Eπ, Ep, C, ga
end

function pf_iter(ct::Plan; optim=true)
	knots = (ct.pgrid, ct.agrid);
	itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()));
	itp_L  = interpolate(knots, ct.L, Gridded(Linear()));
	itp_C  = interpolate(knots, ct.C, Gridded(Linear()));

	new_gπ, new_L, new_y, new_π, new_p, new_C, new_a = optim_step(ct, itp_gπ, itp_L, itp_C, optim)

	return new_gπ, new_L, [new_y, new_π, new_p, new_C, new_a]
end

function update_others!(ct::Plan, new_others, upd_η2)
	new_y, new_π, new_p, new_C, new_a = new_others[:]
	ct.Ey = new_y
	ct.Eπ = new_π
	ct.Ep = new_p
	ct.C  = new_C
	ct.ga = ct.ga + upd_η2 * (new_a - ct.ga)
	nothing
end

function pfi!(ct::Plan; miniter::Int = 2, tol::Float64=1e-4, maxiter::Int64=2000, verbose::Bool=true, accelerate = true)
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
			_, new_L, _ = pf_iter(ct, optim = false)
			ct.L  = upd_η * new_L  + (1.0-upd_η) * ct.L
			end
		end
		
		new_gπ, new_L, new_others = pf_iter(ct, optim = true)
		update_others!(ct, new_others, upd_η)

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

	return (dist <= tol)
end

function solve_all!(mt::MultiType; verbose = true)
	verbose && print("Going over all plans at $(Dates.format(now(), "HH:MM"))\n")
	iter = 0
	tot  = length(mt.ωgrid) * length(mt.χgrid)
	for (jω, ωv) in enumerate(mt.ωgrid), (jχ, χv) in enumerate(mt.χgrid)
		iter += 1
		
		show_ω = @sprintf("%.3g", ωv)
		show_χ = @sprintf("%.3g", annualized(χv))
		
		verbose && print("Plan with (ω, χ) = ($show_ω, $show_χ%)")
		ct = mt.ct

		if jω > 1 && jχ == 1
			ct.L .= mt.L_mat[jω - 1, jχ, :, :]
		end

		ct.ω = ωv
		ct.χ = χv

		flag = pfi!(ct, verbose = false)
		verbose && flag && print(": ✓")
		verbose && !flag && print(": no convergence.")
		
		mt.L_mat[jω, jχ, :, :] .= ct.L
		mt.C_mat[jω, jχ, :, :] .= ct.C
		mt.g_mat[jω, jχ, :, :] .= ct.gπ
		
		perc = 100 * iter / tot
		verbose && print(" $perc% completed.\n")
	end
end

Bayes_plan(ν, z, μ) = z*ν / (z*ν + (1-z)*μ)

function eval_k_to_mu(mt::MultiType, k, itp_L; get_mu::Bool=false, verbose::Bool=false)

	ωgrid, χgrid, L_mat = mt.ωgrid, mt.χgrid, mt.L_mat
	pgrid, agrid = mt.ct.pgrid, mt.ct.agrid

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
				print_save("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)")
			end
			pv = res.minimizer # p0 to get value k in the current type

			νv = mt.ν[jω, jχ, ja] # Get density of the current type
			res = Optim.optimize(
				μ -> (Bayes_plan(νv, mt.z, μ) - pv)^2, 0, 1, GoldenSection())
			disp = res.minimum
			if verbose && disp > 1e-4
				print_save("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)")
			end
			μv = res.minimizer # probability of announcement to start with p0

			μ[jω, jχ, ja] = μv
		end
		p0[jω, jχ, ja] = pv # Save p0
	end
	if get_mu
		knots = (ωgrid[end:-1:1], χgrid, pgrid, agrid)
		itp_C = interpolate(knots, mt.C_mat[end:-1:1,:,:,:], Gridded(Linear()))
		C_eqm = [ itp_C(ωv, χv, p0[jω, jχ, ja], av) for (jω, ωv) in enumerate(ωgrid), (jχ, χv) in enumerate(χgrid), (ja, av) in enumerate(agrid) ]
		C_mat = [ sum( [itp_C(ωv, χv, p0[jω, jχ, ja], av) * sum(μ[jω, jχ, ja]) for (ja, av) in enumerate(agrid)] ) for (jω, ωv) in enumerate(ωgrid), (jχ, χv) in enumerate(χgrid) ]
		return μ, C_eqm, C_mat
	else
		return sum(μ)
	end
end

function find_plan_μ(mt::MultiType; annualize::Bool=false, decay::Bool=false)
	pgrid, agrid = mt.ct.pgrid, mt.ct.agrid
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

function find_equil!(mt::MultiType, z0=mt.ct.pgrid[3])
	mt.z = z0
	pgrid, agrid = mt.ct.pgrid, mt.ct.agrid
	ωgrid, χgrid = mt.ωgrid, mt.χgrid
	L_mat = mt.L_mat

	jp0 = floor(Int, length(mt.ct.pgrid)*0.9)

	k_max = mean(L_mat[:,:,1,:]) # mean L when p = 0 (should be constant across plans)
	k_min = minimum(L_mat[:,:,jp0,:]) # lower loss 
	V = var(L_mat[:,:,1,:])
	if V > 1e-4
		print_save("WARNING: variance of Lᴺ = $(@sprintf("%0.3g",V)), should be 0")
	end

	knots = (ωgrid[end:-1:1], χgrid, pgrid, agrid)
	itp_L = interpolate(knots, L_mat[end:-1:1,:,:,:], Gridded(Linear()))

	res = Optim.optimize(
		k -> (eval_k_to_mu(mt, k, itp_L)-1)^2,
		k_min, k_max, GoldenSection())

	if res.minimum > 1e-4
		print_save("WARNING: Couldn't find μ at z = $z0")
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
