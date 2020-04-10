using Distributed

using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD

include("type_def.jl")
include("reporting_routines.jl")
include("simul.jl")
include("plotting_routines.jl")
include("planner.jl")

function output_bayes(ct::CrazyType, pv, av)
	knots = (ct.pgrid, ct.agrid);
	itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()));

	# exp_π = pv*av + (1-pv)*itp_gπ(pv, av)
	exp_π = itp_gπ(pv, av)

	println("gπ = [$(annualized(itp_gπ(pv, av))-1.96*ct.σ), $(annualized(itp_gπ(pv, av))+1.96*ct.σ)], av = $(annualized(av))")

	println("$(pdf_ϵ(ct, exp_π - av ))")
	println("$(pdf_ϵ(ct, 0.0 ))")

	aprime = ϕ(ct, av)
	π_myopic = pv * aprime + (1.0-pv) * itp_gπ(pv, aprime)

	Nv = 50
	yv = zeros(Nv)
	ym = zeros(Nv)
	πvec = range(av - 1.96*ct.σ, av + 1.96*ct.σ, length=Nv)
	for (jj, πv) in enumerate(πvec)

		pprime = Bayes(ct, πv, exp_π, pv, av)
		exp_π′ = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)
		yv[jj] = PC(ct, πv, exp_π, exp_π′)
		ym[jj] = PC(ct, πv, exp_π, π_myopic)

		# yv[jj] = pdf_ϵ(ct, πv - av)
		# yv[jj] = pprime

	end

	plot([
		scatter(;x=annualized.(πvec), y=yv)
		# scatter(;x=annualized.(πvec), y=ym)
		])
end

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

function cond_Ldev(ct::CrazyType, itp_gπ, itp_L, obs_π, pv, av)
	aprime = ϕ(ct, av)

	πe = pv*av + (1-pv)*exp_π
	exp_π′ = itp_gπ(0.0, aprime)

	y = PC(ct, obs_π, πe, exp_π′) # Automatically uses method for forward or backward
	L′ = itp_L(0.0, aprime)

	L = (ct.ystar-y)^2 + ct.γ * obs_π^2 + ct.β * L′

	return L
end

function cond_L_inner(ct::Plan, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime)
	exp_π  = itp_gπ(pv, av)
	pprime = Bayes(ct, obs_π, exp_π, pv, av)

	πe = pv*av + (1-pv)*exp_π

	if aprime <= minimum(ct.agrid) || aprime >= maximum(ct.agrid)
		itp_L  = extrapolate(itp_L,  Interpolations.Flat())
		itp_gπ = extrapolate(itp_gπ, Interpolations.Flat())
		itp_C  = extrapolate(itp_C, Interpolations.Flat())
	end

	L′::Float64 = itp_L(pprime, aprime)
	exp_π′::Float64 = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)

	y = PC(ct, obs_π, πe, exp_π′) # Automatically uses method for forward or backward
	L = (ct.ystar-y)^2 + ct.γ * obs_π^2 + ct.β * L′
	C′::Float64 = itp_C(pprime, aprime)

	return L, pprime, y, C′
end

function cond_L(ct::Plan, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime)
	L, pprime, y, C′ = cond_L_inner(ct, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime)
	return L
end	
function cond_L_others(ct::Plan, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime)
	L, pprime, y, C′ = cond_L_inner(ct, itp_gπ, itp_L, itp_C, obs_π, pv, av, aprime)
	return y, pprime, C′
end

get_sumprob(ct::Plan) = cdf_ϵ(ct, 3.09*ct.σ) - cdf_ϵ(ct, -3.09*ct.σ)

function exp_L_y(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, aprime)

	sum_prob = get_sumprob(ct)

	f_y(ϵv) = cond_L_others(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime)[1] * pdf_ϵ(ct, ϵv)
	Ey, err = hquadrature(f_y, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	f_p(ϵv) = cond_L_others(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime)[2] * pdf_ϵ(ct, ϵv)
	Ep, err = hquadrature(f_p, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	f_C(ϵv) = cond_L_others(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime)[3] * pdf_ϵ(ct, ϵv)
	Ec, err = hquadrature(f_C, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)

	Ey = Ey / sum_prob
	Ep = Ep / sum_prob
	Ec = Ec / sum_prob

	return Ey, Ep, Ec
end

function exp_L(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, aprime)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, aprime) * pdf_ϵ(ct, ϵv)
	val::Float64, err::Float64 = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	sum_prob = get_sumprob(ct)

	return val/sum_prob
end

function opt_L(ct::DovisKirpalani, itp_gπ, itp_L, itp_C, xguess, pv, av)

	minπ = max(0, xguess[1] - 3.09*ct.σ)
	maxπ = min(1.1*maximum(ct.agrid), xguess[1] + 3.09*ct.σ)
	if maxπ < minπ + 1.1*maximum(ct.agrid) / 10
		maxπ = minπ + 1.1*maximum(ct.agrid) / 10
	end
	mina = minimum(ct.agrid)
	maxa = maximum(ct.agrid)

	obj_f(x) = exp_L(ct, itp_gπ, itp_L, itp_C, x[1], pv, av, x[2])
	res = Optim.optimize(obj_f, [mina, mina], [maxa, maxa], xguess, Fminbox(NelderMead()))

	gπ::Float64, ga::Float64 = res.minimizer
	L::Float64 = res.minimum

	if Optim.converged(res) == false
		resb = Optim.optimize(obj_f, [mina, mina], [maxa, maxa], xguess, Fminbox(LBFGS()))
		if resb.minimum < res.minimum
			gπ, ga = resb.minimizer
			L = resb.minimum
		end
	end

	return gπ, L, ga
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, itp_C, xguess, pv, av)
	π_guess = xguess[1]
	minπ = max(0, π_guess - 3.09*ct.σ)
	maxπ = min(1.1*maximum(ct.agrid), π_guess + 3.09*ct.σ)
	if maxπ < minπ + 1.1*maximum(ct.agrid) / 10
		maxπ = minπ + 1.1*maximum(ct.agrid) / 10
	end
	
	# aprime = ϕ(ct, av)
	aprime = xguess[2]

	obj_f(x) = exp_L(ct, itp_gπ, itp_L, itp_C, x, pv, av, aprime)
	res = Optim.optimize(
		gπ -> obj_f(first(gπ)),
		[π_guess], LBFGS()#, autodiff=:forward#, Optim.Options(f_tol=1e-12)
		)

	gπ::Float64, L::Float64 = first(res.minimizer), res.minimum

	if Optim.converged(res) == false
		resb = Optim.optimize(
				gπ -> exp_L(ct, itp_gπ, itp_L, itp_C, gπ, pv, av, aprime),
				minπ, maxπ, Brent(), rel_tol=1e-12, abs_tol=1e-12#, iterations=100000
				)
		if resb.minimum < res.minimum
			gπ, L = resb.minimizer, resb.minimum
		end
	end

	return gπ, L, aprime
end

function optim_step(ct::Plan, itp_gπ, itp_L, itp_C, gπ_guess; optimize::Bool=true)
	gπ, ga = zeros(size(ct.gπ)), zeros(size(ct.ga))
	L  	   = zeros(size(ct.L))
	Ey, Eπ = zeros(size(ct.Ey)), zeros(size(ct.Eπ))
	Ep, C  = zeros(size(ct.Ep)), zeros(size(ct.C))
	πN 	   = Nash(ct)
	maxa = maximum(ct.agrid)
	mina = minimum(ct.agrid)
	length_a = maxa-mina
	apgrid = gridmake(1:ct.Np, 1:ct.Na)
	Threads.@threads for js in 1:size(apgrid,1)
    # for js in 1:size(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.pgrid[jp], ct.agrid[ja]

		a_guess = max(min(ct.ga[jp, ja], maxa),mina)
		π_guess = max(min(gπ_guess[jp, ja], maxa),mina)
		xguess = [π_guess, a_guess]
		if optimize
			# π_guess = itp_gπ(pv, av)
			gπ[jp, ja], L[jp, ja], aprime = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av)
		else
			aprime = a_guess
			gπ[jp, ja] = π_guess
			L[jp, ja] = exp_L(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av, aprime)
		end
		ga[jp, ja] = aprime
		Ey[jp, ja], Ep[jp, ja], EC′ = exp_L_y(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av, aprime)
		Eπ[jp, ja] = pv * av + (1.0 - pv) * gπ[jp, ja]

		if av >= πN || isapprox(av, πN)
			C[jp, ja] = (1-ct.β)*1 + ct.β * EC′
		else
			C[jp, ja] = (1-ct.β)*(πN - Eπ[jp,ja])/(πN-av) + ct.β * EC′
		end
	end

	return gπ, L, Ey, Eπ, Ep, C, ga
end


function pf_iter(ct::Plan, Egπ, gπ_guess; optimize::Bool=true)
	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, Egπ, Gridded(Linear()))
	itp_L  = interpolate(knots, ct.L, Gridded(Linear()))
	itp_C  = interpolate(knots, ct.C, Gridded(Linear()))

	new_gπ, new_L, new_y, new_π, new_p, new_C, new_a = optim_step(ct, itp_gπ, itp_L, itp_C, gπ_guess; optimize=optimize)

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

function pfi!(ct::Plan, Egπ; tol::Float64=1e-12, maxiter::Int64=300, verbose::Bool=true, reset_guess::Bool=false)
	dist = 10.
	iter = 0
	upd_η2 = 0.75

	rep = "\nStarting PFI (tol = $(@sprintf("%0.3g",tol)))"
	verbose ? print_save(rep,true) : print(rep)

    if reset_guess
	    # ct.gπ = zeros(size(ct.gπ))
		ct.L = ones(ct.Np, ct.Na)
	end

	old_gπ = copy(Egπ)
	new_gπ = zeros(size(old_gπ))

	while dist > tol && iter < maxiter
		iter += 1

		for jj in 1:10
			_, new_L, _ = pf_iter(ct, Egπ, old_gπ; optimize=false)
			ct.L  = upd_η2 * new_L  + (1.0-upd_η2) * ct.L
		end
		# println("Iter $iter step")
		old_L = copy(ct.L)

		new_gπ, new_L, new_others = pf_iter(ct, Egπ, old_gπ)
		update_others!(ct, new_others, upd_η2)

		norm_L = max(sqrt.(sum(old_L .^2)) / length(old_L), 100tol)
		dist = sqrt.(sum( (new_L  - old_L ).^2 ))/length(old_L) / norm_L

		ct.L  = upd_η2 * new_L  + (1.0-upd_η2) * ct.L
		old_gπ = upd_η2 * new_gπ + (1.0-upd_η2) * old_gπ

		# if verbose && iter % 10 == 0
		# 	print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
		# end
	end

	if verbose && dist <= tol
		print("\nConverged in $iter iterations.")
	elseif verbose
		print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
	end

	return (dist <= tol), new_gπ
end

decay_η(ct::Plan, η) = max(0.95*η, 5e-6)

function report_start(ct::CrazyType)
	print_save("\nRun with ω = $(@sprintf("%.3g",ct.ω)), χ = $(@sprintf("%.3g",annualized(ct.χ)))% at $(Dates.format(now(), "HH:MM"))")
	nothing
end
function report_start(dk::DovisKirpalani)
	nothing
end

function reset_L!(pp::Plan)
	nothing
end
function reset_L!(dk::DovisKirpalani)
	dk.L *= 0
	nothing
end

function solve!(dk::DovisKirpalani; tol::Float64=5e-4, maxiter::Int64=2500)
	dist = 10.
	iter = 0

	tol_epfi = 1e-3

	ct = CrazyType(dk)
	Epfi!(ct, maxiter = 500)

	while dist > tol && iter < maxiter
		iter += 1

		old_ga = copy(dk.ga)

		ct = CrazyType(dk)
		dist_π = Epfi!(ct, tol=tol_epfi)

		# dk = DovisKirpalani(ct)
		# dk.ga = ct.ga
		dk.gπ = copy(ct.gπ)
		dk.L  = copy(ct.L)

		Epfi!(dk; maxiter = 1, tol_pfi = tol_epfi)

		norm_ga = max(sqrt.(sum(annualized.(old_ga) .^2)) / length(annualized.(old_ga)), 10tol)
		dist_a = sqrt.(sum( (annualized.(dk.ga)  - annualized.(old_ga) ).^2 ))/length(annualized.(old_ga)) / norm_ga

		# rep_status = "\nAfter $iter iterations, d(π) = $(@sprintf("%0.3g",dist)) at |a| = $(@sprintf("%0.3g",norm_ga)))"
		rep_status = "\nAfter $iter iterations, d(a) = $(@sprintf("%0.3g",dist)) at |a| = $(@sprintf("%0.3g",norm_ga))"

		dist_π <= tol_epfi ? rep_status *= " ✓" : nothing

		print_save(rep_status)

		dist = max(dist_a, dist_π)
		tol_epfi *= 0.95
		tol_epfi = max(5e-4, tol_epfi)
	end
end

function Epfi!(ct::Plan; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, tempplots::Bool=false, upd_η::Float64=0.01, switch_η = 10, tol_pfi = 2e-3 / 0.99)
	dist = 10.
	iter = 0
	
	report_start(ct)
	dists = []

	reset_guess = false
	while dist > tol && iter < maxiter
		iter += 1
		tol_pfi = max(tol_pfi*0.98, 2e-6)

		old_gπ, old_L, old_ga = copy(ct.gπ), copy(ct.L), copy(ct.ga);

		# reset_L!(ct)

		flag, new_gπ = pfi!(ct, old_gπ; verbose=verbose, reset_guess=reset_guess, tol=tol_pfi);
		reset_guess = !flag

		norm_gπ = max(sqrt.(sum(annualized.(ct.gπ) .^2)) / length(annualized.(ct.gπ)), 20tol)
		dist_π = sqrt.(sum( (annualized.(new_gπ)  - annualized.(ct.gπ) ).^2 ))/length(annualized.(ct.gπ)) / norm_gπ
		norm_ga = max(sqrt.(sum(annualized.(old_ga) .^2)) / length(annualized.(old_ga)), 20tol)
		dist_a = sqrt.(sum( (annualized.(ct.ga)  - annualized.(old_ga) ).^2 ))/length(annualized.(old_ga)) / norm_ga
		dist = max(dist_π, dist_a/10)
		push!(dists, dist)
		rep_status = "\nAfter $iter iterations, d(π) = $(@sprintf("%0.3g",dist)) at |π,a| = ($(@sprintf("%0.3g",norm_gπ)), $(@sprintf("%0.3g",norm_ga)))"
		if flag
			rep_status *= "✓ "
		end
		if verbose #&& iter % 10 == 0
			print_save(rep_status*"\n", true)
		else
			print(rep_status)
		end

		ct.gπ = upd_η * new_gπ + (1-upd_η) * ct.gπ;
		ct.ga = upd_η * ct.ga + (1-upd_η) * old_ga;

		if tempplots && (iter % 5 == 0 || dist <= tol)
			p1, pL, pE, pC, pp, _ = makeplots_ct_pa(ct);
			relayout!(p1, title="iter = $iter")
			savejson(p1, pwd()*"/../Graphs/tests/temp.json")
			relayout!(pL, title="iter = $iter")
			savejson(pL, pwd()*"/../Graphs/tests/tempL.json")
			# relayout!(pE, title="iter = $iter")
			# savejson(pE, pwd()*"/../Graphs/tests/tempLpE.json")
			p2 = makeplot_conv(dists; switch_η=switch_η);
			savejson(p2, pwd()*"/../Graphs/tests/tempconv.json")
		end

		if iter == floor(Int, switch_η*0.4)
			upd_η = min(upd_η, 0.01)
		elseif iter % switch_η == 0
			upd_η = decay_η(ct, upd_η) # Automatically uses the updating method for fwd or bwd
		end
		if verbose
			print_save("new upd_η = $(@sprintf("%0.3g", upd_η))", true)
		end

	end

	tolC, maxiterC = 1e-5, 2000
	dist2, iter = 1+tolC, 0
	while dist2 > tolC && iter < maxiterC
		iter += 1

		old_C = copy(ct.C)
		Egπ = copy(ct.gπ)
		_, _, new_others = pf_iter(ct, Egπ, Egπ, optimize = false)

		new_C = new_others[4]
		dist2 = sqrt.(sum( (new_C  - old_C ).^2 )) / sqrt.(sum(old_C .^2))

		update_others!(ct, new_others, 0.5)
	end

	if verbose && dist <= tol
		print_save("\nConverged in $iter iterations.",true)
	elseif verbose
		print_save("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))",true)
	end
	if dist <= tol
		p1, pL, pπ, pC, pp, _ = makeplots_ct_pa(ct);
		savejson(pC, pwd()*"/../Graphs/tests/tempC.json")
		savejson(pπ, pwd()*"/../Graphs/tests/tempg.json")
	end
	
	return dist
end

# function choose_ω!(L_mat, ct::CrazyType{Forward}, Nω=size(L_mat,1); remote::Bool=true, upd_η=0.1)
# 	choose_ω!(L_mat, ct, Forward, Nω; remote=remote, upd_η=upd_η)
# end

# function choose_ω!(L_mat, ct::CrazyType{Simultaneous}, Nω=size(L_mat,1); remote::Bool=true, upd_η=0.1)
# 	choose_ω!(L_mat, ct, Simultaneous, Nω; remote=remote, upd_η=upd_η)
# end

function choose_ω!(L_mat, ct::CrazyType, Nω=size(L_mat,1); upd_η=0.1)
	T = which_PC(ct)
	ct_best = CrazyType(T; γ=ct.γ, κ=ct.κ, σ=ct.σ, β=ct.β, ystar=ct.ystar)

	if T == Simultaneous
		ωmax = 3.0
	elseif T == Forward
		ωmax = 1.5
	end
	ωgrid = cdf.(Beta(1,1), range(1,0,length=Nω))
	move_grids!(ωgrid, xmax = ωmax, xmin = 0.01)

	Na = length(ct.agrid)
	Nχ = size(L_mat, 2)
	χgrid = range(0.0, 0.43*Nash(ct), length = Nχ)

	update_ga!(ct, ω = ωgrid[1], χ = χgrid[1])
	dist = Epfi!(ct)
	print_save("\nDone with initial setup $(ifelse(dist<1e-3, "✓", ""))")

	print_save("\nLooping over behavioral types with ω ∈ [$(minimum(ωgrid)), $(maximum(ωgrid))]")
	print_save("\n")

	L_min = 100.
	ω_min = 1.0
	χ_min = 1.0
	a_min = 1.0
	t0 = time()
	Lplot = []
	aplot = []
	C_mat = NaN * L_mat
	L_mat_ctour = zeros(Nω, Nχ) * NaN
	C_mat_ctour = zeros(Nω, Nχ) * NaN
	Lmin = 1e8
	ja_min = 1
	for (jχ, χv) in enumerate(χgrid)
		L_vec = []
		a_vec = []
		ω_vec = []

		""" tol = 11e-4 """
		function wrap_Epfi!(ct::CrazyType, ωv, L_vec, a_vec, ω_vec, Lplot, L_mat_save, C_mat, aplot, jω, jχ)
			update_ga!(ct, ω = ωv)

			t1 = time()
			tol = 5e-4
			# if length(L_vec) > 0
			# 	upd_η = 0.005
			# end
			dist = Epfi!(ct, verbose = true, tol=tol, tempplots=false, upd_η=upd_η)
			write(pwd()*"/../temp.txt", "")
			
			flag = (dist <= tol)
			Lmin, ja = findmin(ct.L[3,:])
			Cmin = ct.C[3,ja]
			# Cmin = ct.C[3,end]

			for jp in 1:length(ct.pgrid), ja in 1:length(ct.agrid)
				C_mat[jω, jχ, jp, ja] = ct.C[jp, ja]
			end
			
			s = ": done in $(time_print(time()-t1))"
			flag ? s = s*" ✓" : nothing
			print_save(s)

			L_mat_save[:,:] = ct.L

			push!(L_vec, Lmin)
			push!(a_vec, ct.agrid[ja])
			push!(ω_vec, ωv)

			perm_order = sortperm(ω_vec)

			new_L = scatter(;x=ω_vec[perm_order], y=L_vec[perm_order], name = "χ = $(@sprintf("%.3g",annualized(χv)))%", line_shape="spline")
			new_a = scatter(;x=ω_vec[perm_order], y=annualized.(a_vec[perm_order]), name = "χ = $(@sprintf("%.3g",annualized(χv)))%")

			all_Ls = new_L
			all_as = new_a
			if length(Lplot) == 0
			else
				all_Ls = vcat([Lplot[jj] for jj in 1:length(Lplot)], new_L)
				all_as = vcat([aplot[jj] for jj in 1:length(aplot)], new_a)
			end
			p3 = plot(all_Ls)
			relayout!(p3, title="lim_𝑝 min_𝑎 𝓛(𝑝,𝑎,ω,χ)", xaxis=attr(;zeroline=false, title="ω"))
			savejson(p3, pwd()*"/../Graphs/tests/Loss_omega.json")
	
			p4 = plot(all_as)
			relayout!(p4, title="lim_𝑝 arg min_𝑎 𝓛(𝑝,𝑎,ω,χ)", xaxis=attr(;zeroline=false, title="ω"), yaxis_title="%", mode="lines+markers")
			savejson(p4, pwd()*"/../Graphs/tests/a0.json")

			return Lmin, Cmin, ja, flag
		end

		ωmin = 1e8
		amin = 1e8
		for (jω, ωv) in enumerate(ωgrid)
			ωv = ωgrid[jω]
			old_L, old_gπ = copy(ct.L), copy(ct.gπ)
			if jω == 1 && jχ > 1
				old_ct = load("../ct_1_temp.jld", "ct")
				old_L, old_gπ = copy(old_ct.L), copy(old_ct.gπ)
			end

			ct = CrazyType(T; χ=χv, γ=ct.γ, κ=ct.κ, σ=ct.σ, β=ct.β, ystar=ct.ystar)
			
			ct.L, ct.gπ = old_L, old_gπ
			
			L_mat_save = zeros(ct.Np, ct.Na)
			L, C, ja, flag = wrap_Epfi!(ct, ωv, L_vec, a_vec, ω_vec, Lplot, L_mat_save, C_mat, aplot, jω, jχ)

			L_mat[jω, jχ, :, :] = L_mat_save
			L_mat_ctour[jω, jχ] = L

			C_mat_ctour[jω, jχ] = C 

			for slides in [true, false]
				pLct = plot_L_contour(ωgrid, χgrid, L_mat_ctour, name_y="𝓛", slides=slides)
				savejson(pLct, pwd()*"/../Graphs/tests/contour$(ifelse(slides, "_slides", "_paper")).json")
			end

			# pCct = plot_L_contour(ωgrid, χgrid, C_mat_ctour)
			# savejson(pCct, pwd()*"/../Graphs/tests/Ccontour.json")			

			# print_save("\nCurrent L = $L against current min = $Lmin")

			if jω == 1
				save("../../ct_1_temp.jld", "ct", ct)
				save("../ct_1_temp.jld", "ct", ct)
			end


			if jχ == 1 && jω == 2 && flag
				save("../../ct_1.jld", "ct", ct)
				save("../ct_1.jld", "ct", ct)
			end

			if L < L_min
				L_min = L_mat_ctour[jω, jχ]
				ω_min = ωv
				χ_min = χv
				a_min = a_vec[jω]
				ja_min = ja

				save("../../ct_opt.jld", "ct", ct)
				ct_best.ω, ct_best.χ = ωv, χv
				ct_best.L, ct_best.gπ = ct.L, ct.gπ

				_, pL, pπ, _, pp, _ = makeplots_ct_pa(ct);
				savejson(pL, pwd()*"/../Graphs/tests/opt_L.json")
				savejson(pπ, pwd()*"/../Graphs/tests/opt_g.json")
				savejson(pp, pwd()*"/../Graphs/tests/opt_p.json")


				psim, pLsim = plot_simul(ct, T = 40, N = 50000, jp0 = 3)
				savejson(psim, pwd()*"/../Graphs/tests/simul_opt.json")
				savejson(pLsim,pwd()*"/../Graphs/tests/simul_Lopt.json")
			end
			if jω == length(ωgrid) && jχ == 1
				psim, pLsim = plot_simul(ct, T = 40, N = 50000, jp0 = 3)
				savejson(psim, pwd()*"/../Graphs/tests/simul_1.json")
				savejson(pLsim,pwd()*"/../Graphs/tests/simul_L1.json")
				_, pL, pπ, _, pp, _ = makeplots_ct_pa(ct, slides=true);
				savejson(pL, pwd()*"/../Graphs/tests/first_L_slides.json")
				savejson(pπ, pwd()*"/../Graphs/tests/first_g_slides.json")
				savejson(pp, pwd()*"/../Graphs/tests/first_p_slides.json")
				_, pL, pπ, _, pp, _ = makeplots_ct_pa(ct, slides=false);
				savejson(pL, pwd()*"/../Graphs/tests/first_L_paper.json")
				savejson(pπ, pwd()*"/../Graphs/tests/first_g_paper.json")
				savejson(pp, pwd()*"/../Graphs/tests/first_p_paper.json")
				save("../../first_ct.jld", "ct", ct)
			end

			for slides in [true, false]
				pCct = plot_L_contour(ωgrid, χgrid, C_mat[:,:,3,ja_min], name_y="C", slides=slides)
				savejson(pCct, pwd()*"/../Graphs/tests/Ccontour$(ifelse(slides, "_slides", "_paper")).json")
			end

		end

		s = "\nMinimum element is $(@sprintf("%.3g",Lmin)) with a₀ = $(@sprintf("%.3g", annualized(amin)))"
		# Optim.converged(res) ? s = s*" ✓" : nothing
		print_save(s)

		perm_order = sortperm(ω_vec)

		ω_vec = ω_vec[perm_order]
		L_vec = L_vec[perm_order]
		a_vec = a_vec[perm_order]

		new_L = scatter(;x=ω_vec, y=L_vec, name = "χ = $(@sprintf("%.3g",annualized(χv)))%")
		push!(Lplot, new_L)

		new_a = scatter(;x=ω_vec, y=annualized.(a_vec), name = "χ = $(@sprintf("%.3g",annualized(χv)))%")
		push!(aplot, new_a)

		#=
			if remote
				p1 = makeplots_ct_pa(ct)
				relayout!(p1, title="ω = $(@sprintf("%.3g",ct.ω))", width=1200, height=900)
				savejson(p1, pwd()*"/../Graphs/tests/summary_jom_$(jω).json")

				p2 = plot_simul(ct);
				savejson(p2, pwd()*"/../Graphs/tests/simul_jom_$(jω).json");
			end
		=#
	end

	print_save("\nWent through the spectrum of ω's in $(time_print(time()-t0))")
	print_save("\nOverall minimum announcement c = (a₀, ω, χ) = $(annualized(a_min)), $ω_min, $(annualized(χ_min))")

	for slides = [true, false]
		p1 = plot_plans_p(ct, L_mat, ωgrid, χgrid, slides=slides)
		savejson(p1, pwd()*"/../Graphs/tests/plans$(ifelse(slides, "_slides", "_paper")).json")
	end

	ν = ones(length(ωgrid), length(χgrid), length(ct_best.agrid))
	ν *= 1/sum(ν)
	mt = MultiType(ct_best, ωgrid, χgrid, ct.pgrid[3], ν, ν, L_mat, C_mat)

	return annualized(a_min), ω_min, annualized(χ_min), mt
end

Bayes_plan(ν, z, μ) = z*ν / (z*ν + (1-z)*μ)

function eval_k_to_mu(mt::MultiType, k, itp_L; get_mu::Bool=false)

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
			if disp > 1e-4
				print_save("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)")
			end
			pv = res.minimizer

			νv = mt.ν[jω, jχ, ja]
			res = Optim.optimize(
				μ -> (Bayes_plan(νv, mt.z, μ) - pv)^2, 0, 1, GoldenSection())
			disp = res.minimum
			if disp > 1e-4
				print_save("WARNING: Couldn't find p0 at state ($ωv, $χv, $av)")
			end
			μv = res.minimizer

			μ[jω, jχ, ja] = μv
		end
		p0[jω, jχ, ja] = pv
	end
	if get_mu 
		return μ
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
		print_save("WARNING: variance of Lᴺ = $(@sprintf("%0.3g",V))")
	end

	knots = (ωgrid[end:-1:1], χgrid, pgrid, agrid)
	itp_L = interpolate(knots, L_mat[end:-1:1,:,:,:], Gridded(Linear()))

	res = Optim.optimize(
		k -> (eval_k_to_mu(mt, k, itp_L)-1)^2,
		k_min, k_max, GoldenSection())

	if res.minimum > 1e-4
		print_save("WARNING: Couldn't find μ at z = $zv")
	end

	k_star = res.minimizer

	mt.μ = eval_k_to_mu(mt, k_star, itp_L; get_mu = true)

	return k_star
end

function mimic_z(mt::MultiType, N=50; decay::Bool=false)

	zgrid = cdf.(Beta(4,1), range(0,1,length=N))
	move_grids!(zgrid, xmax=0.9, xmin=mt.ct.pgrid[3])

	data = zeros(N,6)
	datanames = ["ω", "a", "χ", "s_ω", "s_a", "s_χ"]

	for (jz, zv) in enumerate(zgrid)
		find_equil!(mt, zv)
		data[jz,:] .= find_plan_μ(mt, decay=decay)
	end

	return data, datanames, zgrid
end
