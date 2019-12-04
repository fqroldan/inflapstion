using Distributed

# @everywhere 
using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD

# @everywhere 
include("type_def.jl")
# @everywhere 
include("reporting_routines.jl")
# @everywhere 
include("simul.jl")
# @everywhere 
include("plotting_routines.jl")

# @everywhere begin

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

function Bayes(ct::CrazyType, obs_π, exp_π, pv, av)

	numer = pv * pdf_ϵ(ct, obs_π - av)
	denomin = numer + (1.0-pv) * pdf_ϵ(ct, obs_π - exp_π)

	p′ = numer / denomin

	p′ = max(0.0, min(1.0, p′))

	if isapprox(denomin, 0.0)
		p′ = 0.0
	end
	# drift = (1.0 - pv) * 0.15
	# drift = -(pv) * 0.15

	return p′
end

PC(ct::CrazyType{Forward}, obs_π, πe, exp_π′) = (1/ct.κ) * (obs_π - ct.β * exp_π′)
PC(ct::CrazyType{Backward}, obs_π, πe, exp_π′) = 1/ct.κ  * (obs_π - πe)

function cond_Ldev(ct::CrazyType, itp_gπ, itp_L, obs_π, pv, av)
	aprime = ϕ(ct, av)

	πe = pv*av + (1-pv)*exp_π
	exp_π′ = itp_gπ(0.0, aprime)

	y = PC(ct, obs_π, πe, exp_π′) # Automatically uses method for forward or backward
	L′ = itp_L(0.0, aprime)

	L = (ct.ystar-y)^2 + ct.γ * obs_π^2 + ct.β * L′

	return L
end

function cond_L(ct::CrazyType, itp_gπ, itp_L, itp_C, obs_π, pv, av; get_y::Bool=false)
	exp_π  = itp_gπ(pv, av)
	if isapprox(pv, 0.0)
		pprime = 0.0
	elseif isapprox(pv, 1.0)
		pprime = 1.0
	else
		pprime = Bayes(ct, obs_π, exp_π, pv, av)
	end
	aprime = ϕ(ct, av)

	πe = pv*av + (1-pv)*exp_π

	L′ = itp_L(pprime, aprime)
	C′ = itp_C(pprime, aprime)
	exp_π′ = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)

	y = PC(ct, obs_π, πe, exp_π′) # Automatically uses method for forward or backward
	L = (ct.ystar-y)^2 + ct.γ * obs_π^2 + ct.β * L′
	if get_y
		return y, pprime, C′
	end
	return L
end

function exp_L(ct::CrazyType, itp_gπ, itp_L, itp_C, control_π, pv, av; get_y::Bool=false)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av) * pdf_ϵ(ct, ϵv)
	(val, err) = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)

	# sum_prob, err = hquadrature(x -> pdf_ϵ(ct, x), -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
	sum_prob = cdf_ϵ(ct, 3.09*ct.σ) - cdf_ϵ(ct, -3.09*ct.σ)

	val = val / sum_prob

	if get_y
		f_y(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av; get_y=true)[1] * pdf_ϵ(ct, ϵv)
		Ey, err = hquadrature(f_y, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
		f_p(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av; get_y=true)[2] * pdf_ϵ(ct, ϵv)
		Ep, err = hquadrature(f_p, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)
		f_C(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av; get_y=true)[3] * pdf_ϵ(ct, ϵv)
		Ec, err = hquadrature(f_p, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-10, atol=0, maxevals=0)

		Ey = Ey / sum_prob
		Ep = Ep / sum_prob
		Ec = Ec / sum_prob

		return Ey, Ep, Ec
	end

	return val
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, itp_C, π_guess, pv, av)

	minπ = max(0, π_guess - 3.09*ct.σ)
	maxπ = min(1.1*maximum(ct.agrid), π_guess + 3.09*ct.σ)
	
#=	res = Optim.optimize(
			gπ -> exp_L(ct, itp_gπ, itp_L, itp_C, gπ, pv, av),
			minπ, maxπ, GoldenSection()#, rel_tol=1e-20, abs_tol=1e-20, iterations=10000
			)
=#
	
	obj_f(x) = exp_L(ct, itp_gπ, itp_L, itp_C, x, pv, av)
	res = Optim.optimize(
		gπ -> obj_f(first(gπ)),
		[π_guess], LBFGS()#, autodiff=:forward#, Optim.Options(f_tol=1e-12)
		)

	gπ, L = first(res.minimizer), res.minimum

	if Optim.converged(res) == false
		# a = Optim.iterations(res)
		print_save("π∈ [$minπ, $maxπ]")
		# println(a)
		resb = Optim.optimize(
				gπ -> exp_L(ct, itp_gπ, itp_L, itp_C, gπ, pv, av),
				minπ, maxπ, Brent(), rel_tol=1e-12, abs_tol=1e-12#, iterations=100000
				)
		if resb.minimum < res.minimum
			gπ, L = resb.minimizer, resb.minimum
		end
	end

	return gπ, L
end

function optim_step(ct::CrazyType, itp_gπ, itp_L, itp_C, gπ_guess; optimize::Bool=true)
	# gπ, L  = SharedArray{Float64}(ct.gπ), SharedArray{Float64}(ct.L)
	# Ey, Eπ = SharedArray{Float64}(ct.Ey), SharedArray{Float64}(ct.Eπ)
	# Ep, C  = SharedArray{Float64}(ct.Ep), SharedArray{Float64}(ct.C)
	gπ, L  = zeros(size(ct.gπ)), zeros(size(ct.L))
	Ey, Eπ = zeros(size(ct.Ey)), zeros(size(ct.Eπ))
	Ep, C  = zeros(size(ct.Ep)), zeros(size(ct.C))
	πN 	   = Nash(ct)
	apgrid = gridmake(1:ct.Np, 1:ct.Na)
	Threads.@threads for js in 1:size(apgrid,1)
	# @sync @distributed  for js in 1:size(apgrid,1)
    # for js in 1:size(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.pgrid[jp], ct.agrid[ja]
		π_guess = gπ_guess[jp, ja]
		if optimize
			# π_guess = itp_gπ(pv, av)
			gπ[jp, ja], L[jp, ja] = opt_L(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av)
		else
			gπ[jp, ja] = π_guess
			L[jp, ja] = exp_L(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av)
		end
		Ey[jp, ja], Ep[jp, ja], EC′ = exp_L(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av; get_y=true)
		Eπ[jp, ja] = pv * av + (1.0 - pv) * gπ[jp, ja]

		if av >= πN || isapprox(av, πN)
			C[jp, ja] = (1-ct.β)*1 + ct.β * EC′
		else
			C[jp, ja] = (1-ct.β)*(πN - Eπ[jp,ja])/(πN-av) + ct.β * EC′
		end
	end

	return gπ, L, Ey, Eπ, Ep, C
end

function pf_iter(ct::CrazyType, Egπ, gπ_guess; optimize::Bool=true)
	#=	
	knots = (ct.pgrid[2:end], ct.agrid)
	itp_gπ_1 = interpolate(knots, Egπ[2:end,:],  Gridded(Linear()))
	itp_gπ_2 = extrapolate(itp_gπ_1, Flat())
	itp_L_1 = interpolate(knots, ct.L[2:end,:],  Gridded(Linear()))
	itp_L_2 = extrapolate(itp_L_1, Flat())

	η = 0.9
	plow = ct.pgrid[2] * η + ct.pgrid[1] * (1-η)
	gπ_lowp = [itp_gπ_2(plow, av) for (ja, av) in enumerate(ct.agrid)]
	L_lowp = [itp_L_2(plow, av) for (ja, av) in enumerate(ct.agrid)]

	pgrid_large = [ct.pgrid[1]; plow; ct.pgrid[2:end]]

	gπ_large = Array{Float64}(undef, ct.Np+1, ct.Na)
	L_large = Array{Float64}(undef, ct.Np+1, ct.Na)
	for jp in 1:ct.Np+1
		for (ja, av) in enumerate(ct.agrid)
			if jp > 2
				gπ_large[jp, ja] = Egπ[jp-1,ja]
				L_large[jp, ja] = ct.L[jp-1,ja]
			elseif jp == 1
				gπ_large[jp, ja] = gπ_lowp[ja]
				L_large[jp, ja] = L_lowp[ja]
			else
				gπ_large[jp, ja] = Egπ[1, ja]
				L_large[jp, ja] = ct.L[1, ja]
			end
		end
	end
	knots = (pgrid_large, ct.agrid)
	itp_gπ = interpolate(knots, gπ_large, Gridded(Linear()))
	itp_L  = interpolate(knots, L_large, Gridded(Linear()))
	=#
	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, Egπ, Gridded(Linear()))
	itp_L  = interpolate(knots, ct.L, Gridded(Linear()))
	itp_C  = interpolate(knots, ct.C, Gridded(Linear()))


	new_gπ, new_L, new_y, new_π, new_p, new_C = optim_step(ct, itp_gπ, itp_L, itp_C, gπ_guess; optimize=optimize)

	return new_gπ, new_L, new_y, new_π, new_p, new_C
end

function pfi!(ct::CrazyType, Egπ; tol::Float64=1e-12, maxiter::Int64=1000, verbose::Bool=true, reset_guess::Bool=false)
	dist = 10.
	iter = 0
	upd_η2 = 0.75

	rep = "\nStarting PFI (tol = $(@sprintf("%0.3g",tol)))"
	verbose ? print_save(rep,true) : print(rep)

    if reset_guess
	    ct.gπ = zeros(size(ct.gπ))
		ct.L = ones(ct.Np, ct.Na)
	end

	old_gπ = copy(Egπ)
	new_gπ = zeros(size(old_gπ))

	while dist > tol && iter < maxiter
		iter += 1

		for jj in 1:10
			_, new_L, _, _, _ = pf_iter(ct, Egπ, old_gπ; optimize=false)
			ct.L  = upd_η2 * new_L  + (1.0-upd_η2) * ct.L
		end

		old_L = copy(ct.L)

		new_gπ, new_L, new_y, new_π, new_p, new_C = pf_iter(ct, Egπ, old_gπ)

		dist = sqrt.(sum( (new_L  - old_L ).^2 )) / sqrt.(sum(old_L .^2))

		ct.L  = upd_η2 * new_L  + (1.0-upd_η2) * ct.L
		old_gπ = upd_η2 * new_gπ + (1.0-upd_η2) * old_gπ
		ct.Ey = new_y
		ct.Eπ = new_π
		ct.Ep = new_p
		ct.C  = new_C

		# if verbose && iter % 10 == 0
		# 	print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
		# end
	end

	dist2 = 10.
	iter2 = 0
	while dist > tol && iter2 < maxiter
		iter2 += 1
		old_C = copy(ct.C)
		_, _, _, _, _, new_C = pf_iter(ct, Egπ, old_gπ; optimize=false)
		dist2 = sqrt.(sum( (new_C  - old_C ).^2 )) / sqrt.(sum(old_C .^2))
		ct.C  = upd_η2 * new_C  + (1.0-upd_η2) * ct.C
	end

	if verbose && dist <= tol
		print("\nConverged in $iter iterations.")
	elseif verbose
		print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
	end

	return (dist <= tol), new_gπ
end

decay_η(ct::CrazyType, η) = max(0.8*η, 5e-5)

function Epfi!(ct::CrazyType; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, tempplots::Bool=false, upd_η::Float64=0.01, switch_η = 40)
	dist = 10.
	iter = 0
	
	print_save("\nRun with ω = $(@sprintf("%.3g",ct.ω)), χ = $(@sprintf("%.3g",annualized(ct.χ)))% at $(Dates.format(now(), "HH:MM"))")

	dists = []

	reset_guess = false
	tol_pfi = 1e-5 / 0.99
	while dist > tol && iter < maxiter
		iter += 1
		tol_pfi = max(tol_pfi*0.9, 1e-8)

		old_gπ, old_L = copy(ct.gπ), copy(ct.L);

		flag, new_gπ = pfi!(ct, old_gπ; verbose=verbose, reset_guess=reset_guess, tol=tol_pfi);
		reset_guess = !flag

		dist = sqrt.(sum( (new_gπ  - ct.gπ ).^2 )) / sqrt.(sum(ct.gπ .^2))
		push!(dists, dist)
		rep_status = "\nAfter $iter iterations, d(π) = $(@sprintf("%0.3g",dist))"
		if flag
			rep_status *= "✓ "
		end
		if verbose #&& iter % 10 == 0
			print_save(rep_status*"\n", true)
		else
			print(rep_status)
		end

		ct.gπ = upd_η * new_gπ + (1.0-upd_η) * ct.gπ;

		if tempplots && (iter % 5 == 0 || dist <= tol)
			p1, pL, pE, pC, pp = makeplots_ct_pa(ct);
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
	if verbose && dist <= tol
		print_save("\nConverged in $iter iterations.",true)
	elseif verbose
		print_save("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))",true)
	end
	p1, pL, pπ, pC, pp = makeplots_ct_pa(ct);
	savejson(pC, pwd()*"/../Graphs/tests/tempC.json")
	savejson(pπ, pwd()*"/../Graphs/tests/tempg.json")
	
	return dist
end

# function choose_ω!(L_mat, ct::CrazyType{Forward}, Nω=size(L_mat,1); remote::Bool=true, upd_η=0.1)
# 	choose_ω!(L_mat, ct, Forward, Nω; remote=remote, upd_η=upd_η)
# end

# function choose_ω!(L_mat, ct::CrazyType{Backward}, Nω=size(L_mat,1); remote::Bool=true, upd_η=0.1)
# 	choose_ω!(L_mat, ct, Backward, Nω; remote=remote, upd_η=upd_η)
# end

function choose_ω!(L_mat, ct::CrazyType, Nω=size(L_mat,1); upd_η=0.1)

	T = which_PC(ct)

	if T == Backward
		ωmax = 3.0
	elseif T == Forward
		ωmax = 1.0
	end
	ωgrid = cdf.(Beta(1,1), range(1,0,length=Nω))
	move_grids!(ωgrid, xmax = ωmax, xmin = 0.01)

	Nχ = size(L_mat, 2)
	χgrid = range(0.0, 0.5*Nash(ct), length = Nχ)

	print_save("\nLooping over behavioral types with ω ∈ [$(minimum(ωgrid)), $(maximum(ωgrid))]")
	print_save("\n")

	L_min = 100.
	ω_min = 1.0
	χ_min = 1.0
	a_min = 1.0
	t0 = time()
	Lplot = []
	aplot = []
	L_mat_ctour = zeros(Nω, Nχ) * NaN
	C_mat_ctour = zeros(Nω, Nχ) * NaN
	Lmin = 1e8
	for (jχ, χv) in enumerate(χgrid)
		L_vec = []
		a_vec = []
		ω_vec = []

		function wrap_Epfi!(ct::CrazyType, ωv, L_vec, a_vec, ω_vec, Lplot, L_mat_save, aplot)
			ct.ω = ωv

			t1 = time()
			tol = 90e-4
			# if length(L_vec) > 0
			# 	upd_η = 0.005
			# end
			dist = Epfi!(ct, verbose = true, tol=tol, tempplots=false, upd_η=upd_η)
			write(pwd()*"/../temp.txt", "")
			
			flag = (dist <= tol)
			Lmin, ja = findmin(ct.L[3,:])
			# Cmin = ct.C[3,ja]
			Cmin = ct.C[3,end]
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

			return Lmin, Cmin
		end

		ωmin = 1e8
		amin = 1e8
		for (jω, ωv) in enumerate(ωgrid)
			old_L, old_gπ = copy(ct.L), copy(ct.gπ)
			ct = CrazyType(T; χ = χv)
			ct.L, ct.gπ = old_L, old_gπ
			
			L_mat_save = zeros(ct.Np, ct.Na)
			L, C = wrap_Epfi!(ct, ωv, L_vec, a_vec, ω_vec, Lplot, L_mat_save, aplot)

			L_mat[jω, jχ, :, :] = L_mat_save
			L_mat_ctour[jω, jχ] = L

			C_mat_ctour[jω, jχ] = C 

			pLct = plot_L_contour(ωgrid, χgrid, L_mat_ctour)
			savejson(pLct, pwd()*"/../Graphs/tests/contour.json")

			pCct = plot_L_contour(ωgrid, χgrid, C_mat_ctour)
			savejson(pCct, pwd()*"/../Graphs/tests/Ccontour.json")			

			# print_save("\nCurrent L = $L against current min = $Lmin")

			if jχ == 1 && jω == 2
				save("../../ct_1.jld", "ct", ct)
			end

			if L < L_min
				L_min = L_mat_ctour[jω, jχ]
				ω_min = ωv
				χ_min = χv
				a_min = a_vec[jω]

				save("../../ct_opt.jld", "ct", ct)

				_, pL, pπ, _, pp = makeplots_ct_pa(ct);
				savejson(pL, pwd()*"/../Graphs/tests/opt_L.json")
				savejson(pπ, pwd()*"/../Graphs/tests/opt_g.json")
				savejson(pp, pwd()*"/../Graphs/tests/opt_p.json")


				psim, pLsim = plot_simul(ct, T = 40, N = 50000, jp0 = 3)
				savejson(psim, pwd()*"/../Graphs/tests/simul_opt.json")
				savejson(pLsim,pwd()*"/../Graphs/tests/simul_Lopt.json")
			end
			if jω == 1 && jχ == 1
				psim, pLsim = plot_simul(ct, T = 40, N = 50000, jp0 = 3)
				savejson(psim, pwd()*"/../Graphs/tests/simul_1.json")
				savejson(pLsim,pwd()*"/../Graphs/tests/simul_L1.json")
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

	p1 = plot_plans_p(ct, L_mat, ωgrid, χgrid)

	nothing
end
# end # everywhere

# ct = CrazyType(; ω = ωmin)
# Epfi!(ct);


# ct = CrazyType(ω = 0.0125);
# Epfi!(ct, tempplots=true)
# p2 = plot_simul(ct, noshocks=true)
# if remote
# 	savejson(p2, pwd()*"/../Graphs/tests/simul.json")
# end

# p1 = makeplots_ct_pa(ct);
# p1



# using JLD
# save("ct.jld", "ct", ct)

# plot_simul(ct, T=50)
