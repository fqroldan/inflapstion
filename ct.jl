using Distributed

@everywhere using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays, Dates

@everywhere begin
function Bayes(ct::CrazyType, obs_π, exp_π, pv, av)

	numer = pv * pdf_ϵ(ct, obs_π - av)
	denomin = numer + (1.0-pv) * pdf_ϵ(ct, obs_π - exp_π)

	p′ = numer / denomin

	p′ = max(0.0, min(1.0, p′))

	if isapprox(denomin, 0.0) || isapprox(numer, 0.0)
		p′ = 0.0
	end
	# drift = (1.0 - pv) * 0.15
	# drift = -(pv) * 0.15

	return p′
end

NKPC(ct::CrazyType, obs_π, exp_π′) = (1.0/ct.κ) * (obs_π - ct.β * exp_π′)
# BLPC(ct::CrazyType, obs_π, exp_π)  = ct.κ * (obs_π - exp_π)

function cond_L(ct::CrazyType, itp_gπ, itp_L, obs_π, pv, av; get_y::Bool=false)
	exp_π  = itp_gπ(pv, av)
	if isapprox(pv, 0.0)
		pprime = 0.0
	elseif isapprox(pv, 1.0)
		pprime = 1.0
	else
		pprime = Bayes(ct, obs_π, exp_π, pv, av)
	end
	aprime = ϕ(ct, av)

	#=
	σ_η = 0.05
	η_vec = range(-1.96*σ_η, 1.96*σ_η, length = 9)
	pη = pdf.(Normal(0,σ_η), η_vec)
	pη = pη / sum(pη)

	ap_vec = aprime .* (1.0 .+ η_vec)
	L′ = 0.0
	for (jap, apv) in enumerate(ap_vec)
		apv = max(min(apv, maximum(ct.agrid)), minimum(ct.agrid))
		L′ += itp_L(pprime, apv)# * pη[jap]
	end
	=#

	L′ = itp_L(pprime, aprime)
	exp_π′ = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)

	y = NKPC(ct, obs_π, exp_π′)
	L = (ct.ystar-y)^2 + ct.γ * obs_π^2 + ct.β * L′
	if get_y
		return y, pprime
	end
	return L
end

function exp_L(ct::CrazyType, itp_gπ, itp_L, control_π, pv, av; get_y::Bool=false)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, pv, av) * pdf_ϵ(ct, ϵv)
	(val, err) = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-12, atol=0, maxevals=0)

	# sum_prob, err = hquadrature(x -> pdf_ϵ(ct, x), -3.09*ct.σ, 3.09*ct.σ, rtol=1e-12, atol=0, maxevals=0)
	sum_prob = cdf_ϵ(ct, 3.09*ct.σ) - cdf_ϵ(ct, -3.09*ct.σ)

	val = val / sum_prob

	if get_y
		f_y(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, pv, av; get_y=true)[1] * pdf_ϵ(ct, ϵv)
		Ey, err = hquadrature(f_y, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-12, atol=0, maxevals=0)
		f_p(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, pv, av; get_y=true)[2] * pdf_ϵ(ct, ϵv)
		Ep, err = hquadrature(f_p, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-12, atol=0, maxevals=0)

		Ey = Ey / sum_prob
		Ep = Ep / sum_prob

		return Ey, Ep
	end

	return val
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, π_guess, pv, av)

	minπ, maxπ = -0.25, 1.1*maximum(ct.agrid)
	#=
	res = Optim.optimize(
			gπ -> exp_L(ct, itp_gπ, itp_L, gπ, pv, av),
			minπ, maxπ, Brent()#, rel_tol=1e-20, abs_tol=1e-20, iterations=10000
			)
	=#
	obj_f(x) = exp_L(ct, itp_gπ, itp_L, x, pv, av)
	res = Optim.optimize(
		gπ -> obj_f(first(gπ)),
		[π_guess], LBFGS(), Optim.Options(f_tol=1e-12)
		)

	gπ, L = first(res.minimizer), res.minimum

	if Optim.converged(res) == false
		# a = Optim.iterations(res)
		# println(a)
		resb = Optim.optimize(
				gπ -> exp_L(ct, itp_gπ, itp_L, gπ, pv, av),
				minπ, maxπ, Brent(), rel_tol=1e-18, abs_tol=1e-18#, iterations=100000
				)
		if resb.minimum < res.minimum
			gπ, L = resb.minimizer, resb.minimum
		end
	end
	return gπ, L
end

function optim_step(ct::CrazyType, itp_gπ, itp_L, gπ_guess; optimize::Bool=true)
	gπ, L  = SharedArray{Float64}(ct.gπ), SharedArray{Float64}(ct.L)
	Ey, Eπ = SharedArray{Float64}(ct.Ey), SharedArray{Float64}(ct.Eπ)
	Ep 	   = SharedArray{Float64}(ct.Ep)
	# gπ, L = Array{Float64}(undef, size(ct.gπ)), Array{Float64}(undef, size(ct.L))
	apgrid = gridmake(1:ct.Np, 1:ct.Na)
	@sync @distributed for js in 1:size(apgrid,1)
    # for js in 1:size(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.pgrid[jp], ct.agrid[ja]
		π_guess = gπ_guess[jp, ja]
		if optimize
			# π_guess = itp_gπ(pv, av)
			gπ[jp, ja], L[jp, ja] = opt_L(ct, itp_gπ, itp_L, π_guess, pv, av)
		else
			gπ[jp, ja] = π_guess
			L[jp, ja] = exp_L(ct, itp_gπ, itp_L, π_guess, pv, av)
		end
		Ey[jp, ja], Ep[jp, ja] = exp_L(ct, itp_gπ, itp_L, π_guess, pv, av; get_y=true)
		Eπ[jp, ja] = pv * av + (1.0 - pv) * gπ[jp, ja]
	end

	return gπ, L, Ey, Eπ, Ep
end

function pf_iter(ct::CrazyType, Egπ, gπ_guess; optimize::Bool=true)
	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, Egπ,  Gridded(Linear()))
	itp_L  = interpolate(knots, ct.L, Gridded(Linear()))

	# itp = interpolate(ct.L, BSpline(Cubic(Line(OnGrid()))))
	# itp_L = Interpolations.scale(itp, ct.pgrid, ct.agrid)
	# itp = interpolate(ct.gπ, BSpline(Cubic(Line(OnGrid()))))
	# itp_gπ = Interpolations.scale(itp, ct.pgrid, ct.agrid)

	new_gπ, new_L, new_y, new_π, new_p = optim_step(ct, itp_gπ, itp_L, gπ_guess; optimize=optimize)

	return new_gπ, new_L, new_y, new_π, new_p
end

function pfi!(ct::CrazyType, Egπ; tol::Float64=1e-12, maxiter::Int64=1000, verbose::Bool=true, reset_guess::Bool=false)
	dist = 10.
	iter = 0
	upd_η = 0.75

	rep = "\nStarting PFI (tol = $(@sprintf("%0.3g",tol)))"
	verbose ? print_save(rep) : print(rep)

    if reset_guess
	    ct.gπ = zeros(size(ct.gπ))
		ct.L = ones(ct.Np, ct.Na)
	end
	while dist > tol && iter < maxiter
		iter += 1

		for jj in 1:5
			_, new_L, _, _, _ = pf_iter(ct, Egπ, ct.gπ; optimize=false)
			ct.L  = upd_η * new_L  + (1.0-upd_η) * ct.L
		end

		old_gπ, old_L = copy(ct.gπ), copy(ct.L)

		new_gπ, new_L, new_y, new_π, new_p = pf_iter(ct, Egπ, ct.gπ)

		dist = sqrt.(sum( (new_L  - old_L ).^2 )) / sqrt.(sum(old_L .^2))

		ct.L  = upd_η * new_L  + (1.0-upd_η) * ct.L
		ct.gπ = upd_η * new_gπ + (1.0-upd_η) * ct.gπ
		ct.Ey = new_y
		ct.Eπ = new_π
		ct.Ep = new_p

		if verbose && iter % 10 == 0
			print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
		end
	end
	if verbose && dist <= tol
		print("\nConverged in $iter iterations.")
	elseif verbose
		print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
	end
	return (dist <= tol)
end

function Epfi!(ct::CrazyType; tol::Float64=5e-4, maxiter::Int64=2000, verbose::Bool=true, tempplots::Bool=false)
	dist = 10.
	iter = 0
	upd_η = 0.025
	
	print_save("\nStarting run with ω = $(@sprintf("%.3g",ct.ω)) at $(Dates.format(now(), "HH:MM"))")

	dists = []

	reset_guess = false
	tol_pfi = 1e-8 / 0.99
	while dist > tol && iter < maxiter
		iter += 1
		tol_pfi = max(min(tol_pfi*0.99, dist * 1e-7), 1e-12)

		old_gπ, old_L = copy(ct.gπ), copy(ct.L);

		flag = pfi!(ct, old_gπ; verbose = verbose, reset_guess=reset_guess, tol=tol_pfi);
		reset_guess = !flag

		dist = sqrt.(sum( (ct.gπ  - old_gπ ).^2 )) / sqrt.(sum(old_gπ .^2))
		push!(dists, dist)
		rep_status = "\nAfter $iter iterations, d(π) = $(@sprintf("%0.3g",dist))"
		if flag
			rep_status *= "✓ "
		end
		if verbose #&& iter % 10 == 0
			print_save(rep_status*"\n")
		else
			print(rep_status)
		end

		ct.gπ = upd_η * ct.gπ + (1.0-upd_η) * old_gπ;

		if tempplots
			p1 = makeplots_ct_pa(ct);
			relayout!(p1, title="iter = $iter")
			savejson(p1, pwd()*"/../Graphs/tests/temp.json")
			p2 = makeplot_conv(dists);
			savejson(p2, pwd()*"/../Graphs/tests/tempconv.json")
		end

	end
	if verbose && dist <= tol
		print("\nConverged in $iter iterations.")
	elseif verbose
		print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
	end

	return dist
end

function choose_ω!(L_mat, initial_guess::CrazyType=CrazyType(), Nω=size(L_mat,1); remote::Bool=true)

	ωgrid = cdf.(Beta(5,1), range(1,0,length=Nω))
	move_grids!(ωgrid, xmax = 0.2)

	ct = initial_guess

	print_save("\nLooping over behavioral types with ω ∈ [$(minimum(ωgrid)), $(maximum(ωgrid))]")
	print_save("\n")

	L_min = 100.
	ωmin = 1.0
	amin_min = 1.0
	t0 = time()
	jamin_vec = Vector{Int64}(undef, Nω)
	for (jω, ωv) in enumerate(ωgrid)
		Lguess, πguess = ct.L, ct.gπ
		ct = CrazyType(; ω = ωv)
		ct.L = Lguess
		ct.gπ = πguess
		t1 = time()
		tol = 5e-4
		dist = Epfi!(ct, verbose = false, tol=tol, tempplots=true)
		flag = (dist <= tol)

		# p1, p2, p3 = makeplots_ct(ct)
		# if remote
		# 	savejson(p1, pwd()*"/../Graphs/tests/ct_1_jomega_$(jω).json")
		# 	savejson(p2, pwd()*"/../Graphs/tests/ct_2_jomega_$(jω).json")
		# 	savejson(p3, pwd()*"/../Graphs/tests/ct_3_jomega_$(jω).json")
		# end

		if remote
			p1 = makeplots_ct_pa(ct)
			relayout!(p1, title="ω = $(@sprintf("%.3g",ωv))")
			savejson(p1, pwd()*"/../Graphs/tests/summary_jom_$(jω).json")

			p2 = plot_simul(ct)
			savejson(p2, pwd()*"/../Graphs/tests/simul_jom_$(jω).json")
		end

		# Save the corresponding value function
		L_mat[jω, :, :] = ct.L[:, :]
		lmin, ja = findmin(ct.L[2,:])
		amin = ct.agrid[ja]
		jamin_vec[jω] = ja
		print_save(": done in $(time_print(time()-t1))")
		s = "\nMinimum element is $(@sprintf("%.3g",lmin)) with a₀ = $(@sprintf("%.3g", amin))"
		flag ? s = s*" ✓" : nothing
		print_save(s)
		if lmin < L_min
			L_min = lmin
			ωmin = ωv
			amin_min = amin
		end

		if remote
			Lplot = [L_mat[jj, 2, jamin_vec[jj]] for jj in 1:jω]
			p3 = plot([
				scatter(;x=ωgrid[1:jω], y=Lplot)
				])
			savejson(p3, pwd()*"/../Graphs/tests/Loss_omega.json")
		end
	end

	print_save("\nWent through the spectrum of ω's in $(time_print(time()-t0))")
	Lplot = [L_mat[jj, 2, jamin_vec[jj]] for jj in 1:Nω]
	p1 = plot([
		scatter(;x=ωgrid, y=Lplot)
		])

	print_save("\nOverall minimum announcement a₀ = $amin_min with ω = $ωmin")

	return ωmin, p1
end
end # everywhere

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
