using Distributed

@everywhere using Distributions, Interpolations, Optim, HCubature, QuantEcon, LaTeXStrings, Printf, PlotlyJS, Distributed, SharedArrays
@everywhere include("reporting_routines.jl")

@everywhere begin

mutable struct CrazyType
	β::Float64
	γ::Float64
	α::Float64
	σ::Float64
	ystar::Float64
	ω::Float64

	pgrid::Vector{Float64}
	agrid::Vector{Float64}

	Np::Int64
	Na::Int64

	gπ::Array{Float64, 2}
	L::Array{Float64, 2}
	
	Ey::Array{Float64, 2}
	Eπ::Array{Float64, 2}
	Ep::Array{Float64, 2}
end
function CrazyType(;
		β = 0.96,
		γ = 60.0,
		α = 0.17,
		# α = 0.8,
		# α = 0.02,
		σ = 0.01,
		ystar = 0.1,
		ω = 0.271,
		# ω = 0.05,
		# ω = 0.1,
		Np = 90,
		Na = 60
		)

	y = 1.0 / (1.0+2.0*γ*α^2) * ystar - 2.0*γ*α / (1.0+2.0*γ*α^2) * 0.0

	A = α / (1.0 - β + α^2*γ) * ystar

	curv = 0.25
	pgrid = range(0, 1, length=Np).^(1.0/curv)
	curv = 0.5
	agrid = range(0, (1.1*A)^curv, length=Na).^(1.0/curv)

	gπ = zeros(Np, Na)
	L = ones(Np, Na)
	for (jp, pv) in enumerate(pgrid)
		for (ja, av) in enumerate(agrid)
			y_Nash = 1.0 / (1.0+2.0*γ*α^2) * ystar - 2.0*γ*α / (1.0+2.0*γ*α^2) * av
			# gπ[jp, ja] = A
			# yN = (1.0/α) * (gπ[jp, ja] - ct.β * exp(-ω) * A)
			# L[jp, ja] = (1.0/(1.0-β)) * ((ystar-yN)^2 + γ*gπ[jp,ja]^2)
		end
	end

	Ey = zeros(Np, Na)
	Eπ = zeros(Np, Na)
	Ep = zeros(Np, Na)

	return CrazyType(β, γ, α, σ, ystar, ω, pgrid, agrid, Np, Na, gπ, L, Ey, Eπ, Ep)
end

ϕ(ct::CrazyType, a::Float64) = exp(-ct.ω) * a

dist_ϵ(ct) = Normal(0, ct.σ)
pdf_ϵ(ct, ϵv) = pdf.(dist_ϵ(ct), ϵv)
cdf_ϵ(ct, ϵv) = cdf.(dist_ϵ(ct), ϵv)

function Bayes(ct::CrazyType, obs_π, exp_π, pv, av)

	numer = pv * pdf_ϵ(ct, obs_π - av)
	denomin = numer + (1.0-pv) * pdf_ϵ(ct, obs_π - exp_π)

	# drift = (1.0 - pv) * 0.15
	# drift = -(pv) * 0.15

	return numer / denomin# + drift
end

NKPC(ct::CrazyType, obs_π, exp_π′) = (1.0/ct.α) * (obs_π - ct.β * exp_π′)
# BLPC(ct::CrazyType, obs_π, exp_π)  = ct.α * (obs_π - exp_π)

function cond_L(ct::CrazyType, itp_gπ, itp_L, obs_π, pv, av; get_y::Bool=false)
	exp_π  = itp_gπ(pv, av)
	if isapprox(pv, 0.0)
		pprime = 0.0
	elseif isapprox(pv, 1.0)
		pprime = 1.0
	else
		pprime = Bayes(ct, obs_π, exp_π, pv, av)
	end

#=	σ_η = 0.05
	η_vec = range(-1.96*σ_η, 1.96*σ_η, length = 9)
	pη = pdf.(Normal(0,σ_η), η_vec)
	pη = pη / sum(pη)
=#
	aprime = ϕ(ct, av)

	#=ap_vec = aprime# .* (1.0 .+ η_vec)

	L′ = 0.0
	for (jap, apv) in enumerate(ap_vec)
		apv = max(min(apv, maximum(ct.agrid)), minimum(ct.agrid))
		L′ += itp_L(pprime, apv)# * pη[jap]
	end
	=#

	L′ = itp_L(pprime, aprime)

	exp_π′ = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)

	#=
	gπ′ = itp_gπ(pv, aprime)
	exp_π′ = pv * aprime + (1.0-pv) * gπ′
	=#

	y = NKPC(ct, obs_π, exp_π′)
	# y = BLPC(ct, obs_π, exp_π)

	L = (ct.ystar-y)^2 + ct.γ * obs_π^2 + ct.β * L′
	if get_y
		return y, pprime
	end
	return L
end

function exp_L(ct::CrazyType, itp_gπ, itp_L, control_π, pv, av; get_y::Bool=false)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, pv, av) * pdf_ϵ(ct, ϵv)
	(val, err) = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-32, atol=0, maxevals=0)

	# sum_prob, err = hquadrature(x -> pdf_ϵ(ct, x), -3.09*ct.σ, 3.09*ct.σ, rtol=1e-32, atol=0, maxevals=0)
	sum_prob = cdf_ϵ(ct, 3.09*ct.σ) - cdf_ϵ(ct, -3.09*ct.σ)

	val = val / sum_prob

	if get_y
		f_y(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, pv, av; get_y=true)[1] * pdf_ϵ(ct, ϵv)
		Ey, err = hquadrature(f_y, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-32, atol=0, maxevals=0)
		f_p(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, pv, av; get_y=true)[2] * pdf_ϵ(ct, ϵv)
		Ep, err = hquadrature(f_p, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-32, atol=0, maxevals=0)

		Ey = Ey / sum_prob
		Ep = Ep / sum_prob

		return Ey, Ep
	end

	return val
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, π_guess, pv, av)

	minπ, maxπ = -0.25, 1.1*maximum(ct.agrid)
	res = Optim.optimize(
			gπ -> exp_L(ct, itp_gπ, itp_L, gπ, pv, av),
			minπ, maxπ, GoldenSection(), rel_tol=1e-20, abs_tol=1e-20, iterations=10000
			)
	gπ = res.minimizer
	L = res.minimum

	if Optim.converged(res) == false
		# a = Optim.iterations(res)
		# println(a)
		resb = Optim.optimize(
				gπ -> exp_L(ct, itp_gπ, itp_L, gπ, pv, av),
				minπ, maxπ, Brent(), rel_tol=1e-18, abs_tol=1e-18#, iterations=100000
				)
		if resb.minimum < res.minimum
			gπ = resb.minimizer
			L = resb.minimum
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

function pfi!(ct::CrazyType, Egπ; tol::Float64=1e-12, maxiter::Int64=150, verbose::Bool=true, reset_guess::Bool=false)
	dist = 10.
	iter = 0
	upd_η = 0.75
    if verbose
        print_save("\nStarting PFI")
    end

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
			print_save("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
		end
	end
	if verbose && dist <= tol
		print_save("\nConverged in $iter iterations.")
	elseif verbose
		print_save("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
	end
	return (dist <= tol)
end

function Epfi!(ct::CrazyType; tol::Float64=1e-6, maxiter::Int64=100, verbose::Bool=true)
	dist = 10.
	iter = 0
	upd_η = 0.5

	reset_guess = false
	while dist > tol && iter < maxiter
		iter += 1

		old_gπ, old_L = copy(ct.gπ), copy(ct.L);

		flag = pfi!(ct, old_gπ; verbose = verbose, reset_guess=reset_guess);
		reset_guess = !flag

		dist = sqrt.(sum( (ct.gπ  - old_gπ ).^2 )) / sqrt.(sum(old_gπ .^2))
		if verbose #&& iter % 10 == 0
			print_save("\nAfter $iter iterations, d(π) = $(@sprintf("%0.3g",dist))\n")
		end

		ct.gπ = upd_η * ct.gπ + (1.0-upd_η) * old_gπ;

		# makeplots_ct(ct; make_png=true, tempplot=true)
	end
	nothing
end

function plot_ct(ct::CrazyType, y_tuple, n_tuple; make_pdf::Bool=false, make_png::Bool=false, tempplot::Bool=true)

	if length(y_tuple) != length(n_tuple)
		throw(error("Make sure # y's = # n's"))
	end

	col = [	"#1f77b4",  # muted blue
		"#ff7f0e",  # safety orange
		"#2ca02c",  # cooked asparagus green
		"#d62728",  # brick red
		"#9467bd",  # muted purple
		"#8c564b",  # chestnut brown
		"#e377c2",  # raspberry yogurt pink
		"#7f7f7f",  # middle gray
		"#bcbd22",  # curry yellow-green
		"#17becf"   # blue-teal
		]

	function lines(ct::CrazyType, y_mat; dim::Int64=0, title::String="", showleg::Bool=false)
		if dim == 1
			xgrid = ct.pgrid
			zgrid = ct.agrid
			xtitle= "𝑝"
		elseif dim == 2
			xgrid = ct.agrid
			zgrid = ct.pgrid
			xtitle= "𝑎"
		else
			throw(error("wrong dim"))
		end
		Nz = length(zgrid)
		l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, Nz)
		for (jz, zv) in enumerate(zgrid)
			if dim == 1
				y_vec = y_mat[:, jz]
				name = "𝑎"
			elseif dim == 2
				y_vec = y_mat[jz, :]
				name = "𝑝"
			end
			name = name * " = $(@sprintf("%.2g", zv))"
			jz % 2 == 0 ? showleg_i = showleg : showleg_i = false
			l_new = scatter(;x=xgrid, y=y_vec, name = name, showlegend = showleg_i, marker_color=col[ceil(Int,10*jz/Nz)])
			l[jz] = l_new
		end
		p = plot([l[jz] for jz in 1:Nz], Layout(;title=title, xaxis_title=xtitle))
		return p
	end
	N = length(y_tuple)
	pl = Array{PlotlyJS.SyncPlot,2}(undef, N, 2)
	for jj in 1:N
		pl[jj, 1] = lines(ct, y_tuple[jj], dim = 1, title=n_tuple[jj], showleg = (jj==1))
		pl[jj, 2] = lines(ct, y_tuple[jj], dim = 2, title=n_tuple[jj], showleg = (jj==1))
	end
	p1a = lines(ct, y_tuple[1], dim = 1, title=n_tuple[1], showleg = true)
	p1p = lines(ct, y_tuple[1], dim = 2, title=n_tuple[1], showleg = true)
	p2a = lines(ct, y_tuple[2], dim = 1, title=n_tuple[2])
	p2p = lines(ct, y_tuple[2], dim = 2, title=n_tuple[2])

	p = [p1a p1p; p2a p2p]

	if N == 2
		p = [pl[1,1] pl[1,2]; pl[2,1] pl[2,2]]
	elseif N == 3
		p = [pl[1,1] pl[1,2]; pl[2,1] pl[2,2]; pl[3,1] pl[3,2]]
	end

	relayout!(p, font_family = "Fira Sans Light", font_size = 12, height = 600, width = 950)

	function makeplot(p, ext::String, tempplot::Bool)
		tempplot ? add_temp = "_temp" : add_temp = ""
		savefig(p, pwd() * "/../Graphs/ct" * ext)
	end

	if make_pdf
		makeplot(p, ".pdf", tempplot)
	end
	if make_png
		makeplot(p, ".png", tempplot)
	end

	return p
end

end # everywhere

function makeplots_ct(ct::CrazyType; make_pdf::Bool=false, make_png::Bool=false, tempplot::Bool=true)

	gπ_over_a = zeros(size(ct.gπ))
	Ep_over_p = zeros(size(ct.Ep))
	for (jp, pv) in enumerate(ct.pgrid), (ja,av) in enumerate(ct.agrid)
		gπ_over_a[jp, ja] = ct.gπ[jp, ja] - av
		Ep_over_p[jp, ja] = ct.Ep[jp, ja] - pv
	end

	p1 = plot_ct(ct, (ct.gπ, ct.L), ("gπ", "𝓛"); make_pdf=make_pdf, make_png=make_png, tempplot=tempplot)

	p2 = plot_ct(ct, (ct.Ey, ct.Eπ), ("𝔼y", "𝔼π"); make_pdf=make_pdf, make_png=make_png, tempplot=tempplot)

	p3 = plot_ct(ct, (gπ_over_a, Ep_over_p), ("gπ-a", "𝔼p'-p"); make_pdf=make_pdf, make_png=make_png, tempplot=tempplot)

	return p1, p2, p3
end


function choose_ω()
	Nω = 25
	ωgrid = range(0.0, 0.75, length=Nω)

	ct = CrazyType()
	L_mat = zeros(Nω, ct.Np, ct.Na)

	L_min = 100.
	ωmin = 1.0
	amin_min = 1.0

	jamin_vec = Vector{Int64}(undef, Nω)
	for (jω, ωv) in enumerate(ωgrid)
		Lguess, πguess = ct.L, ct.gπ
		ct = CrazyType(; ω = ωv)
		# ct.L = Lguess
		# ct.gπ = πguess

		flag = Epfi!(ct, verbose = false)

		# Save the corresponding value function
		L_mat[jω, :, :] = ct.L[:, :]
		lmin, ja = findmin(ct.L[2,:])
		amin = ct.agrid[ja]
		jamin_vec[jω] = ja
		s = "\nMinimum element at ω = $(@sprintf("%.3g",ωv)) is $(@sprintf("%.3g",lmin)) with a₀ = $(@sprintf("%.3g", amin))"
		flag ? s = s*" ✓" : nothing
		print_save(s)
		if lmin < L_min
			L_min = lmin
			ωmin = ωv
			amin_min = amin
		end
	end

	Lplot = [L_mat[jj, 2, jamin_vec[jj]] for jj in 1:Nω]
	p1 = plot([
		scatter(;x=ωgrid, y=Lplot)
		])

	print_save("\nOverall minimum announcement a₀ = $amin_min with ω = $ωmin")

	return L_mat, ωmin, p1
end

function iter_simul(ct::CrazyType, itp_gπ, pv, av)
	ϵ = rand(dist_ϵ(ct))

	exp_π = itp_gπ(pv, av)
	obs_π = exp_π+ϵ
	pprime = Bayes(ct, obs_π, exp_π, pv, av)
	aprime = ϕ(ct, av)
	gπ′ = itp_gπ(pprime, aprime)
	exp_π′ = pprime * aprime + (1.0-pprime) * gπ′

	y = NKPC(ct, obs_π, exp_π′)

	return pprime, aprime, obs_π, y
end

function simul(ct::CrazyType; T::Int64=50, jp0::Int64=2)
	p0 = ct.pgrid[jp0]

	_, ind_a0 = findmin(ct.L[jp0, :])
	a0 = ct.agrid[ind_a0]

	p, a = p0, a0

	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()))

	p_vec, a_vec, π_vec, y_vec = zeros(T), zeros(T), zeros(T), zeros(T)
	for tt = 1:T
		p_vec[tt], a_vec[tt] = p, a
		pp, ap, πt, yt = iter_simul(ct, itp_gπ, p, a)
		π_vec[tt], y_vec[tt] = πt, yt
		if tt == T
			p_vec[end], a_vec[end] = pp, ap
		end
		p, a = pp, ap
	end

	return p_vec, a_vec, π_vec, y_vec
end

function plot_simul(ct::CrazyType; T::Int64=50, jp0::Int64=2)
	p_vec, a_vec, π_vec, y_vec = simul(ct, T=T, jp0=jp0)

	pp = plot(scatter(;x=1:T, y=p_vec), Layout(;title="Reputation"))
	pa = plot(scatter(;x=1:T, y=a_vec), Layout(;title="Target"))
	pπ = plot(scatter(;x=1:T, y=π_vec), Layout(;title="Inflation"))
	py = plot(scatter(;x=1:T, y=y_vec), Layout(;title="Output"))

	p = [pp pa; py pπ]
	relayout!(p, font_family = "Fira Sans Light", height = 600, width = 950, font_size = 12)

    return p
end
write(pwd()*"/../output.txt", "")


L_mat, ωmin, p1 = choose_ω()
p1

ct = CrazyType(; ω = ωmin)
Epfi!(ct);
p1, p2, p3 = makeplots_ct(ct);
p1


#=
ct = CrazyType()
Epfi!(ct, maxiter = 500)

p1, p2, p3 = makeplots_ct(ct);
p1
=#

# using JLD
# save("ct.jld", "ct", ct)

# plot_simul(ct, T=50)
