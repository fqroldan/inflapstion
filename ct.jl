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
end
function CrazyType(;
		β = 0.96,
		γ = 1.0,
		α = 0.17,
		σ = 0.3,
		ystar = 0.05,
		#ω = 0.271,
		ω = 0.15,
		Np = 45,
		Na = 35
		)

	y = 1.0 / (1.0+2.0*γ*α^2) * ystar - 2.0*γ*α / (1.0+2.0*γ*α^2) * 0.0

	A = α / (1.0-β*exp(-ω)) * y

	curv = 0.25
	pgrid = range(0, 1, length=Np).^(1.0/curv)
	agrid = range(0, 1.0*A, length=Na)

	gπ = zeros(Np, Na)
	L = zeros(Np, Na)
	for jp in 1:Np
		for (ja, av) in enumerate(agrid)
			y_Nash = 1.0 / (1.0+2.0*γ*α^2) * ystar - 2.0*γ*α / (1.0+2.0*γ*α^2) * av
			gπ[jp, ja] = α * y_Nash + exp(-ω) * av
		end
	end

	return CrazyType(β, γ, α, σ, ystar, ω, pgrid, agrid, Np, Na, gπ, L)
end

ϕ(ct::CrazyType, a::Float64) = exp(-ct.ω) * a

dist_ϵ(ct) = Normal(0, ct.σ)
pdf_ϵ(ct, ϵv) = pdf.(dist_ϵ(ct), ϵv)
function Bayes(ct::CrazyType, obs_π, exp_π, av, pv)

	numer = pv * pdf_ϵ(ct, obs_π - av)
	denomin = numer + (1.0-pv) * pdf_ϵ(ct, obs_π - exp_π)

	return numer / denomin
end

NKPC(ct::CrazyType, obs_π, exp_π′) = (1.0/ct.α) * (obs_π - ct.β * exp_π′)

function cond_L(ct::CrazyType, itp_gπ, itp_L, obs_π, av, pv)
	exp_π  = itp_gπ(pv, av)
	if pv == ct.pgrid[1]
		pprime = 0.0
	elseif pv == ct.pgrid[end]
		pprime = 1.0
	else
		pprime = Bayes(ct, obs_π, exp_π, av, pv)
	end
	aprime = ϕ(ct, av)
	gπ′ = itp_gπ(pprime, aprime)
	exp_π′ = pprime * aprime + (1.0-pprime) * gπ′

	y = NKPC(ct, obs_π, exp_π′)

	L′ = itp_L(pprime, aprime)

	L = (y-ct.ystar)^2 + ct.γ * obs_π^2 + ct.β * L′

	return L
end

function exp_L(ct::CrazyType, itp_gπ, itp_L, control_π, av, pv)

	f(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, av, pv) * pdf_ϵ(ct, ϵv)
	(val, err) = hquadrature(f, -3.09*ct.σ, 3.09*ct.σ, rtol=1e-32, atol=0, maxevals=0)

	sum_prob, err = hquadrature(x -> pdf_ϵ(ct, x), -3.09*ct.σ, 3.09*ct.σ, rtol=1e-32, atol=0, maxevals=0)

	val = val / sum_prob

	return val
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, av, pv)

	minπ, maxπ = -0.1, 1.1*maximum(ct.agrid)
	res = Optim.optimize(
			gπ -> exp_L(ct, itp_gπ, itp_L, gπ, av, pv),
			minπ, maxπ, Brent()#, reltol=1e-32, abstol=1e-32
			)
	gπ = res.minimizer
	L = res.minimum

	return gπ, L
end

function optim_step(ct::CrazyType, itp_gπ, itp_L; optimize::Bool=true)
	gπ, L = SharedArray{Float64}(ct.gπ), SharedArray{Float64}(ct.L)
	# gπ, L = Array{Float64}(undef, size(ct.gπ)), Array{Float64}(undef, size(ct.L))
	apgrid = gridmake(1:ct.Np, 1:ct.Na)
	@sync @distributed for js in 1:size(apgrid,1)
    # for js in 1:size(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.pgrid[jp], ct.agrid[ja]
		if optimize
			gπ[jp, ja], L[jp, ja] = opt_L(ct, itp_gπ, itp_L, av, pv)
		else
			gπ[jp, ja] = ct.gπ[jp, ja]
			L[jp, ja] = exp_L(ct, itp_gπ, itp_L, gπ[jp, ja], av, pv)
		end
	end

	return gπ, L
end

function pf_iter(ct::CrazyType; optimize::Bool=true)
	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()))
	itp_L  = interpolate(knots, ct.L , Gridded(Linear()))

	new_gπ, new_L = optim_step(ct, itp_gπ, itp_L; optimize=optimize)

	return new_gπ, new_L
end

function pfi!(ct::CrazyType; tol::Float64=1e-6, maxiter::Int64=2500, verbose::Bool=true)
	dist = 10.
	iter = 0
	upd_η = 0.75
    if verbose
        print_save("\nStarting PFI")
    end

	while dist > tol && iter < maxiter
		iter += 1
		old_gπ, old_L = copy(ct.gπ), copy(ct.L)

		new_gπ, new_L = pf_iter(ct)

		dist_π = sqrt.(sum( (new_gπ - old_gπ).^2 )) / sqrt.(sum(old_gπ.^2))
		dist_L = sqrt.(sum( (new_L  - old_L ).^2 )) / sqrt.(sum(old_L .^2))

		dist = max(dist_π, dist_L)

		for jj in 1:5
			new_gπ, new_L = pf_iter(ct; optimize=false)
		end

		ct.gπ = upd_η*0.2 * new_gπ + (1.0-upd_η*0.2) * old_gπ
		ct.L  = upd_η * new_L  + (1.0-upd_η) * old_L

		if verbose && iter % 10 == 0
			print_save("\nAfter $iter iterations, d(π, L) = ($(@sprintf("%0.3g",dist_π)), $(@sprintf("%0.3g",dist_L)))")
		end
	end
	return (dist <= tol)
end

function plot_ct(ct::CrazyType; make_pdf::Bool=false, make_png::Bool=false)
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
	pπa = lines(ct, ct.gπ, dim = 1, title="gπ", showleg = true)
	pπp = lines(ct, ct.gπ, dim = 2, title="gπ", showleg = true)
	pLa = lines(ct, ct.L , dim = 1, title="𝓛")
	pLp = lines(ct, ct.L , dim = 2, title="𝓛")

	p = [pπa pπp; pLa pLp]
	relayout!(p, font_family = "Fira Sans Light", font_size = 12, height = 600, width = 950)

	function makeplot(p, ext)
		savefig(p, pwd() * "/../Graphs/ct" * ext)
	end

	if make_pdf
		makeplot(p, ".pdf")
	end
	if make_png
		makeplot(p, ".png")
	end

	return p
end

end # everywhere

function choose_ω()
	Nω = 25
	ωgrid = range(0.0, 0.75, length=Nω)

	ct = CrazyType()
	L_mat = zeros(Nω, ct.Np, ct.Na)

	L_min = 100.
	ωmin = 1.0

	jamin_vec = zeros(Nω)
	for (jω, ωv) in enumerate(ωgrid)
		Lguess, πguess = ct.L, ct.gπ
		ct = CrazyType(; ω = ωv)
		ct.L = Lguess
		ct.gπ = πguess

		flag = pfi!(ct, verbose = false)

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
		end
	end

	p1 = plot([
		scatter(;x=ωgrid, y=L_min[:, 2, jamin_vec[:]])
		])

	return L_mat, ωmin, p1
end

function iter_simul(ct::CrazyType, itp_gπ, pv, av)
	ϵ = rand(dist_ϵ(ct))

	exp_π = itp_gπ(pv, av)
	obs_π = exp_π+ϵ
	pprime = Bayes(ct, obs_π, exp_π, av, pv)
	aprime = ϕ(ct, av)
	gπ′ = itp_gπ(pprime, aprime)
	exp_π′ = pprime * aprime + (1.0-pprime) * gπ′

	y = NKPC(ct, obs_π, exp_π′)

	return pprime, aprime, obs_π, y
end

function simul(ct::CrazyType; T::Int64=50)
	p0 = ct.pgrid[2]

	_, ind_a0 = findmin(ct.L[2, :])
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

function plot_simul(ct::CrazyType; T::Int64=50)
	p_vec, a_vec, π_vec, y_vec = simul(ct, T=T)

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
# ct = CrazyType(; ω = ωmin)

# pfi!(ct)
# plot_ct(ct)

# using JLD
# save("ct.jld", "ct", ct)

# plot_simul(ct, T=50)
