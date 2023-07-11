using Distributions

abstract type PhillipsCurve
end

abstract type Forward <: PhillipsCurve
end

abstract type Simultaneous <: PhillipsCurve
end

abstract type SemiForward <: PhillipsCurve
end

abstract type Fwd_strategy <: Forward
end

abstract type Fwd_literal <: Forward
end

abstract type Fwd_GP <: Forward
end

abstract type Plan{T<:PhillipsCurve}
end

mutable struct CrazyType{T<:PhillipsCurve} <: Plan{T}
	pars::Dict{Symbol, Float64}
	gr::Dict{Symbol, Vector{Float64}}

	opt::Dict{Symbol, Bool}
	
	gπ::Array{Float64, 2}
	ga::Array{Float64, 2}
	L::Array{Float64, 2}
	C::Array{Float64, 2}
	
	Ey::Array{Float64, 2}
	Eπ::Array{Float64, 2}
	Ep::Array{Float64, 2}
end

mutable struct DovisKirpalani{T<:PhillipsCurve} <: Plan{T}
	β::Float64
	γ::Float64
	κ::Float64
	σ::Float64
	ystar::Float64

	pgrid::Vector{Float64}
	agrid::Vector{Float64}

	Np::Int64
	Na::Int64

	gπ::Array{Float64, 2}
	ga::Array{Float64, 2}
	L::Array{Float64, 2}

	C::Array{Float64, 2}
	
	Ey::Array{Float64, 2}
	Eπ::Array{Float64, 2}
	Ep::Array{Float64, 2}
end
DovisKirpalani(ct::CrazyType{T}) where T <: PhillipsCurve = DovisKirpalani{T}(ct.β, ct.γ, ct.κ, ct.σ, ct.ystar, ct.pgrid, ct.agrid, ct.Np, ct.Na, ct.gπ, ct.ga, ct.L, ct.C, ct.Ey, ct.Eπ, ct.Ep)
CrazyType(dk::DovisKirpalani{T}; ω=0.0, χ=0.0, use_a = false, ψ = 0.0) where T <: PhillipsCurve = CrazyType{T}(dk.β, dk.γ, dk.κ, dk.σ, dk.ystar, ω, χ, use_a, ψ, dk.pgrid, dk.agrid, dk.Np, dk.Na, dk.gπ, dk.ga, dk.L, dk.C, dk.Ey, dk.Eπ, dk.Ep)

function switch_PC(pp::DovisKirpalani{T}, T2::PhillipsCurve) where T<:PhillipsCurve
	if T == T2
		return pp
	else
		return DovisKirpalani{T2}(pp.β, pp.γ, pp.κ, pp.σ, pp.ystar, pp.pgrid, pp.agrid, pp.Np, pp.Na, pp.gπ, pp.ga, pp.L, pp.C, pp.Ey, pp.Eπ, pp.Ep)
	end
end

mutable struct MultiType
	ct::CrazyType

	ωgrid::Vector{Float64}
	χgrid::Vector{Float64}

	z::Float64
	ν::Array{Float64, 3}

	μ::Array{Float64, 3}

	L_mat::Array{Float64, 4}
	C_mat::Array{Float64, 4}
	g_mat::Array{Float64, 4}
end

mutable struct Multiψ
	ct::CrazyType
	
	ωgrid::Vector{Float64}
	χgrid::Vector{Float64}
	ψgrid::Vector{Float64}

	L::Array{Float64, 5} # ω, χ, ψ, p, a
	C::Array{Float64, 5}
	g::Array{Float64, 5}
end

struct Prequel
    ωgrid::Vector{Float64}
    χgrid::Vector{Float64}
    agrid::Vector{Float64}
    pgrid::Vector{Float64}
    Agrid::Vector{Float64}
    L::Array{Float64, 5}    # ω, χ, a, p, A
    G::Array{Float64, 5}    # ω, χ, a, p, A
end

mutable struct Ramsey{T<:PhillipsCurve} <: Plan{T}
	β::Float64
	γ::Float64
	κ::Float64
	ystar::Float64

	Nθ::Int64
	θgrid::Vector{Float64}

	g::Array{Float64, 2}
	v::Array{Float64, 1}
end

mutable struct Sustainable{T<:PhillipsCurve} <: Plan{T}
	β::Float64
	γ::Float64
	κ::Float64
	ystar::Float64
	
	θ::Float64
	D::Float64
	σ::Float64

	ξ::Float64
	b::Float64

	Na::Int64
	agrid::Vector{Float64}

	g::Array{Float64, 2}
	v::Array{Float64, 1}
end

function Sustainable(ct::CrazyType{T}, Na = 100; ξ = Nash(ct), pc::DataType=Forward) where T <: PhillipsCurve
	β, γ, κ, ystar = ct.β, ct.γ, ct.κ, ct.ystar

	A = Nash(T, β, γ, κ, ystar)
	agrid = cdf.(Beta(2,2), range(0,1,length=Na))
	move_grids!(agrid, xmax=A, xmin=0.0)

	g = zeros(Na, 3)
	v = zeros(Na)

	yξ = (ystar - β*κ*γ*ξ) / (1+κ^2*γ)
	πξ = κ*yξ + β*ξ

	b = ((yξ - ystar)^2 + γ*πξ^2) / (1-β)
	θ = 0.25
	D = 0.005 * Nash(ct)
	σ = 0.01/4
	
	return Sustainable{pc}(β, γ, κ, ystar, θ, D, σ, ξ, b, Na, agrid, g, v)
end

function set_ξ!(sp::Sustainable, ξ)
	β, γ, κ, ystar = sp.β, sp.γ, sp.κ, sp.ystar
	yξ = (ystar - β*κ*γ*ξ) / (1+κ^2*γ)
	πξ = κ*yξ + β*ξ

	sp.ξ = ξ
	sp.b = ((yξ - ystar)^2 + γ*πξ^2) / (1-β)
	nothing
end

function Ramsey(ct::CrazyType{T}, Nθ=1000) where T<:PhillipsCurve
	β, γ, κ, ystar = ct.β, ct.γ, ct.κ, ct.ystar

	θgrid = range(-10, 10, length=Nθ)
	g = zeros(Nθ,3)

	g[:,3] .= (minimum(θgrid) + maximum(θgrid)) / 2

	v = zeros(Nθ)

	return Ramsey{T}(β, γ, κ, ystar, Nθ, θgrid, g, v)
end
function make_itp(rp::Ramsey, y)
	knots = (rp.θgrid,)
	itp = interpolate(knots, y, Gridded(Linear()))
	return itp
end
function make_itp(sp::Sustainable, y)
	knots = (sp.agrid,)
	itp = interpolate(knots, y, Gridded(Linear()))
	return itp
end

function move_grids!(xgrid; xmin=0.0, xmax=1.0)
	xgrid[:] = xgrid[:] * (xmax-xmin) .+ xmin
 	nothing
 end

function CrazyType(T::DataType;
		β = 1.02^(-0.25),
		# γ = 40.0,
		γ = 60.0,
		κ = 0.17,
		# κ = 0.8,
		# κ = 0.02,
		# σ = 0.01/3,
		# σ = 0.003,
		use_a = true,
		ψ = 0.01,
		σ = 0.01/4,
		ystar = 0.05,
		# ω = 0.271,
		# ω = 0.05,
		ω = 0.1,
		χ = 0.0,
		Np = 50,
		Na = 50,
		cut_0 = false,
		)

	if T == Simultaneous
		# γ = 1.75
	end

	pars = Dict(:β => β, :γ => γ, :κ => κ, :ψ => ψ, :σ => σ, :ystar => ystar, :ω => ω, :χ => χ)
	opt = Dict(:use_a => use_a)

	A = Nash(T, β, γ, κ, ystar)

	pgrid = cdf.(Beta(5,3), range(0,1,length=Np))
	if cut_0
		pgrid = cdf.(Beta(3,3), range(0,1,length=Np+1))[2:end]
	end
	agrid = cdf.(Beta(2,2), range(0,1,length=Na))
	move_grids!(agrid, xmax=A, xmin=0.0)

	gr = Dict(:p => pgrid, :a => agrid)

	gπ, ga = [zeros(Np, Na) for jj in 1:2]
	for jp in 1:Np, (ja, av) in enumerate(agrid)
		gπ[jp, ja] = av
		ga[jp, ja] = ϕ(av, ω, χ)
	end

	L = ones(Np, Na)
	C = ones(Np, Na)

	Ey = zeros(Np, Na)
	Eπ = zeros(Np, Na)
	Ep = zeros(Np, Na)

	return CrazyType{T}(pars, gr, opt, gπ, ga, L, C, Ey, Eπ, Ep)
end

N(ct, k::Symbol) = length(ct.gr[k])

function MultiType(ct::CrazyType;
	Nω = 40,
	Nχ = 30,
	χmin = 0.0,
	χmax = 0.6 * Nash(ct),
	)

	z = 1e-4
	
	Na = length(ct.gr[:a])
	Np = length(ct.gr[:p])

	ωgrid = -log.(range(0, 1, length=2+Nω))[2:end-1]
	χgrid = range(χmin, χmax, length = Nχ)
	
	ν = ones(Nω, Nχ, Na)
	μ = ones(Nω, Nχ, Na)

	L_mat = zeros(Nω, Nχ, Np, Na)
	C_mat = zeros(Nω, Nχ, Np, Na)
	g_mat = zeros(Nω, Nχ, Np, Na)

	return MultiType(ct, ωgrid, χgrid, z, ν, μ, L_mat, C_mat, g_mat)
end

function best_ct(mt::MultiType; jp = 2)

	_, jj = findmin(mt.L_mat[:,:,jp,:])

	jω = jj[1]
	jχ = jj[2]
	ja = jj[3]

	ωv = mt.ωgrid[jω]
	χv = mt.χgrid[jχ]

	ct = mt.ct

	ct.pars[:ω] = ωv
	ct.pars[:χ] = χv
	update_ga!(ct)
	ct.L .= mt.L_mat[jω, jχ, :, :]
	ct.gπ .= mt.g_mat[jω, jχ, :, :]

	return ct
end

function Multiψ(ct::CrazyType;
	Nω = 15,
	Nχ = 35,
	Nψ = 25,
	χmin = 0.0,
	χmax = 0.6 * Nash(ct),
	ψmax = 0.2
	)
	
	Na = length(ct.gr[:a])
	Np = length(ct.gr[:p])

	ωgrid = -log.(range(0, 1, length=2+Nω))[2:end-1]
	χgrid = range(χmin, χmax, length = Nχ)
	
	ψgrid = range(0, ψmax, length=Nψ)

	L = zeros(Nω, Nχ, Nψ, Np, Na)
	C = zeros(Nω, Nχ, Nψ, Np, Na)
	g = zeros(Nω, Nχ, Nψ, Np, Na)

	return Multiψ(ct, ωgrid, χgrid, ψgrid, L, C, g)
end

ϕ(a::Number, ω::Number, χ::Number) = exp(-ω) * (a-χ) + χ
ϕ(ct::Plan, a::Number) = ϕ(a, ct.pars[:ω], ct.pars[:χ])

function update_ga!(ct::CrazyType; ω = ct.pars[:ω], χ = ct.pars[:χ], ψ = ct.pars[:ψ])
	ct.pars[:ω] = ω
	ct.pars[:χ] = χ
	ct.pars[:ψ] = ψ
	for jp in eachindex(ct.gr[:p]), (ja, av) in enumerate(ct.gr[:a])
		ct.ga[jp, ja] = ϕ(ct, av)
	end
	nothing
end

which_PC(ct::Plan{T}) where T <: PhillipsCurve = T

Nash(T::DataType, β, γ, κ, ystar) = ifelse(T==Forward || T==SemiForward, κ / (1.0 - β + κ^2*γ) * ystar, ystar / (κ*γ))

# Nash(ct::CrazyType) = Nash(which_PC(ct), ct.β, ct.γ, ct.κ, ct.ystar)

Nash(ct::Plan{T}) where T <: PhillipsCurve = Nash(T, ct.pars[:β], ct.pars[:γ], ct.pars[:κ], ct.pars[:ystar])


dist_ϵ(ct) = Normal(0, ct.pars[:σ])
pdf_ϵ(ct, ϵv) = pdf.(dist_ϵ(ct), ϵv)
cdf_ϵ(ct, ϵv) = cdf.(dist_ϵ(ct), ϵv)

annualized(π::Real) = 100*((1.0 .+ π).^4 .- 1)
deannual(x::Real) = (x*0.01 + 1.0)^0.25 - 1.0

perc_rate(x) = 100 * (1 .- exp.(-x))


PC(ct::Plan{Forward}, obs_π, πe, exp_π′, πe′) = (1/ct.pars[:κ]) * (obs_π - ct.pars[:β] * exp_π′)
PC(ct::Plan{Simultaneous}, obs_π, πe, exp_π′, πe′) = 1/ct.pars[:κ]  * (obs_π - πe)
PC(ct::Plan{SemiForward}, obs_π, πe, exp_π′, πe′) = (1/ct.pars[:κ]) * (obs_π - ct.pars[:β] * πe′)

# PC(ct::Plan{Forward}, obs_π, pv, av, pprime, aprime, itp_gπ) = (1/ct.κ) * (obs_π - ct.β * (pprime * aprime + (1-pprime) * itp_gπ(pprime, aprime)))
# PC(ct::Plan{Simultaneous}, obs_π, pv, av, pprime, aprime, itp_gπ) = 1/ct.κ * (obs_π - (pv*av+(1-pv)*itp_gπ(pv,av)))
# PC(ct::Plan{SemiForward}, obs_π, pv, av, pprime, aprime, itp_gπ) = (1/ct.κ) * (obs_π - ct.β * (pv * aprime + (1-pv) * itp_gπ(pv, aprime)))
