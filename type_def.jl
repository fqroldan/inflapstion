using Distributions
mutable struct CrazyType
	β::Float64
	γ::Float64
	κ::Float64
	σ::Float64
	ystar::Float64
	ω::Float64
	χ::Float64

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

function move_grids!(xgrid; xmin=0.0, xmax=1.0)
	xgrid[:] = xgrid[:] * (xmax-xmin) .+ xmin
 	nothing
 end

function CrazyType(;
		β = 1.02^(-0.25),
		γ = 60.0,
		κ = 0.17,
		# κ = 0.8,
		# κ = 0.02,
		σ = 0.01/3,
		ystar = 0.05,
		# ω = 0.271,
		# ω = 0.05,
		ω = 0.1,
		χ = 0.0,
		Np = 20,
		Na = 30
		)

	A = κ / (1.0 - β + κ^2*γ) * ystar

	pgrid = cdf.(Beta(5,4), range(0,1,length=Np))
	agrid = cdf.(Beta(2,2), range(0,1,length=Na))
	move_grids!(agrid, xmax=A, xmin=0.0)

	gπ = ones(Np, Na) * A
	L = ones(Np, Na)

	Ey = zeros(Np, Na)
	Eπ = zeros(Np, Na)
	Ep = zeros(Np, Na)

	return CrazyType(β, γ, κ, σ, ystar, ω, χ, pgrid, agrid, Np, Na, gπ, L, Ey, Eπ, Ep)
end

ϕ(ct::CrazyType, a::Float64) = exp(-ct.ω) * (a-ct.χ) + ct.χ
Nash(ct::CrazyType) = ct.κ / (1.0 - ct.β + ct.κ^2*ct.γ) * ct.ystar

dist_ϵ(ct) = Normal(0, ct.σ)
pdf_ϵ(ct, ϵv) = pdf.(dist_ϵ(ct), ϵv)
cdf_ϵ(ct, ϵv) = cdf.(dist_ϵ(ct), ϵv)

annualized(π::Float64) = 100*((1.0 .+ π).^4 .- 1)
