using Distributions

abstract type PhillipsCurve
end

abstract type Forward <: PhillipsCurve
end

abstract type Backward <: PhillipsCurve
end

mutable struct CrazyType{T<:PhillipsCurve}
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
	C::Array{Float64, 2}
	
	Ey::Array{Float64, 2}
	Eπ::Array{Float64, 2}
	Ep::Array{Float64, 2}
end

function move_grids!(xgrid; xmin=0.0, xmax=1.0)
	xgrid[:] = xgrid[:] * (xmax-xmin) .+ xmin
 	nothing
 end

function CrazyType(T::DataType;
		β = 1.02^(-0.25),
		γ = 40.0,
		# γ = 60.0,
		κ = 0.17,
		# κ = 0.8,
		# κ = 0.02,
		# σ = 0.01/3,
		# σ = 0.003,
		σ = 0.01/4,
		ystar = 0.05,
		# ω = 0.271,
		# ω = 0.05,
		ω = 0.1,
		χ = 0.0,
		Np = 30,
		Na = 30
		)

	if T == Backward
		# γ = 1.75
	end

	A = Nash(T, β, γ, κ, ystar)

	pgrid = cdf.(Beta(5,2), range(0,1,length=Np))
	agrid = cdf.(Beta(2,2), range(0,1,length=Na))
	move_grids!(agrid, xmax=A, xmin=0.0)

	gπ = ones(Np, Na) * A
	L = ones(Np, Na)
	C = ones(Np, Na)

	Ey = zeros(Np, Na)
	Eπ = zeros(Np, Na)
	Ep = zeros(Np, Na)

	return CrazyType{T}(β, γ, κ, σ, ystar, ω, χ, pgrid, agrid, Np, Na, gπ, L, C, Ey, Eπ, Ep)
end

ϕ(ct::CrazyType, a::Float64) = exp(-ct.ω) * (a-ct.χ) + ct.χ

which_PC(ct::CrazyType{T}) where T <: PhillipsCurve = T

Nash(T::DataType, β, γ, κ, ystar) = ifelse(T==Forward, κ / (1.0 - β + κ^2*γ) * ystar, ystar / (κ*γ))

Nash(ct::CrazyType) = Nash(which_PC(ct), ct.β, ct.γ, ct.κ, ct.ystar)

dist_ϵ(ct) = Normal(0, ct.σ)
pdf_ϵ(ct, ϵv) = pdf.(dist_ϵ(ct), ϵv)
cdf_ϵ(ct, ϵv) = cdf.(dist_ϵ(ct), ϵv)

annualized(π::Real) = 100*((1.0 .+ π).^4 .- 1)
deannual(x::Real) = (x*0.01 + 1.0)^0.25 - 1.0