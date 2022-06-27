using Distributed

using Distributions, Interpolations, Optim, HCubature, QuantEcon, Printf, PlotlyJS, Distributed, SharedArrays, Dates, JLD2

abstract type PhillipsCurve end

abstract type Forward <: PhillipsCurve end
abstract type Simultaneous <: PhillipsCurve end
abstract type SemiForward <: PhillipsCurve end
abstract type Plan{T<:PhillipsCurve} end

mutable struct CrazyType{T<:PhillipsCurve} <: Plan{T}
    β::Float64
    γ::Float64
    κ::Float64
    σ::Float64
    ystar::Float64
    ω::Float64
    χ::Float64

    use_a::Bool
    ψ::Float64

    pgrid::Vector{Float64}
    agrid::Vector{Float64}

    Np::Int64
    Na::Int64

    gπ::Array{Float64,2}
    ga::Array{Float64,2}
    L::Array{Float64,2}
    C::Array{Float64,2}

    Ey::Array{Float64,2}
    Eπ::Array{Float64,2}
    Ep::Array{Float64,2}
end

function move_grids!(xgrid; xmin=0.0, xmax=1.0)
    xgrid[:] = xgrid[:] * (xmax - xmin) .+ xmin
    nothing
end

function CrazyType(T::DataType;
    β=1.02^(-0.25),
    # γ = 40.0,
    γ=60.0,
    κ=0.17,
    # κ = 0.8,
    # κ = 0.02,
    # σ = 0.01/3,
    # σ = 0.003,
    use_a=true,
    ψ=0.01,
    σ=0.01 / 4,
    ystar=0.05,
    # ω = 0.271,
    # ω = 0.05,
    ω=0.1,
    χ=0.0,
    Np=50,
    Na=50
)

    if T == Simultaneous
        # γ = 1.75
    end

    A = Nash(T, β, γ, κ, ystar)

    pgrid = cdf.(Beta(5, 3), range(0, 1, length=Np))
    agrid = cdf.(Beta(2, 2), range(0, 1, length=Na))
    move_grids!(agrid, xmax=A, xmin=0.0)

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

    return CrazyType{T}(β, γ, κ, σ, ystar, ω, χ, use_a, ψ, pgrid, agrid, Np, Na, gπ, ga, L, C, Ey, Eπ, Ep)
end

Nash(T::DataType, β, γ, κ, ystar) = ifelse(T == Forward || T == SemiForward, κ / (1.0 - β + κ^2 * γ) * ystar, ystar / (κ * γ))
Nash(ct::Plan{T}) where {T<:PhillipsCurve} = Nash(T, ct.β, ct.γ, ct.κ, ct.ystar)

ϕ(a::Float64, ω::Float64, χ::Float64) = exp(-ω) * (a - χ) + χ
ϕ(ct::Plan, a::Float64) = ϕ(a, ct.ω, ct.χ)

function update_ga!(ct::CrazyType; ω=ct.ω, χ=ct.χ)
    ct.ω = ω
    ct.χ = χ
    for jp in 1:ct.Np, (ja, av) in enumerate(ct.agrid)
        ct.ga[jp, ja] = ϕ(ct, av)
    end
    nothing
end

dist_ϵ(ct) = Normal(0, ct.σ)
pdf_ϵ(ct, ϵv) = pdf.(dist_ϵ(ct), ϵv)
cdf_ϵ(ct, ϵv) = cdf.(dist_ϵ(ct), ϵv)

annualized(π::Real) = 100 * ((1.0 .+ π) .^ 4 .- 1)
deannual(x::Real) = (x * 0.01 + 1.0)^0.25 - 1.0

perc_rate(x) = 100 * (1 .- exp.(-x))

PC(ct::Plan{Forward}, obs_π, πe, exp_π′, πe′) = (1 / ct.κ) * (obs_π - ct.β * exp_π′)
PC(ct::Plan{Simultaneous}, obs_π, πe, exp_π′, πe′) = 1 / ct.κ * (obs_π - πe)
PC(ct::Plan{SemiForward}, obs_π, πe, exp_π′, πe′) = (1 / ct.κ) * (obs_π - ct.β * πe′)

function Bayes(ct::Plan, obs_π, exp_π, pv, av)

    if isapprox(pv, 0.0)
        p′ = 0.0
    elseif isapprox(pv, 1.0)
        p′ = 1.0
    else
        numer = pv * pdf_ϵ(ct, obs_π - av)
        denomin = numer + (1.0 - pv) * pdf_ϵ(ct, obs_π - exp_π)
        if isapprox(denomin, 0.0)
            p′ = 0.0
        else
            p′ = numer / denomin
        end
    end

    return p′
end

next_a(ct::Plan, av, π) = ifelse(ct.use_a, ϕ(ct, av), ϕ(ct, av) + ct.ψ * (π - av))

function cond_L(ct::Plan{T}, itp_gπ, itp_L, itp_C, obs_π, pv, av, ge, πe′) where {T<:PhillipsCurve}
    # ge = itp_gπ(pv, av)
    pprime = Bayes(ct, obs_π, ge, pv, av)

    πe = pv * av + (1 - pv) * ge

    a_min, a_max = extrema(ct.agrid)

    π_today = max(a_min, min(a_max, obs_π))
    aprime = next_a(ct, av, π_today)
    aprime = max(a_min, min(a_max, aprime))

    # if aprime <= minimum(ct.agrid) || aprime >= maximum(ct.agrid)
    itp_L = extrapolate(itp_L, Interpolations.Flat())
    itp_gπ = extrapolate(itp_gπ, Interpolations.Flat())
    itp_C = extrapolate(itp_C, Interpolations.Flat())
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

get_sumprob(ct::Plan) = cdf_ϵ(ct, 3.09 * ct.σ) - cdf_ϵ(ct, -3.09 * ct.σ)

function exp_L_y(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, ge, πe′)

    sum_prob = get_sumprob(ct)

    f_p(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, ge, πe′)[2] * pdf_ϵ(ct, ϵv)
    Ep, err = hquadrature(f_p, -3.09 * ct.σ, 3.09 * ct.σ, rtol=1e-10, atol=0, maxevals=0)
    f_y(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, ge, πe′)[3] * pdf_ϵ(ct, ϵv)
    Ey, err = hquadrature(f_y, -3.09 * ct.σ, 3.09 * ct.σ, rtol=1e-10, atol=0, maxevals=0)
    f_C(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, ge, πe′)[4] * pdf_ϵ(ct, ϵv)
    Ec, err = hquadrature(f_C, -3.09 * ct.σ, 3.09 * ct.σ, rtol=1e-10, atol=0, maxevals=0)

    Ey = Ey / sum_prob
    Ep = Ep / sum_prob
    Ec = Ec / sum_prob

    return Ey, Ep, Ec
end

function exp_L(ct::Plan, itp_gπ, itp_L, itp_C, control_π, pv, av, ge, πe′)

    f(ϵv) = cond_L(ct, itp_gπ, itp_L, itp_C, control_π + ϵv, pv, av, ge, πe′)[1] * pdf_ϵ(ct, ϵv)
    val, err = hquadrature(f, -3.09 * ct.σ, 3.09 * ct.σ, rtol=1e-10, atol=0, maxevals=0)
    sum_prob = get_sumprob(ct)

    return val / sum_prob
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′)
    π_guess = xguess[1]
    minπ = max(0, π_guess - 3.09 * ct.σ)
    maxπ = min(1.1 * maximum(ct.agrid), π_guess + 3.09 * ct.σ)
    if maxπ < minπ + 1.1 * maximum(ct.agrid) / 10
        maxπ = minπ + 1.1 * maximum(ct.agrid) / 10
    end

    # aprime = ϕ(ct, av)
    aprime = xguess[2]

    gπ, L = π_guess, itp_L(pv, av)

    obj_f(x) = exp_L(ct, itp_gπ, itp_L, itp_C, first(x), pv, av, ge, πe′)

    try
        res = Optim.optimize(
            gπ -> obj_f(first(gπ)),
            [π_guess], LBFGS()
        )
        # od = OnceDifferentiable(obj_f, [π_guess]; autodiff = :forward)
        # # res = Optim.optimize(od, [π_guess], BFGS())
        # res = Optim.optimize(od, [minπ], [maxπ], [π_guess], Fminbox(BFGS()))

        gπ::Float64, L::Float64 = first(res.minimizer), res.minimum
        if Optim.converged(res)
            return gπ, L, aprime
        end
    catch
    end

    resb = Optim.optimize(obj_f, minπ, maxπ, Brent())
    gπ, L = resb.minimizer, resb.minimum

    return gπ, L, aprime
end

function π_prime(ct, ϵv, pv, av, itp_gπ, ge, aprime)
    exp_π = pv * av + (1 - pv) * ge
    pprime = Bayes(ct, ge + ϵv, exp_π, pv, av)

    ge′ = itp_gπ(pprime, aprime)

    return pv * aprime + (1 - pv) * ge′
end

function exp_π_prime(ct::Plan{SemiForward}, pv, av, itp_gπ, ge, aprime)
    f(ϵv) = π_prime(ct, ϵv, pv, av, itp_gπ, ge, aprime) * pdf_ϵ(ct, ϵv)
    val::Float64, err::Float64 = hquadrature(f, -3.09 * ct.σ, 3.09 * ct.σ, rtol=1e-10, atol=0, maxevals=0)
    sum_prob = get_sumprob(ct)

    return val / sum_prob
end

function exp_π_prime(ct::Plan{T}, pv, av, itp_gπ, ge, aprime) where {T<:PhillipsCurve}
	return 0.0
end

function optim_step(ct::Plan, itp_gπ, itp_L, itp_C, gπ_guess; optimize::Bool=true)
    gπ, ga = zeros(size(ct.gπ)), zeros(size(ct.ga))
    L = zeros(size(ct.L))
    Ey, Eπ = zeros(size(ct.Ey)), zeros(size(ct.Eπ))
    Ep, C = zeros(size(ct.Ep)), zeros(size(ct.C))
    πN = Nash(ct)
    maxa = maximum(ct.agrid)
    mina = minimum(ct.agrid)
    length_a = maxa - mina
    apgrid = gridmake(1:ct.Np, 1:ct.Na)
    Threads.@threads for js in 1:size(apgrid, 1)
        # for js in 1:size(apgrid,1)
        jp, ja = apgrid[js, :]
        pv, av = ct.pgrid[jp], ct.agrid[ja]

        a_guess = max(min(ct.ga[jp, ja], maxa), mina)
        π_guess = max(min(gπ_guess[jp, ja], maxa), mina)
        xguess = [π_guess, a_guess]
        ge = itp_gπ(pv, av)
        πe′ = exp_π_prime(ct, pv, av, itp_gπ, ge, a_guess)
        if optimize
            # π_guess = itp_gπ(pv, av)
            gπ[jp, ja], L[jp, ja], aprime = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′)
        else
            aprime = a_guess
            gπ[jp, ja] = π_guess
            L[jp, ja] = exp_L(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av, ge, πe′)
        end
        ga[jp, ja] = aprime
        Ey[jp, ja], Ep[jp, ja], EC′ = exp_L_y(ct, itp_gπ, itp_L, itp_C, π_guess, pv, av, ge, πe′)
        Eπ[jp, ja] = pv * av + (1.0 - pv) * gπ[jp, ja]

        if av >= πN || isapprox(av, πN)
            C[jp, ja] = (1 - ct.β) * 1 + ct.β * EC′
        else
            C[jp, ja] = (1 - ct.β) * (πN - Eπ[jp, ja]) / (πN - av) + ct.β * EC′
        end
    end

    return gπ, L, Ey, Eπ, Ep, C, ga
end

function pf_iter(ct::Plan, Egπ, gπ_guess; optimize::Bool=true)
    knots = (ct.pgrid, ct.agrid)
    itp_gπ = interpolate(knots, Egπ, Gridded(Linear()))
    itp_L = interpolate(knots, ct.L, Gridded(Linear()))
    itp_C = interpolate(knots, ct.C, Gridded(Linear()))

    new_gπ, new_L, new_y, new_π, new_p, new_C, new_a = optim_step(ct, itp_gπ, itp_L, itp_C, gπ_guess; optimize=optimize)

    return new_gπ, new_L, [new_y, new_π, new_p, new_C, new_a]
end

function update_others!(ct::Plan, new_others, upd_η2)
    new_y, new_π, new_p, new_C, new_a = new_others[:]
    ct.Ey = new_y
    ct.Eπ = new_π
    ct.Ep = new_p
    ct.C = new_C
    ct.ga = ct.ga + upd_η2 * (new_a - ct.ga)
    nothing
end

function pfi!(ct::Plan, Egπ; tol::Float64=1e-12, maxiter::Int64=250, verbose::Bool=true, reset_guess::Bool=false)
    dist = 10.0
    iter = 0
    upd_η2 = 0.75

    rep = "\nStarting PFI (tol = $(@sprintf("%0.3g",tol)))"
    verbose ? print(rep) : print(rep)

    if reset_guess
        # ct.gπ = zeros(size(ct.gπ))
        ct.L = ones(ct.Np, ct.Na)
    end

    old_gπ = copy(Egπ)
    old_L = copy(ct.L)
    new_gπ = zeros(size(old_gπ))

    while dist > tol && iter < maxiter
        iter += 1

        for jj in 1:10
            _, new_L, _ = pf_iter(ct, Egπ, old_gπ; optimize=false)
            ct.L = upd_η2 * new_L + (1.0 - upd_η2) * ct.L
        end
        old_L .= ct.L

        new_gπ, new_L, new_others = pf_iter(ct, Egπ, old_gπ)
        update_others!(ct, new_others, upd_η2)

        norm_L = max(sqrt.(sum(old_L .^ 2)) / length(old_L), 100tol)
        dist = sqrt.(sum((new_L - old_L) .^ 2)) / length(old_L) / norm_L

        ct.L = upd_η2 * new_L + (1.0 - upd_η2) * ct.L
        old_gπ = upd_η2 * new_gπ + (1.0 - upd_η2) * old_gπ
    end

    if verbose && dist <= tol
        print("\nConverged in $iter iterations.")
    elseif verbose
        print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))")
    end

    return (dist <= tol), new_gπ
end

decay_η(ct::Plan, η) = max(0.95 * η, 1e-6)

function report_start(ct::CrazyType)
    print("\nRun with ω = $(@sprintf("%.3g",ct.ω)), χ = $(@sprintf("%.3g",annualized(ct.χ)))% at $(Dates.format(now(), "HH:MM"))\n")
    nothing
end

function Epfi!(ct::Plan; tol::Float64=5e-4, maxiter::Int64=2500, verbose::Bool=true, tempplots::Bool=false, upd_η::Float64=0.01, switch_η=10, tol_pfi=2e-3 / 0.99)
    dist = 10.0
    iter = 0

    report_start(ct)
    dists = []
    old_gπ, old_L, old_ga = similar(ct.gπ), similar(ct.L), similar(ct.ga)

    reset_guess = false
    while dist > tol && iter < maxiter
        iter += 1
        tol_pfi = max(tol_pfi * 0.98, 2e-6)

        old_gπ .= ct.gπ
        old_L .= ct.L
        old_ga .= ct.ga

        # reset_L!(ct)

        flag, new_gπ = pfi!(ct, old_gπ; verbose=verbose, reset_guess=reset_guess, tol=tol_pfi)
        reset_guess = !flag

        norm_gπ = sqrt.(sum(annualized.(ct.gπ) .^ 2)) / length(annualized.(ct.gπ))
        dist_π = sqrt.(sum((annualized.(new_gπ) - annualized.(ct.gπ)) .^ 2)) / length(annualized.(ct.gπ)) / max(norm_gπ, 20tol)
        norm_ga = sqrt.(sum(annualized.(old_ga) .^ 2)) / length(annualized.(old_ga))
        dist_a = sqrt.(sum((annualized.(ct.ga) - annualized.(old_ga)) .^ 2)) / length(annualized.(old_ga)) / max(norm_ga, 20tol)
        dist = max(dist_π, dist_a / 10)

        push!(dists, dist)
        rep_status = "\nAfter $iter iterations, d(π) = $(@sprintf("%0.3g",dist)) at |π,a| = ($(@sprintf("%0.3g",norm_gπ)), $(@sprintf("%0.3g",norm_ga)))"
        if flag
            rep_status *= "✓ "
        end

        verbose && print(rep_status * "\n")

        ct.gπ = upd_η * new_gπ + (1 - upd_η) * ct.gπ
        ct.ga = upd_η * ct.ga + (1 - upd_η) * old_ga

        if iter == floor(Int, switch_η * 0.4)
            upd_η = min(upd_η, 0.01)
        elseif iter % switch_η == 0
            upd_η = decay_η(ct, upd_η) # Automatically uses the updating method for fwd or bwd
        end
        if verbose
            print("new upd_η = $(@sprintf("%0.3g", upd_η))")
        end

    end

    # Update credibility
    tolC, maxiterC = 1e-5, 2000
    dist2, iter = 1 + tolC, 0
    while dist2 > tolC && iter < maxiterC
        iter += 1

        old_C = copy(ct.C)
        Egπ = copy(ct.gπ)
        _, _, new_others = pf_iter(ct, Egπ, Egπ, optimize=false)

        new_C = new_others[4]
        dist2 = sqrt.(sum((new_C - old_C) .^ 2)) / sqrt.(sum(old_C .^ 2))

        update_others!(ct, new_others, 0.5)
    end

    if verbose && dist <= tol
        print("\nConverged in $iter iterations.", true)
    elseif verbose
        print("\nAfter $iter iterations, d(L) = $(@sprintf("%0.3g",dist))", true)
    end

    return dist
end

