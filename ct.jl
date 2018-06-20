@everywhere begin

using Distributions, Interpolations, Optim, Cubature, PlotlyJS, Rsvg, QuantEcon, LaTeXStrings

type CrazyType
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
        β = 0.9,
        γ = 1.,
        α = 0.17,
        σ = 0.125,
        ystar = 0.05,
        ω = 0.3,
        Np = 25,
        Na = 25
        )

    A = 1/(α*γ) * ystar

    curv = 0.6
    pgrid = linspace(0, 1, Np).^(1./curv)
    agrid = linspace(0, 1.25*A, Na)

    gπ = zeros(Np, Na)
    L = zeros(Np, Na)
    for jp in 1:Np
        for (ja, av) in enumerate(agrid)
            gπ[jp, ja] = av
        end
    end

    return CrazyType(β, γ, α, σ, ystar, ω, pgrid, agrid, Np, Na, gπ, L)
end

ϕ(ct::CrazyType, a::Float64) = exp(-ct.ω) * a

dist_ϵ(ct, ϵv) = pdf.(Normal(0.,ct.σ), ϵv)
function Bayes(ct::CrazyType, obs_π, exp_π, av, pv)

    numerator = pv * dist_ϵ(ct, obs_π - av)
    denominator = numerator + (1.-pv) * dist_ϵ(ct, obs_π - exp_π)

    return numerator / denominator
end

function cond_L(ct::CrazyType, itp_gπ, itp_L, obs_π, av, pv)
    exp_π  = itp_gπ[pv, av]
    pprime = Bayes(ct, obs_π, exp_π, av, pv)
    aprime = ϕ(ct, av)
    gπ′ = itp_gπ[pprime, aprime]
    exp_π′ = pprime * aprime + (1.-pprime) * gπ′

    y = (1./ct.α) * (obs_π - exp_π′)

    L′ = itp_L[pprime, aprime]

    L = (y-ct.ystar)^2 + ct.γ * obs_π^2 + ct.β * L′

    return L
end

function exp_L(ct::CrazyType, itp_gπ, itp_L, control_π, av, pv)

    f(ϵv) = cond_L(ct, itp_gπ, itp_L, control_π + ϵv, av, pv) * dist_ϵ(ct, ϵv)

    (val, err) = hquadrature(f, -1.96*ct.σ, 1.96*ct.σ, reltol=1e-6, abstol=0, maxevals=0)

    return val
end

function opt_L(ct::CrazyType, itp_gπ, itp_L, av, pv)

    minπ, maxπ = -0.9, 1.25 / (ct.α*ct.γ) * ct.ystar
    res = Optim.optimize(
            gπ -> exp_L(ct, itp_gπ, itp_L, gπ, av, pv),
            minπ, maxπ, GoldenSection()
            )
    gπ = res.minimizer
    L = res.minimum

    return gπ, L
end

function optim_step(ct::CrazyType, itp_gπ, itp_L)

    gπ, L = SharedArray(ct.gπ), SharedArray(ct.L)
    apgrid = gridmake(1:ct.Np, 1:ct.Na)
    @sync @parallel for js in 1:size(apgrid,1)
        jp, ja = apgrid[js, :]
        pv, av = ct.pgrid[jp], ct.agrid[ja]
        gπ[jp, ja], L[jp, ja] = opt_L(ct, itp_gπ, itp_L, av, pv)
    end

    return gπ, L
end

function pf_iter(ct::CrazyType)
    knots = (ct.pgrid, ct.agrid)
    itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()))
    itp_L  = interpolate(knots, ct.L , Gridded(Linear()))

    new_gπ, new_L = optim_step(ct, itp_gπ, itp_L)

    return new_gπ, new_L
end

function pfi!(ct::CrazyType; tol::Float64=1e-6, maxiter::Int64=1000, verbose::Bool=true)
    dist = 10.
    iter = 0
    upd_η = 0.33

    while dist > tol && iter < maxiter
        iter += 1
        old_gπ, old_L = copy(ct.gπ), copy(ct.L)

        new_gπ, new_L = pf_iter(ct)

        dist_π = sqrt.(sum( (new_gπ - old_gπ).^2 )) / sqrt.(sum(old_gπ.^2))
        dist_L = sqrt.(sum( (new_L  - old_L ).^2 )) / sqrt.(sum(old_L .^2))

        dist = max(dist_π, dist_L)

        ct.gπ = upd_η * new_gπ + (1.-upd_η) * old_gπ
        ct.L  = upd_η * new_L  + (1.-upd_η) * old_L

        if verbose && iter % 10 == 0
            print("\nAfter $iter iterations, d(π, L) = ($(@sprintf("%0.3g",dist_π)), $(@sprintf("%0.3g",dist_L)))")
        end
    end
end

function plot_ct(ct::CrazyType)
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
            xtitle= "p"
        elseif dim == 2
            xgrid = ct.agrid
            zgrid = ct.pgrid
            xtitle= "a"
        else
            throw(error("wrong dim"))
        end
        Nz = length(zgrid)
        l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(Nz)
        for (jz, zv) in enumerate(zgrid)
            if dim == 1
                y_vec = y_mat[:, jz]
                name = "a"
            elseif dim == 2
                y_vec = y_mat[jz, :]
                name = "p"
            end
            name = name * " = $(round(zv,2))"
            jz % 2 == 0? showleg_i = showleg: showleg_i = false
            l_new = scatter(;x=xgrid, y=y_vec, name = name, showlegend = showleg_i, marker_color=col[ceil(Int,10*jz/Nz)])
            l[jz] = l_new
        end
        p = plot([l[jz] for jz in 1:Nz], Layout(;title=title, xaxis_title=xtitle))
        return p
    end
    pπa = lines(ct, ct.gπ, dim = 1, title="gπ", showleg = true)
    pπp = lines(ct, ct.gπ, dim = 2, title="gπ", showleg = true)
    pLa = lines(ct, ct.L , dim = 1, title="L")
    pLp = lines(ct, ct.L , dim = 2, title="L")

    p = [pπa pπp; pLa pLp]
    p.plot.layout["font_family"] = "Fira Sans Light"
    p.plot.layout["height"] = 600
    p.plot.layout["width"]  = 950
    p.plot.layout["font_size"] = 12

    savefig(p, pwd() * "/../Graphs/ct.png")
    Void
end

end # everywhere
ct = CrazyType()

pfi!(ct)
plot_ct(ct)
