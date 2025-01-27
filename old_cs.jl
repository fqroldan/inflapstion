function plot_cs(k::Symbol; T = 11, jp = 2, showLmat = false, showL = false, slides = true, dark = false, share = false)

    s = sk(k)

    K = ifelse(k == :β, 19, 20)

    Xdict = Dict(k => zeros(K)*NaN for k in (:β, :σ, :κ))
    ωvec = zeros(K) * NaN
    χvec = zeros(K) * NaN
    avec = zeros(K) * NaN

    Lvec = zeros(K)
    Lmat = Vector{Vector{Float64}}()
    pgrid = Vector{Float64}()
    πN = 0.

    plans = zeros(T, K)
    for jj in 1:K
        if k == :κ
            name = "mt"
        else
            name = "mt$jj"
        end
        mt = load("Output/CompStats/$s/run$jj/$name.jld2", "mt")
                    
        Xdict[:σ][jj] = mt.ct.σ
        Xdict[:β][jj] = mt.ct.β
        Xdict[:κ][jj] = mt.ct.κ

        L, jmin = findmin(mt.L_mat[:, :, jp, :])

        lm = [findmin(mt.L_mat[:, :, jpx, :])[1] for jpx in 2:length(mt.ct.pgrid)-1]
        push!(Lmat, lm)
        pgrid = mt.ct.pgrid[2:end-1]

        jω = jmin[1]
        jχ = jmin[2]
        ja = jmin[3]

        Ω = exp(-mt.ωgrid[jω])
        χ = annualized(mt.χgrid[jχ])
        a = annualized(mt.ct.agrid[ja])
        ωvec[jj] = mt.ωgrid[jω]
        χvec[jj] = mt.χgrid[jχ]
        avec[jj] = mt.ct.agrid[ja]

        πN = Nash(Forward, mt.ct.β, mt.ct.γ, mt.ct.κ, mt.ct.ystar)
        
        for tt in 1:T
            plans[tt, jj] = χ + Ω^(tt - 1) * (a - χ)
        end

        Lvec[jj] = L
    end

    Lmat = [Lmat[jj][jpx] for jpx in eachindex(Lmat[1]), jj in eachindex(Lmat)]

    xvec = reformat_x(Xdict[k], k)

    yaxis_title = "%"
    if share
        plans .*= 1 / annualized(πN)
        yaxis_title = "Share of Nash inflation"
    end

    if showL
        return plot(scatter(x=xvec, y=Lvec), Layout(template=qtemplate(; slides, dark), xaxis_title="<i>$(string(k))", title="lim<sub><i>p→0</i></sub> min<sub>(a,ω,χ)</sub> 𝓛(<i>p,a,ω,χ)"))
    end

    if showLmat
        @assert size(Lmat) == (length(pgrid), length(xvec))
        return plot(
                contour(z = Lmat', x = pgrid, y = xvec, colorscale = [[jj, get(ColorSchemes.davos, 1-jj, :clamp)] for jj in range(0,1,length=50)]
            ), Layout(template=qtemplate(;slides,dark), title="lim<sub><i>p→0</i></sub> min<sub>(a,ω,χ)</sub> 𝓛(<i>p,a,ω,χ)", xaxis_title="<i>p", yaxis_title = "<i>$(string(k))"))
    end

    plot(scatscol(plans, 0:T-1, xvec, name=string(k)), Layout(;
        template=qtemplate(; slides, dark), xaxis_title="<i>Quarters", yaxis_range=[-0.01, maximum(plans) * 1.05], yaxis_title
    )
    )


end
