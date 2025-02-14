function consolidate(k::Symbol)
    s = sk(k)

    xvec = Float64[]
    Mvec = Vector{MultiType}(undef, 0)

    for jj in 1:20
        if k == :κ
            name = "mt"
        else
            name = "mt$jj"
        end
        try
            mt = load("Output/CompStats/$s/run$jj/$name.jld2", "mt")

            pars = Dict(:β => mt.ct.β, :σ => mt.ct.σ, :κ => mt.ct.κ, :ystar => mt.ct.ystar, :ω => mt.ct.ω, :χ => mt.ct.χ, :γ => mt.ct.γ, :ψ => 0.)
            
            gr = Dict(:a => mt.ct.agrid, :p => mt.ct.pgrid)
            opt = Dict(:use_a => true)

            ct = CrazyType{Forward}(pars, gr, opt, mt.ct.gπ, mt.ct.ga, mt.ct.L, mt.ct.C, mt.ct.Ey, mt.ct.Eπ, mt.ct.Ep)

            G = zeros(size(mt.L_mat))
            mt_new = MultiType(ct, mt.ωgrid, mt.χgrid, mt.z, mt.ν, mt.μ, mt.L_mat, mt.C_mat, G)

            push!(xvec, mt_new.ct.pars[k])
            push!(Mvec, mt_new)
        catch
            print("Could not find mt $name at 'Output/CompStats/$s/run$jj/\n")
        end
    end
    save("Output/CompStats/$s/CS_consol.jld2", "Mvec", Mvec, "$(s)vec", xvec)
end

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


function plot_consol(k::Symbol; T=11, jp=2, showLmat=false, showL=false, slides=true, dark=false, share=false)
    s = sk(k)

    Mvec, xvec = load("Output/CompStats/$s/CS_consol.jld2", "Mvec", "$(s)vec")
    
    K = length(xvec)

    plans = zeros(T, K)
    Lvec = zeros(K)
    Lmat = Vector{Vector{Float64}}()
    pgrid = Vector{Float64}()
    πN = 0.

    for jj in eachindex(xvec)
        mt = Mvec[jj]
        L, jmin = findmin(mt.L_mat[:, :, jp, :])
        jω = jmin[1]
        jχ = jmin[2]
        ja = jmin[3]

        lm = [findmin(mt.L_mat[:, :, jpx, :])[1] for jpx in 2:N(mt.ct, :p)-1]
        push!(Lmat, lm)
        pgrid = mt.ct.gr[:p][2:end-1]

        Ω = exp(-mt.ωgrid[jω])
        χ = annualized(mt.χgrid[jχ])
        a = annualized(mt.ct.gr[:a][ja])
        # ωvec[jj] = mt.ωgrid[jω]
        # χvec[jj] = mt.χgrid[jχ]
        # avec[jj] = mt.ct.gr[:a][ja]

        πN = Nash(mt)

        for tt in 1:T
            plans[tt, jj] = χ + Ω^(tt - 1) * (a - χ)
        end

        Lvec[jj] = L
    end

    Lmat = [Lmat[jj][jpx] for jpx in eachindex(Lmat[1]), jj in eachindex(Lmat)]
    
    xvec = reformat_x(xvec, k)
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
            contour(z=Lmat', x=pgrid, y=xvec, colorscale=[[jj, get(ColorSchemes.davos, 1 - jj, :clamp)] for jj in range(0, 1, length=50)]
            ), Layout(template=qtemplate(; slides, dark), title="lim<sub><i>p→0</i></sub> min<sub>(a,ω,χ)</sub> 𝓛(<i>p,a,ω,χ)", xaxis_title="<i>p", yaxis_title="<i>$(string(k))"))
    end

    suff = ""
    name = string(k)
    if k == :β
        suff = "%"
        name = "β<sup>-1</sup>-1"
    end


    plot(scatscol(plans, 0:T-1, xvec; name, suff), Layout(;
        template=qtemplate(; slides, dark), xaxis_title="<i>Quarters", yaxis_range=[-0.01, maximum(plans) * 1.05], yaxis_title
    )
    )
end