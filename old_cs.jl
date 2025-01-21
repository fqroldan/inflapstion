function plot_cs(k::Symbol; T = 11, showL = false, slides = true, dark = false)

    s = sk(k)

    K = ifelse(k == :β, 19, 20)

    Xdict = Dict(k => zeros(K)*NaN for k in (:β, :σ, :κ))
    ωvec = zeros(K) * NaN
    χvec = zeros(K) * NaN
    avec = zeros(K) * NaN

    Lvec = zeros(K)

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

        L, jmin = findmin(mt.L_mat[:, :, 2, :])

        jω = jmin[1]
        jχ = jmin[2]
        ja = jmin[3]

        Ω = exp(-mt.ωgrid[jω])
        χ = annualized(mt.χgrid[jχ])
        a = annualized(mt.ct.agrid[ja])
        ωvec[jj] = mt.ωgrid[jω]
        χvec[jj] = mt.χgrid[jχ]
        avec[jj] = mt.ct.agrid[ja]
        
        for tt in 1:T
            plans[tt, jj] = χ + Ω^(tt - 1) * (a - χ)
        end

        Lvec[jj] = L
    end

    xvec = reformat_x(Xdict[k], k)

    if showL
        return plot(scatter(x=xvec, y=Lvec), Layout(template=qtemplate(; slides, dark), xaxis_title="<i>$(string(k))", title="lim<sub><i>p→0</i></sub> min<sub>(a,ω,χ)</sub> 𝓛(<i>p,a,ω,χ)"))
    end

    plot(scatscol(plans, 0:T-1, xvec, name=string(k)), Layout(;
        template=qtemplate(; slides, dark), xaxis_title="<i>Quarters", yaxis_range=[-0.01, maximum(plans) * 1.05]
    )
    )


end
