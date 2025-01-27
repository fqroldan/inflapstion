function compstats(mt::MultiType, xvec::AbstractVector, k::Symbol)

    ωvec = similar(xvec)
    χvec = similar(xvec)
    avec = similar(xvec)
    Lmat = zeros(length(xvec), length(mt.ct.gr[:p]))

    iter = enumerate(xvec)
    for (jx, xv) in iter

        mt.ct.pars[k] = xv

        println("Starting with $(string(k)) = $(@sprintf("%.3g", xv)) at $(Dates.format(now(), "HH:MM"))")
        solve_all!(mt, check = true, tol = 5e-4, tinyreport = true)

        L, jj = findmin(mt.L_mat[:,:,2,:])

        jω = jj[1]
	    jχ = jj[2]
	    ja = jj[3]

        ωvec[jx] = mt.ωgrid[jω]
        χvec[jx] = mt.χgrid[jχ]
        avec[jx] = mt.ct.gr[:a][ja]

        for jp in axes(Lmat, 2)
            Lvec[jx, jp] = findmin(mt.L_mat[:,:,jp,:])[1]
        end
    end

    return ωvec, χvec, avec, Lmat, xvec
end

function compstats(mt::MultiType, k::Symbol, K=15)
    if k == :σ
        smin, smax = 0.75/400, 1.25/400
    elseif k == :β
        smin, smax = 1.01^(-0.25), 1.1^(-0.25)
    elseif k == :κ
        smin, smax = 0.1, 0.25
    else
        throw(error("Wrong parameter name"))
    end

    svec = range(smin, smax, length=K)
    compstats(mt, svec, k)
end

function sk(k::Symbol)
    s = string(k)
    if s == "σ"
        s = "sigma"
    elseif s == "β"
        s = "beta"
    elseif s == "κ"
        s = "kappa"
    end
    return s
end

function save_compstats(k::Symbol, K = 15)
    mt = load("Output/JET/mt.jld2", "mt")

    s = sk(k)

    ωvec, χvec, avec, Lmat, xvec, pvec = compstats(mt, k, K)

    save("Output/JET/compstats_$s.jld2", "ωvec", ωvec, "χvec", χvec, "avec", avec, "Lmat", Lmat, "$(s)vec", xvec, "pvec", pvec)
end