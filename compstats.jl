function compstats(mt::MultiType, xvec::AbstractVector, k::Symbol; fname, saveprog, tol)

    Lm   = Vector{typeof(mt.L_mat)}(undef, length(xvec))

    compstats!(xvec, Lm, k, mt; fname, saveprog, tol)
end

function compstats!(xvec, Lm, k::Symbol, mt::MultiType; fname, saveprog, tol)

    ωvec = similar(xvec)
    χvec = similar(xvec)
    avec = similar(xvec)
    Lmat = zeros(length(xvec), length(mt.ct.gr[:p]))
    
    iter = enumerate(xvec)
    for (jx, xv) in iter

        if isassigned(Lm, jx)
            mt.L_mat .= Lm[jx]
        end

        mt.ct.pars[k] = xv

        println("Starting with $(string(k)) = $(@sprintf("%.3g", xv)) on $(Dates.format(now(), "dd-u at HH:MM"))")
        solve_all!(mt; check = true, tol, tinyreport = true, maxiter = 500)

        L, jj = findmin(mt.L_mat[:,:,2,:])

        jω = jj[1]
	    jχ = jj[2]
	    ja = jj[3]

        ωvec[jx] = mt.ωgrid[jω]
        χvec[jx] = mt.χgrid[jχ]
        avec[jx] = mt.ct.gr[:a][ja]

        Lm[jx] = copy(mt.L_mat)
        for jp in axes(Lmat, 2)
            Lmat[jx, jp] = findmin(mt.L_mat[:,:,jp,:])[1]
        end

        if saveprog
            save(fname, "xvec", xvec, "Lm", Lm, "k", k, "mt", mt)
        end
    end

    return ωvec, χvec, avec, Lmat, xvec, Lm
end

function compstats(mt::MultiType, k::Symbol; K=15, fname="Output/JET/temp.jld2", saveprog=length(fname) > 0, tol=5e-4)
    if k == :σ
        smin, smax = 0.75/400, 1.25/400
    elseif k == :β
        smin, smax = 1.01^(-0.25), 1.09^(-0.25)
    elseif k == :κ
        smin, smax = 0.1, 0.25
    else
        throw(error("Wrong parameter name"))
    end

    svec = range(smin, smax, length=K)
    compstats(mt, svec, k; fname, saveprog, tol)
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

function save_compstats(k::Symbol; K=15, folder="Output/JET/", tol=5e-4)
    mt = load(folder * "mt.jld2", "mt")

    s = sk(k)
    fname = ifelse(length(folder)>0, folder * "temp_$(s).jld2", "")

    ωvec, χvec, avec, Lmat, xvec, pvec = compstats(mt, k; K, fname, tol)

    save(folder * "compstats_$s.jld2", "ωvec", ωvec, "χvec", χvec, "avec", avec, "Lmat", Lmat, "$(s)vec", xvec, "pvec", pvec)
end