function compstats(mt::MultiType, xvec::AbstractVector, k::Symbol)

    ωvec = similar(xvec)
    χvec = similar(xvec)
    avec = similar(xvec)

    iter = ProgressBar(enumerate(xvec))
    for (jx, xv) in iter

        mt.ct.pars[k] = xv

        println(iter, "Starting with $(string(k)) = $(@sprintf("%.3g", xv))")
        solve_all!(mt, iter, check = true, tol = 5e-4)

        _, jj = findmin(mt.L_mat[:,:,2,:])

        jω = jj[1]
	    jχ = jj[2]
	    ja = jj[3]

        ωvec[jx] = mt.ωgrid[jω]
        χvec[jx] = mt.χgrid[jχ]
        avec[jx] = mt.ct.gr[:a][ja]
    end

    return ωvec, χvec, avec, xvec
end

function compstats(mt::MultiType, k::Symbol, K=15)
    if k == :σ
        smin, smax = 0.85/400, 1.25/400
    elseif k == :β
        smin, smax = 1.01^(-0.25), 1.05^(-0.25)
    elseif k == :κ
        smin, smax = 0.1, 0.25
    else
        throw(error("Wrong parameter name"))
    end

    svec = range(smin, smax, length=K)
    compstats(mt, svec, k)
end


function save_compstats(k::Symbol, K = 15)
    mt = load("Output/JET/mt.jld2", "mt")

    s = string(k)
    if s == "σ"
        s = "sigma"
    elseif s == "β"
        s = "beta"
    elseif s == "κ"
        s = "kappa"
    end

    ωvec, χvec, avec, xvec = compstats(mt, k, K)

    save("Output/JET/compstats_$s.jld2", "ωvec", ωvec, "χvec", χvec, "avec", avec, "$(s)vec", xvec)
end
