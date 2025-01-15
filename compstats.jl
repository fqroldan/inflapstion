function compstats(mt::MultiType, xvec::AbstractVector, k::Symbol)

    ωvec = similar(xvec)
    χvec = similar(xvec)
    avec = similar(xvec)

    for (jx, xv) in enumerate(xvec)

        mt.ct.pars[k] = xv

        print("Starting with $(string(k)) = $(@sprintf("%.3g", xv))\n")
        solve_all!(mt, check = true, tol = 5e-4)

        _, jj = findmin(mt.L_mat[:,:,2,:])

        jω = jj[1]
	    jχ = jj[2]
	    ja = jj[3]

        ωvec[jx] = mt.ωgrid[jω]
        χvec[jx] = mt.χgrid[jχ]
        avec[jx] = mt.ct.gr[:a][ja]
    end

    return ωvec, χvec, avec
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
