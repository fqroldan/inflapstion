function Prequel(mt::MultiType)
    ωgrid = mt.ωgrid
    χgrid = mt.χgrid
    agrid = mt.ct.agrid
    pgrid = mt.ct.pgrid
    Agrid = agrid

    L = zeros(length(ωgrid), length(χgrid), length(agrid), length(pgrid), length(Agrid))
    G = zeros(size(L))

    return Prequel(ωgrid, χgrid, agrid, pgrid, Agrid, L, G)
end

function optim_step(ct::Plan, itp_gπ, itp_L, itp_C, aprime::Number)
	gπ = zeros(size(ct.gπ))
	L  = zeros(size(ct.L))
	πN = Nash(ct)

	h = 0.05 * πN
	Gmin, Gmax = -h, πN + h

    broken = 0
    tot = 0
	
	apgrid = gridmake(1:ct.Np, 1:ct.Na)
	Threads.@threads for js in axes(apgrid,1)
		jp, ja = apgrid[js, :]
		pv, av = ct.pgrid[jp], ct.agrid[ja]
        
		ge = ct.gπ[jp, ja]
        
		xguess = [ge, aprime]
        
		πe′ = exp_π_prime(ct, pv, av, itp_gπ, ge, aprime)
        
        obj(G) = (G - opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, G, πe′, use_ϕ = false)[1])^2
        
        res = Optim.optimize(obj, Gmin, Gmax, GoldenSection())
        
        tot += 1
        if sqrt(res.minimum) < 5e-4
            Gc = res.minimizer
        else
            broken += 1
            Gc = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′, use_ϕ = false)[1]
        end

        _, Lc, _ = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, Gc, πe′, use_ϕ = false)

		gπ[jp, ja] = Gc
		L[jp, ja] = Lc

	end
    share = broken / tot
	return gπ, L, share
end

function pf_iter(ct::Plan, aprime::Number)
	knts = (ct.pgrid, ct.agrid);
	itp_gπ = interpolate(knts, ct.gπ, Gridded(Linear()));
	itp_L  = interpolate(knts, ct.L, Gridded(Linear()));
	itp_C  = interpolate(knts, ct.C, Gridded(Linear()));

	new_gπ, new_L, share = optim_step(ct, itp_gπ, itp_L, itp_C, aprime)

	return new_gπ, new_L, share
end

function solve_t0(mt::MultiType)

    m0 = Prequel(mt)
    ct = mt.ct

    finalshare = 0.0
    
    iter = 0
	tot  = length(m0.ωgrid) * length(m0.χgrid) * length(m0.agrid)
    for (jω, ωv) in enumerate(m0.ωgrid), (jχ, χv) in enumerate(m0.χgrid), (ja, av) in enumerate(m0.agrid)
        iter += 1

        ct.L  .= mt.L_mat[jω, jχ, :, :]
        ct.gπ .= mt.g_mat[jω, jχ, :, :]
        
        ct.ω = ωv
        ct.χ = χv
        aprime = av
        
		update_ga!(ct)
        G, L, share = pf_iter(ct, aprime)

        finalshare = finalshare * (iter-1)/iter + share / iter
        
        m0.L[jω, jχ, ja, :, :] .= L
        m0.G[jω, jχ, ja, :, :] .= G
        
        print("Plan $iter of $tot: share > tol = $(@sprintf("%.3g", 100*share))%. So far share = $(@sprintf("%.3g", 100*finalshare))%\n")
    end

    return m0
end