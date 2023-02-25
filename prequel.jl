function Prequel_s(mt::MultiType; jp = 2)

    _, jj = findmin(mt.L_mat[:,:,jp,:])

    χ_star = mt.χgrid[jj[2]]

    Agrid = [0., χ_star]

    return Prequel(mt, Agrid = Agrid)
end

function Prequel(mt::MultiType; Agrid = mt.ct.agrid)
    ωgrid = mt.ωgrid
    χgrid = mt.χgrid
    agrid = mt.ct.agrid
    pgrid = mt.ct.pgrid

    L = zeros(length(ωgrid), length(χgrid), length(agrid), length(pgrid), length(Agrid))
    G = zeros(size(L))

    return Prequel(ωgrid, χgrid, agrid, pgrid, Agrid, L, G)
end

Np(m::Prequel) = length(m.pgrid)
NA(m::Prequel) = length(m.Agrid)

function optim_step(ct::Plan, m0::Prequel, itp_gπ, itp_L, itp_C, aprime::Number)
	gπ = zeros(Np(m0), NA(m0))
    L = similar(gπ)
    # gπ = zeros(size(ct.gπ))
	# L  = zeros(size(ct.L))
	πN = Nash(ct)

	h = 0.05 * πN
	Gmin, Gmax = -h, πN + h

    broken = 0
    tot = 0
	
	apgrid = gridmake(1:Np(m0), 1:NA(m0))
	Threads.@threads for js in axes(apgrid,1)
		jp, jA = apgrid[js, :]
		pv, av = m0.pgrid[jp], m0.agrid[jA]
        
		ge = πN
        
		xguess = [ge, aprime]
        
		πe′ = exp_π_prime(ct, pv, av, itp_gπ, ge, aprime)
        
        obj(G) = (G - opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, G, πe′, use_ϕ = false)[1])^2
        
        res = Optim.optimize(obj, Gmin, Gmax, GoldenSection())
        
        tot += 1
        if Optim.converged(res) && sqrt(res.minimum) < 5e-4
            Gc = res.minimizer
        else
            broken += 1
            Gc = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′, use_ϕ = false)[1]
        end

        _, Lc, _ = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, Gc, πe′, use_ϕ = false)

		gπ[jp, jA] = Gc
		L[jp, jA] = Lc

	end
    share = broken / tot
	return gπ, L, share
end


function solve_t0(mt::MultiType)

    m0 = Prequel_s(mt)
    ct = mt.ct

    finalshare = 0.0
    
    iter = 0
	tot  = length(m0.ωgrid) * length(m0.χgrid) * length(m0.agrid)
    for (jω, ωv) in enumerate(m0.ωgrid), (jχ, χv) in enumerate(m0.χgrid), (ja, av) in enumerate(m0.agrid)
        iter += 1

        knts   = (mt.ct.pgrid, mt.ct.agrid)
        itp_L  = interpolate(knts, mt.L_mat[jω, jχ, :, :], Gridded(Linear()))
        itp_gπ = interpolate(knts, mt.g_mat[jω, jχ, :, :], Gridded(Linear()))
    	itp_C  = interpolate(knts, mt.C_mat[jω, jχ, :, :], Gridded(Linear()));
    
        ct.ω = ωv
        ct.χ = χv
        aprime = av
        
		update_ga!(ct)
        G, L, share = optim_step(ct, m0, itp_gπ, itp_L, itp_C, aprime)

        finalshare = finalshare * (iter-1)/iter + share / iter
        
        m0.L[jω, jχ, ja, :, :] .= L
        m0.G[jω, jχ, ja, :, :] .= G
        
        print("Plan $iter of $tot: share > tol = $(@sprintf("%.3g", 100*share))%. So far share = $(@sprintf("%.3g", 100*finalshare))%\n")
    end

    return m0
end