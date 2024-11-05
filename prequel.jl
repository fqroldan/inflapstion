function Prequel_s(mt::MultiType; jp = 2, Np = 10)

    _, jj = findmin(mt.L_mat[:,:,jp,:])

    χ_star = mt.χgrid[jj[2]]

    Agrid = [0., χ_star]

	pgrid = cdf.(Beta(6,1), range(0,1,length=100))[1:Np]

    return Prequel(mt, Agrid = Agrid)#, pgrid = pgrid)
end

function Prequel(mt::MultiType; Agrid = mt.ct.agrid, pgrid = mt.ct.pgrid)
    ωgrid = mt.ωgrid
    χgrid = mt.χgrid
    agrid = mt.ct.agrid

    L = zeros(length(ωgrid), length(χgrid), length(agrid), length(pgrid), length(Agrid))
    G = zeros(size(L))

    return Prequel(ωgrid, χgrid, agrid, pgrid, Agrid, L, G)
end

function makeZero(mt::MultiType; Na = length(mt.ct.gr[:a]), pgrid = mt.ct.gr[:p], T = 0)


    A = .8 * Nash(mt.ct)
    agrid = cdf.(Beta(0.5,1.5), range(0,1,length=Na))
	move_grids!(agrid, xmax=A, xmin=0.0)

    pars = Dict(k => v for (k,v) in mt.ct.pars)
    gr = Dict(:p => pgrid, :a => agrid)
    
    K = T + 2

    L = zeros(length(gr[:p]), length(gr[:a]), (length(gr[:a]) for _ in 1:T)...)

    G = similar(L)

    return Zero{K}(pars, gr, L, G)
end

function makeZero(z0::Zero{T}) where T

    pars = Dict(k => v for (k, v) in z0.pars)
    gr = Dict(k => copy(z0.gr[k]) for k in keys(z0.gr))

    K = T + 1

    L = zeros(length(gr[:p]), length(gr[:a]), (length(gr[:a]) for _ in 1:(K-2))...)

    G = similar(L)

    return Zero{K}(pars, gr, L, G)
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
		pv, av = m0.pgrid[jp], m0.Agrid[jA]
        
		ge = πN
        
		xguess = [ge, aprime]
        
		πe′ = exp_π_prime(ct, pv, av, itp_gπ, ge, aprime)
        
        obj(G) = (G - opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, G, πe′, use_ϕ = false)[1])^2
        
        res = Optim.optimize(obj, Gmin, Gmax, GoldenSection())
        
        tot += 1
        if Optim.converged(res) && sqrt(res.minimum) < 5e-4
            Gc = res.minimizer
        else
            res = Optim.optimize(obj, Gmin, Gmax, Brent())
            if Optim.converged(res) && sqrt(res.minimum) < 5e-4
                Gc = res.minimizer
            else
                broken += 1
                Gc = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, ge, πe′, use_ϕ = false)[1]
            end
        end

        _, Lc, _ = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, av, Gc, πe′, use_ϕ = false)

		gπ[jp, jA] = Gc
		L[jp, jA] = Lc

	end
    share = broken / tot
	return gπ, L, share
end


function initZero(ct::CrazyType; Np = 80, Na = 45)

    pars = ct.pars
    gr = Dict(:p => ct.gr[:p], :a => ct.gr[:a])

    L = zeros(length(gr[:p]), length(gr[:a]))
    G = similar(L)

    return Zero{2}(pars, gr, L, G)
end

preperiods(m0::Zero) = length(size(m0.L)) - 2

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

periods(z::Zero{T}) where T = T

QuantEcon.gridmake(z::Zero{2}) = gridmake(1:length(z.gr[:a]))
QuantEcon.gridmake(z::Zero{3}) = gridmake(1:length(z.gr[:a]), 1:length(z.gr[:a]))
QuantEcon.gridmake(z::Zero{4}) = gridmake(1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]))
# QuantEcon.gridmake(z::Zero{5}) = gridmake(1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]))
# QuantEcon.gridmake(z::Zero{6}) = gridmake(1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]))
# QuantEcon.gridmake(z::Zero{7}) = gridmake(1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]), 1:length(z.gr[:a]))

getnext(y::Array{<:Number,2}, _) = y[:, :]
getnext(y::Array{<:Number,3}, js) = y[:, js[1], :]
getnext(y::Array{<:Number,4}, js) = y[:, js[1], js[2], :]
# getnext(y::Array{<:Number,5}, js) = y[:, js[1], js[2], js[3], :]
# getnext(y::Array{<:Number,6}, js) = y[:, js[1], js[2], js[3], js[4], :]
# getnext(y::Array{<:Number,7}, js) = y[:, js[1], js[2], js[3], js[4], js[5], :]

fillz!(y::Array{<:Number,3}, jp, ja, j0, x) = (y[jp, ja[1], j0] = x)
fillz!(y::Array{<:Number,4}, jp, ja, j0, x) = (y[jp, ja[1], ja[2], j0] = x)
fillz!(y::Array{<:Number,5}, jp, ja, j0, x) = (y[jp, ja[1], ja[2], ja[3], j0] = x)

function solve!(z::Zero{2})

    
    πN = Nash(z)

    Np = length(z.gr[:p])
    for (ja, a) in enumerate(z.gr[:a])

        print("Running ct with a/πN = $(@sprintf("%.3g", 100*a/πN))%: ")

        ct = CCT(z.pars; a, Np)

        flag = pfi!(ct, verbose = false)

        if flag
            print("✓\n")
        else
            print("WARNING: no convergence\n")
        end

        z.L[:, ja] .= ct.L[:, 1]
        z.G[:, ja] .= ct.gπ[:, 1]
    end
    nothing
end

function extend(z0::Zero)
    z = makeZero(z0)
    πN = Nash(z0)
    
    gguess = 0.25 * πN
    h = 0.05 * πN
    Gmin, Gmax = -h, πN + h
    
    ct = CCT(z.pars)

    S = gridmake(z0)
    for jj in axes(S,1)

        ja = S[jj, :]
        avec = [z0.gr[:a][jj] for jj in ja]
        
        print("Out loop $jj with (a) = $([round(x/Nash(z), sigdigits=3) for x in avec])")

        aprime = avec[end]
        
        Lp = getnext(z0.L, ja)
        gp = getnext(z0.G, ja)

        knts = (z0.gr[:p], z0.gr[:a])
        itp_L  = interpolate(knts, Lp, Gridded(Linear()))
        itp_gπ = interpolate(knts, gp, Gridded(Linear()))
        itp_C  = interpolate(knts, 0*gp, Gridded(Linear()))

        counter = 0
        for (jp, pv) in enumerate(z.gr[:p]), (j0, a0) in enumerate(z.gr[:a])

            xguess = [gguess, aprime]

            πe′ = 0

            obj(G) = (G - opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, a0, G, πe′)[1])^2

            res = Optim.optimize(obj, Gmin, Gmax, GoldenSection())

            if Optim.converged(res) && sqrt(res.minimum) < 5e-4
            else
                counter += 1
            end

            Gc = res.minimizer
            _, Lc, _ = opt_L(ct, itp_gπ, itp_L, itp_C, xguess, pv, a0, Gc, πe′, use_ϕ=false)

            fillz!(z.G, jp, ja, j0, Gc)
            fillz!(z.L, jp, ja, j0, Lc)
        end
        success_rate = 1 - counter / (length(z.gr[:p]) * length(z.gr[:a]))
        print(": $(@sprintf("%.3g",100 * success_rate))% successful\n")
    end
    return z
end
