function iter_simul(ct::Plan, itp_gπ, itp_ga, pv, av; noshocks::Bool=false)
	ϵ = rand(dist_ϵ(ct))
	if noshocks
		ϵ = 0.0
	end

	exp_π = itp_gπ(pv, av)
	obs_π = exp_π+ϵ
	
	pprime = Bayes(ct, obs_π, exp_π, pv, av)

	# aprime = ϕ(ct, av)
	aprime = itp_ga(pv, av)
	exp_π′ = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)

	y = PC(ct, obs_π, exp_π, exp_π′)

	return pprime, aprime, obs_π, y, exp_π, ϵ
end

function simul(ct::Plan; T::Int64=50, jp0::Int64=3, noshocks::Bool=false)
	p0 = ct.pgrid[jp0]

	_, ind_a0 = findmin(ct.L[jp0, :])

	# ind_a0 = floor(Int, ct.Na*0.25)

	a0 = ct.agrid[ind_a0]
	p, a = p0, a0

	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()))
	itp_ga = interpolate(knots, ct.ga, Gridded(Linear()))
	itp_L  = interpolate(knots, ct.L, Gridded(Linear()))

	p_vec, a_vec, π_vec, y_vec, g_vec, ϵ_vec, L_vec = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
	for tt = 1:T
		p_vec[tt], a_vec[tt] = p, a
		L_vec[tt] = itp_L(p, a)
		pp, ap, πt, yt, gt, ϵt = iter_simul(ct, itp_gπ, itp_ga, p, a; noshocks=noshocks)
		π_vec[tt], y_vec[tt], g_vec[tt], ϵ_vec[tt] = πt, yt, gt, ϵt

		p, a = pp, ap
	end

	σϵ = std(ϵ_vec)
	# print("\nStd of shocks = $(@sprintf("%.3g", σϵ))")

	return p_vec, a_vec, π_vec, y_vec, g_vec, L_vec
end


function plot_announcement(pp::DovisKirpalani, T::Int64 = 4*5; K::Int64=16)

	pvec = range(minimum(pp.pgrid), maximum(pp.pgrid), length=K)

	avec = zeros(K, T)
	πvec = zeros(K, T)

	knots = (pp.pgrid, pp.agrid)
	itp_L  = interpolate(knots, pp.L, Gridded(Linear()))
	itp_gπ = interpolate(knots, pp.gπ, Gridded(Linear()))
	itp_ga = interpolate(knots, pp.ga, Gridded(Linear()))

	a0 = [findmin(pp.L[jp,:])[1] for (jp,pv) in enumerate(pp.pgrid)]

	for jt in 1:T

		avec[:, jt] .= a0

		# WHAT DO YOU DO WITH MOVEMENTS IN P ???

	end

	return avec, πvec
end		
