function iter_simul(ct::CrazyType, itp_gπ, pv, av; noshocks::Bool=false)
	ϵ = rand(dist_ϵ(ct))
	if noshocks
		ϵ = 0.0
	end

	exp_π = itp_gπ(pv, av)
	obs_π = exp_π+ϵ
	
	pprime = Bayes(ct, obs_π, exp_π, pv, av)

	aprime = ϕ(ct, av)
	exp_π′ = pprime * aprime + (1.0-pprime) * itp_gπ(pprime, aprime)

	y = NKPC(ct, obs_π, exp_π′)

	return pprime, aprime, obs_π, y, exp_π, ϵ
end

function simul(ct::CrazyType; T::Int64=50, jp0::Int64=2, noshocks::Bool=false)
	p0 = ct.pgrid[jp0]

	_, ind_a0 = findmin(ct.L[3, :])
	a0 = ct.agrid[ind_a0]

	p, a = p0, a0

	knots = (ct.pgrid, ct.agrid)
	itp_gπ = interpolate(knots, ct.gπ, Gridded(Linear()))
	itp_L  = interpolate(knots, ct.L, Gridded(Linear()))

	p_vec, a_vec, π_vec, y_vec, g_vec, ϵ_vec, L_vec = zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T), zeros(T)
	for tt = 1:T
		p_vec[tt], a_vec[tt] = p, a
		L_vec[tt] = itp_L(p, a)
		pp, ap, πt, yt, gt, ϵt = iter_simul(ct, itp_gπ, p, a; noshocks=noshocks)
		π_vec[tt], y_vec[tt], g_vec[tt], ϵ_vec[tt] = πt, yt, gt, ϵt

		p, a = pp, ap
	end

	σϵ = std(ϵ_vec)
	# print("\nStd of shocks = $(@sprintf("%.3g", σϵ))")

	return p_vec, a_vec, π_vec, y_vec, g_vec, L_vec
end
