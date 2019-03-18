using PlotlyJS

col = [	"#1f77b4",  # muted blue
		"#ff7f0e",  # safety orange
		"#2ca02c",  # cooked asparagus green
		"#d62728",  # brick red
		"#9467bd",  # muted purple
		"#8c564b",  # chestnut brown
		"#e377c2",  # raspberry yogurt pink
		"#7f7f7f",  # middle gray
		"#bcbd22",  # curry yellow-green
		"#17becf"   # blue-teal
		]

function lines(ct::CrazyType, y_mat; dim::Int64=0, title::String="", showleg::Bool=false)
	if dim == 1
		xgrid = ct.pgrid
		zgrid = ct.agrid
		xtitle= "𝑝"
	elseif dim == 2
		xgrid = ct.agrid
		zgrid = ct.pgrid
		xtitle= "𝑎"
	else
		throw(error("wrong dim"))
	end
	Nz = length(zgrid)
	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, Nz)
	for (jz, zv) in enumerate(zgrid)
		if dim == 1
			y_vec = y_mat[:, jz]
			name = "𝑎"
		elseif dim == 2
			y_vec = y_mat[jz, :]
			name = "𝑝"
		end
		name = name * " = $(@sprintf("%.2g", zv))"
		jz % 2 == 0 ? showleg_i = showleg : showleg_i = false
		l_new = scatter(;x=xgrid, y=y_vec, name = name, showlegend = showleg_i, marker_color=col[ceil(Int,10*jz/Nz)])
		l[jz] = l_new
	end
	p = plot([l[jz] for jz in 1:Nz], Layout(;title=title, xaxis_title=xtitle))
	return p
end

function plot_ct(ct::CrazyType, y_tuple, n_tuple; make_pdf::Bool=false, make_png::Bool=false)
	if length(y_tuple) != length(n_tuple)
		throw(error("Make sure # y's = # n's"))
	end

	N = length(y_tuple)
	pl = Array{PlotlyJS.SyncPlot,2}(undef, N, 2)
	for jj in 1:N
		pl[jj, 1] = lines(ct, y_tuple[jj], dim = 1, title=n_tuple[jj], showleg = (jj==1))
		pl[jj, 2] = lines(ct, y_tuple[jj], dim = 2, title=n_tuple[jj], showleg = (jj==1))
	end

	# p = hvcat(2, pl[:])

	relayout!(p, font_family = "Fira Sans Light", font_size = 12, plot_bgcolor="rgba(250, 250, 250, 1.0)", paper_bgcolor="rgba(250, 250, 250, 1.0)")

	function makeplot(p, ext::String)
		savefig(p, pwd() * "/../Graphs/ct" * ext)
	end

	if make_pdf
		makeplot(p, ".pdf")
	end
	if make_png
		makeplot(p, ".png")
	end

	return p
end

annualized(π::Float64) = 100*((1.0 .+ π).^4 .- 1)

function plot_ct_pa(ct::CrazyType, y=ct.L, name="𝓛"; ytitle="")

	a_max = Nash(ct)
	jamax = findfirst(ct.agrid.>=a_max)

	function set_col(ja, agrid, rel::Bool=false)
		if rel
			return ceil(Int,1+9*(agrid[ja])/a_max)
		else
			return ceil(Int,10*ja/jamax)
		end
	end

	p1 = plot([
		scatter(;x=ct.pgrid, y=y[:,ja], marker_color=col[set_col(ja,ct.agrid)], name = "a=$(@sprintf("%.3g", annualized(av)))") for (ja,av) in enumerate(ct.agrid) if av <= a_max
		], Layout(;title=name, fontsize=20,font_family="Fira Sans Light", xaxis_zeroline=false, xaxis_title= "𝑝", yaxis_title=ytitle))
	return p1
end

function makeplots_ct(ct::CrazyType; make_pdf::Bool=false, make_png::Bool=false)

	gπ_over_a = zeros(size(ct.gπ))
	Ep_over_p = zeros(size(ct.Ep))
	for (jp, pv) in enumerate(ct.pgrid), (ja,av) in enumerate(ct.agrid)
		gπ_over_a[jp, ja] = ct.gπ[jp, ja] - av
		Ep_over_p[jp, ja] = ct.Ep[jp, ja] - pv
	end

	p1 = plot_ct(ct, (ct.gπ, ct.L), ("gπ", "𝓛"); make_pdf=make_pdf, make_png=make_png)

	p2 = plot_ct(ct, (ct.Ey, ct.Eπ), ("𝔼y", "𝔼π"); make_pdf=make_pdf, make_png=make_png)

	p3 = plot_ct(ct, (gπ_over_a, Ep_over_p), ("gπ-a", "𝔼p'-p"); make_pdf=make_pdf, make_png=make_png)

	return p1, p2, p3
end

function makeplots_ct_pa(ct::CrazyType)

	gπ_over_a = zeros(size(ct.gπ))
	Ep_over_p = zeros(size(ct.Ep))
	for (jp, pv) in enumerate(ct.pgrid), (ja,av) in enumerate(ct.agrid)
		gπ_over_a[jp, ja] = ct.gπ[jp, ja] - av
		Ep_over_p[jp, ja] = ct.Ep[jp, ja] - pv
	end

	annual_π = (1 .+ ct.gπ).^4 .- 1

	pL = plot_ct_pa(ct, ct.L, "𝓛")
	pπ = plot_ct_pa(ct, 100*annual_π, "gπ", ytitle="%")
	py = plot_ct_pa(ct, ct.Ey, "𝔼y")
	pp = plot_ct_pa(ct, Ep_over_p, "𝔼p'-p")

	p = [pL pπ; py pp]

	relayout!(p, font_family = "Fira Sans Light", font_size = 12, plot_bgcolor="rgba(250, 250, 250, 1.0)", paper_bgcolor="rgba(250, 250, 250, 1.0)")

	return p
end
