using PlotlyJS, ColorSchemes, Printf
include("type_def.jl")

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


""" Define styles """
def_style = let
	axis = attr(showgrid = true, gridcolor="#e2e2e2", gridwidth=0.5, zeroline=false)
	layout = Layout(xaxis = axis, yaxis=axis)
	Style(layout=layout)
end

slides_def = let
	layout = Layout(plot_bgcolor="#fafafa", paper_bgcolor="#fafafa",
		width=1920*0.45, height=1080*0.45, font_size=16, font_family="Lato",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

dark_bg = let
	axis = attr(gridcolor="#353535")
	layout = Layout(plot_bgcolor="#1e1e1e", paper_bgcolor="#1e1e1e", font_color="white", xaxis=axis,yaxis=axis)
	Style(layout=layout)
end
slides_dark = Style(slides_def, dark_bg)

paper = let
	layout = Layout(width = 1920 * 0.5, height = 1080 * 0.35, font_size=16, font_family = "Linux Libertine",
		legend = attr(orientation = "h", x=0.05))
	Style(def_style, layout=layout)
end

function style_plot!(pl; slides::Bool=true)
	if slides
		relayout!(pl,
			font_family = "Lato", font_size=14,
			paper_bgcolor="#fafafa", plot_bgcolor="#fafafa"
			)
	else
		relayout!(pl,
			font_family = "Linux Libertine", font_size=14,
			paper_bgcolor="white", plot_bgcolor="white"
			)
	end
	nothing
end

function lines(ct::Plan, y_mat; dim::Int64=0, title::String="", showleg::Bool=false)
	if dim == 1
		xgrid = ct.pgrid
		zgrid = ct.agrid
		xtitle= "<i>p</i>"
	elseif dim == 2
		xgrid = ct.agrid
		zgrid = ct.pgrid
		xtitle= "<i>a"
	else
		throw(error("wrong dim"))
	end
	Nz = length(zgrid)
	cvec = range(col[1], stop=col[2], length=Nz)
	l = Array{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, Nz)
	for (jz, zv) in enumerate(zgrid)
		if dim == 1
			y_vec = y_mat[:, jz]
			name = "<i>a"
		elseif dim == 2
			y_vec = y_mat[jz, :]
			name = "<i>p</i>"
		end
		name = name * " = $(@sprintf("%.2g", zv))"
		jz % 2 == 0 ? showleg_i = showleg : showleg_i = false
		l_new = scatter(;x=xgrid, y=y_vec, name = name, showlegend = showleg_i, marker_color=cvec[jz])
		l[jz] = l_new
	end
	p = plot([l[jz] for jz in 1:Nz], Layout(;title=title, xaxis_title=xtitle))
	return p
end

function plot_ct(ct::Plan, y_tuple, n_tuple; make_pdf::Bool=false, make_png::Bool=false)
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

	relayout!(p, font_family = "Fira Sans Light", font_size = 12, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa")

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

function plot_ct_pa(ct::Plan, y=ct.L, name="𝓛"; ytitle="", reverse_draw::Bool=false, positive_p::Bool=false, few_lines::Bool=false, f_p = false)

	a_max = Nash(ct)
	if maximum(ct.agrid) > a_max
		jamax = findfirst(ct.agrid.>=a_max)
	else
		jamax = length(ct.agrid)
	end

	positive_p ? xvec = ct.pgrid[2:end] : xvec = ct.pgrid
	positive_p ? y = y[2:end, :] : nothing

	few_lines ? step_a = 2 : step_a = 1

	colorpal = ColorSchemes.lapaz

	function set_col(ja, agrid, rel::Bool=false)
		weight = min(1,max(0,(ja-1)/(jamax-1)))
		if rel
			weight = min(1, agrid[ja] / a_max)
		end
		# return weighted_color_mean(weight, parse(Colorant,col[4]), parse(Colorant,col[1]))
		return get(colorpal, weight)
	end

	xs = [0.5 for ja in length(ct.agrid):-1:1 if ct.agrid[ja] <= a_max]
	ys = [mean(y[:, ja]) for ja in length(ct.agrid):-1:1 if ct.agrid[ja] <= a_max]
	colscale = [[vv, get(colorpal, 0.9*vv)] for vv in range(0,1,length=100)]
	cols = range(0,1,length=length(xs))
	colnames = round.(range(0,1,length=6), digits=2)


	p1 = plot([
		scatter(mode = "markers", marker_opacity = 0,
				x = xs, y = ys, showlegend=false,
				marker = attr(color=cols, reversescale=false, colorscale=colscale, colorbar = attr(tickvals=range(0,1,length=length(colnames)), title="&nbsp;<i>a/π<sup>N", ticktext=colnames))
				)
		[scatter(;x=xvec, y=y[:,ja], marker_color=set_col(ja,ct.agrid), name = "a=$(@sprintf("%.3g", annualized(ct.agrid[ja])))") for ja in 1:step_a:length(ct.agrid) if ct.agrid[ja] <= a_max]
		], Layout(;title=name, xaxis_zeroline=false, xaxis_title= "<i>p</i>", yaxis_title=ytitle));

	if reverse_draw

		p1 = plot([
			scatter(mode = "markers", marker_opacity = 0,
				x = xs, y = ys, showlegend=false,
				marker = attr(color=cols, reversescale=false, colorscale=colscale, colorbar = attr(tickvals=range(0,1,length=length(colnames)), title="&nbsp;<i>a/π<sup>N", ticktext=colnames))
				)
			[scatter(;x=xvec, y=y[:,ja], marker_color=set_col(ja,ct.agrid), showlegend=false, name = "a=$(@sprintf("%.3g", annualized(ct.agrid[ja])))") for ja in length(ct.agrid):-1:1 if ct.agrid[ja] <= a_max]
			], Layout(;title=name, fontsize=16,font_family="Fira Sans Light", xaxis_zeroline=false, xaxis_title= "<i>p</i>", yaxis_title=ytitle))
	end

	if f_p
		xvec = annualized.(ct.agrid)
		xs = [annualized(mean(ct.agrid)) for jp in eachindex(ct.pgrid)]
		ys = [mean(y[jp, :]) for jp in eachindex(ct.pgrid)]
		cols = range(0,1,length=length(xs))

		p1 = plot([
			scatter(mode = "markers", marker_opacity = 0,
					x = xs, y = ys, showlegend=false,
					marker = attr(color=cols, reversescale=false, colorscale=colscale, colorbar = attr(tickvals=range(0,1,length=length(colnames)), title="&nbsp;&nbsp;<i>p", ticktext=colnames))
					)
			[scatter(;x=xvec, y=y[jp,:], showlegend=false, marker_color=set_col(jp,ct.pgrid), name = "p=$(@sprintf("%.3g", pv))") for (jp, pv) in enumerate(ct.pgrid)]
			], Layout(;title=name, xaxis_zeroline=false, xaxis_title= "<i>a</i>", yaxis_title=ytitle));
	end

	return p1
end

function makeplots_ct(ct::Plan; make_pdf::Bool=false, make_png::Bool=false)

	gπ_minus_a = zeros(size(ct.gπ))
	Ep_minus_p = zeros(size(ct.Ep))
	for (jp, pv) in enumerate(ct.pgrid), (ja,av) in enumerate(ct.agrid)
		gπ_minus_a[jp, ja] = ct.gπ[jp, ja] - av
		Ep_minus_p[jp, ja] = ct.Ep[jp, ja] - pv
	end

	p1 = plot_ct(ct, (ct.gπ, ct.L), ("gπ", "𝓛"); make_pdf=make_pdf, make_png=make_png)

	p2 = plot_ct(ct, (ct.Ey, ct.Eπ), ("𝔼y", "𝔼π"); make_pdf=make_pdf, make_png=make_png)

	p3 = plot_ct(ct, (gπ_minus_a, Ep_minus_p), ("gπ-a", "𝔼[<i>p'-p</i>]"); make_pdf=make_pdf, make_png=make_png)

	return p1, p2, p3
end

function makeplots_ct_pa(ct::Plan; slides::Bool=true)
	""" Currently run for paper on ct.ω = 0.01, ct.χ = 0 """
	gπ_minus_a = zeros(size(ct.gπ))
	Eπ_minus_a = zeros(size(ct.gπ))
	Ep_minus_p = zeros(size(ct.Ep))
	for (jp, pv) in enumerate(ct.pgrid), (ja,av) in enumerate(ct.agrid)
		gπ_minus_a[jp, ja] = ct.gπ[jp, ja] - av
		Eπ_minus_a[jp, ja] = pv*av + (1.0-pv)*ct.gπ[jp, ja] - av
		Ep_minus_p[jp, ja] = ct.Ep[jp, ja] - pv
	end

	if slides
		font = "Lato"
		bgcol = "#fafafa"
		heights = 500 * ones(4)
	else
		font = "Linux Libertine"
		bgcol = "white"
		heights = [400, 350, 400, 500]
	end

	annual_π = annualized.(gπ_minus_a)
	Eπ_a 	 = annualized.(Eπ_minus_a)

	pL = plot_ct_pa(ct, ct.L, "𝓛"; reverse_draw=true)
	pπ = plot_ct_pa(ct, annual_π, "<i>g<sup>⋆</sup> - a", ytitle="%")
	pa = plot_ct_pa(ct, annualized.(ct.ga), "<i>g<sub>a", ytitle="%")
	pE = plot_ct_pa(ct, Eπ_a, "𝔼π-a", ytitle="%")
	py = plot_ct_pa(ct, ct.Ey, "𝔼y")
	pp = plot_ct_pa(ct, Ep_minus_p, "𝔼[<i>p'-p</i>]")
	pC = plot_ct_pa(ct, ct.C, "𝓒")

	p = [pL pπ; py pp]

	relayout!(p, font_family = font, font_size = 16, plot_bgcolor=bgcol, paper_bgcolor=bgcol)
	relayout!(pL, font_family = font, font_size = 16, plot_bgcolor=bgcol, paper_bgcolor=bgcol, width=900, height = heights[1])

	relayout!(pπ, font_family=font, xaxis_title="<i>p", yaxis_title="%", font_size=16, width=900, height=heights[2], plot_bgcolor=bgcol, paper_bgcolor=bgcol, yaxis_range = [-0.01, 0.6])
	restyle!(pπ, showlegend=false)
	relayout!(pa, font_family=font, xaxis_title="<i>p", yaxis_title="%", font_size=16, width=900, height=heights[2], plot_bgcolor=bgcol, paper_bgcolor=bgcol, yaxis_zeroline=false, xaxis_zeroline=false)
	restyle!(pa, showlegend=false)

	relayout!(pp, font_family=font, xaxis_title="<i>p", font_size=16, width=900, height=heights[3], plot_bgcolor=bgcol, paper_bgcolor=bgcol)
	restyle!(pp, showlegend=false)

	return p, pL, pπ, pC, pp, pa
end


function plot_simul(ct::Plan; T::Int64=50, N=10000, jp0::Int64=3, noshocks::Bool=false, CIs::Bool=false)
	# Update simulations codes
	# include("simul.jl")

	p_mat, a_mat, π_mat, y_mat, g_mat, L_mat = zeros(T,N), zeros(T,N), zeros(T,N), zeros(T,N), zeros(T,N), zeros(T,N)

	for jn in 1:N
	    p_vec, a_vec, π_vec, y_vec, g_vec, L_vec = simul(ct; jp0=jp0, T=T, noshocks=noshocks)
	    p_mat[:,jn] = p_vec
	    L_mat[:,jn] = L_vec
	    a_mat[:,jn], π_mat[:,jn], y_mat[:,jn], g_mat[:,jn] = annualized.(a_vec), annualized.(π_vec), annualized.(y_vec), annualized.(g_vec)
	end

	# k = 2
	# quantiles = linspace(0,1, k+2)
	quantiles = [0.25; 0.75]
	k = length(quantiles)
	p_qnt, a_qnt, π_qnt, y_qnt, g_qnt, L_qnt = zeros(T,k), zeros(T,k), zeros(T,k), zeros(T,k), zeros(T,k), zeros(T,k)
	for jk in 1:k
	    for jt in 1:T
	        qnt = quantiles[jk]
	        p_qnt[jt,jk], a_qnt[jt,jk], π_qnt[jt,jk], y_qnt[jt,jk], g_qnt[jt,jk], L_qnt[jt,jk] = quantile(p_mat[jt, :], qnt), quantile(a_mat[jt, :], qnt), quantile(π_mat[jt, :], qnt), quantile(y_mat[jt, :], qnt), quantile(g_mat[jt, :], qnt), quantile(L_mat[jt, :], qnt)
	    end
	end
	p_med, a_med, π_med, y_med, g_med, L_med = vec(median(p_mat, dims=2)), vec(median(a_mat, dims=2)), vec(median(π_mat, dims=2)), vec(median(y_mat, dims=2)), vec(median(g_mat, dims=2)), vec(median(L_mat, dims=2))
	p_avg, a_avg, π_avg, y_avg, g_avg, L_avg = vec(mean(p_mat, dims = 2)), vec(mean(a_mat, dims = 2)), vec(mean(π_mat, dims = 2)), vec(mean(y_mat, dims = 2)), vec(mean(g_mat, dims = 2)), vec(mean(L_mat, dims = 2))

	prep = plot([
			[scatter(;x=(1:T)/4, y=p_qnt[:,jk], showlegend=false, opacity=0.25, line_color=col[1]) for jk in 1:k if CIs]
			scatter(;x=(1:T)/4, y=p_avg, showlegend=false, line_color=col[1])
			scatter(;x=(1:T)/4, y=p_med, showlegend=false, line_color=col[1], line_dash="dashdot")
			], Layout(;title="Reputation", yaxis_zeroline=false, font_family = "Fira Sans Light", font_size = 16))
	ptar = plot([
			[scatter(;x=(1:T)/4, y=a_qnt[:,jk], showlegend=false, opacity=0.25, line_color=col[2]) for jk in 1:k if CIs]
			scatter(;x=(1:T)/4, y=a_avg, showlegend=false, line_color=col[2])
			scatter(;x=(1:T)/4, y=g_avg, showlegend=false, line_color=col[5], line_dash="dot")
			# scatter(;x=(1:T)/4, y=a_med, showlegend=false, line_color=col[2], line_dash="dashdot")
			], Layout(;title="Target", yaxis_zeroline=false, font_family = "Fira Sans Light", font_size = 16))
	pinf = plot([
			[scatter(;x=(1:T)/4, y=π_qnt[:,jk], showlegend=false, opacity=0.25, line_color=col[3]) for jk in 1:k if CIs]
			scatter(;x=(1:T)/4, y=π_avg, showlegend=false, line_color=col[3])
			scatter(;x=(1:T)/4, y=π_med, showlegend=false, line_color=col[3], line_dash="dashdot")
			scatter(;x=(1:T)/4, y=g_avg, showlegend=false, line_color=col[5], line_dash="dot")
			], Layout(;title="Inflation", yaxis_zeroline=false, font_family = "Fira Sans Light", font_size = 16))
	pout = plot([
			[scatter(;x=(1:T)/4, y=y_qnt[:,jk], showlegend=false, opacity=0.25, line_color=col[4]) for jk in 1:k if CIs]
			scatter(;x=(1:T)/4, y=y_avg, showlegend=false, line_color=col[4])
			scatter(;x=(1:T)/4, y=y_med, showlegend=false, line_color=col[4], line_dash="dashdot")
			], Layout(;title="Output", yaxis_zeroline=false, font_family = "Fira Sans Light", font_size = 16))
	p = [prep ptar; pinf pout]

	relayout!(p, font_family = "Fira Sans Light", font_size = 14, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa")

	pL = plot([
		scatter(;x=(1:T)/4, y=L_avg, showlegend=false, line_color=col[4])
		scatter(;x=(1:T)/4, y=L_med, showlegend=false, line_color=col[4], line_dash="dashdot")
		], Layout(;title="𝓛", font_family = "Fira Sans Light", font_size = 16))

    return p, pL, ptar
end

function makeplot_conv(dists::Vector; switch_η=25)
	T = length(dists)

	function MA_t(t::Int64)
		return [100*mean(dists[jt-t:jt]) for jt in (t+1):T]
	end

	shapes = [vline(xchange, line_dash="dot", line_width=1, line_color="black") for xchange in 1:T if xchange%switch_η==0]

	push!(shapes, hline(25e-4*100, line_dash="dash", line_width=1, line_color="black"))
	
	yvec = MA_t(0)
	ls = [scatter(;x=1:T, y=yvec, showlegend=false)]
	push!(shapes, hline(minimum(yvec), line_dash="dash", line_width=1, line_color=col[1]) )
	
	if T > 11
		yvec = MA_t(10)
		push!(ls, scatter(;x=5:T-5, y=yvec, showlegend=false))
		push!(shapes, hline(minimum(yvec), line_dash="dash", line_width=1, line_color=col[2]))
		if T > 51
			yvec = MA_t(50)
			push!(ls, scatter(;x=25:T-25, y=yvec, showlegend=false))
			push!(shapes, hline(minimum(yvec), line_dash="dash", line_width=1, line_color=col[3]))
			if T > 101
				yvec = MA_t(100)
				push!(ls, scatter(;x=50:T-50, y=yvec, showlegend=false))
				push!(shapes, hline(minimum(yvec), line_dash="dash", line_width=1, line_color=col[4]))
			end
		end
	end
	p1 = plot(ls, Layout(;shapes = shapes))

	relayout!(p1, yaxis_type="log", title="‖g′-g‖/‖g‖", xaxis_title="iterations", yaxis_title="%")
	return p1
end

function plot_C_contour_a(mt::MultiType; slides=false)
	jω_min = findmin(mt.L_mat[:,:,3,:])[2][1]

	C = zeros(length(mt.ct.agrid), length(mt.χgrid))
	for ja in 1:length(mt.ct.agrid), jχ in 1:length(mt.χgrid)
		C[ja, jχ] = mt.C_mat[jω_min, jχ, 3, ja]
	end
	return plot_L_contour(mt.ct.agrid, mt.χgrid, C; name_x = "a", name_y="C", slides=slides)
end

function plot_L_contour_a(mt::MultiType; slides=false)
	L = zeros(length(mt.ct.agrid), length(mt.χgrid))
	for ja in 1:length(mt.ct.agrid), jχ in 1:length(mt.χgrid)
		minL, jω = findmin(mt.L_mat[:, jχ, 3, ja])
		L[ja, jχ] = minL
	end
	return plot_L_contour(mt.ct.agrid, mt.χgrid, L; name_x = "a", name_y="𝓛", slides=slides)
end

function plot_C_contour(mt::MultiType; slides=false)
	ja_min = findmin(mt.L_mat[:,:,3,:])[2][3]

	C = zeros(length(mt.ωgrid), length(mt.χgrid))
	for jω in 1:length(mt.ωgrid), jχ in 1:length(mt.χgrid)
		C[jω, jχ] = mt.C_mat[jω, jχ, 3, ja_min]
	end
	return plot_L_contour(mt.ωgrid, mt.χgrid, C; name_y="C", slides=slides)
end

function plot_L_contour(mt::MultiType; slides=false)
	L = zeros(length(mt.ωgrid), length(mt.χgrid))
	for jω in 1:length(mt.ωgrid), jχ in 1:length(mt.χgrid)
		minL, ja = findmin(mt.L_mat[jω, jχ, 3, :])
		L[jω, jχ] = minL
	end
	return plot_L_contour(mt.ωgrid, mt.χgrid, L; name_y="𝓛", slides=slides)
end

function plot_L_contour(xgrid, ygrid, L_mat; name_x="ω", name_y="𝓛", slides::Bool=false)

	L_filled, temp = findmin(L_mat[.!isnan.(L_mat)])
	jjxy = findfirst(L_mat.==L_filled)

	# _, jjxy = findmin(L_mat)

	ygrid = annualized.(ygrid)
	if name_x == "ω"
		xgrid = perc_rate(xgrid)
		xtitle = "Decay rate  (<i>%</i>)"
	elseif name_x == "a"
		xgrid = annualized.(xgrid)
		xtitle = "Initial inflation  (<i>a<sub>0</sub></i>)"
	end

	
	xmin = xgrid[jjxy[1]]
	ymin = ygrid[jjxy[2]]

	if name_y == "𝓛"
		title = "lim<sub><i>p→0</i></sub> min<sub><i>a</i></sub> 𝓛(<i>p,a,ω,χ</i>)"
		shape_vec = [attr(;x0=xmin-0.001, x1 = xmin+0.001, y0 = ymin-0.002, y1=ymin+0.002, line_color="#08282e", type="circle")]
	elseif name_y == "C"
		title = "lim<sub><i>p→0</i></sub> C(<i>p,a*,ω,χ</i>)"
		shape_vec = []
	end

	colpal = ColorSchemes.lapaz

	ctχω = contour(;
		x = xgrid, y = ygrid,
		z = L_mat,
		# colorscale = vcat([[jj, get(colpal, jj)] for jj in range(0,1,length=50)][1:49]
		# 	# ,[[1, "fafafa"]]
		# 	, [[1, get(ColorSchemes.lajolla, 0.1)]]
		# 	), reversescale = true
		colorscale = [[jj, get(colpal, 1-jj)] for jj in range(0,1,length=50)]
		)
	p1 = plot(ctχω, Layout(;title=title, xaxis_title=xtitle, yaxis_title="Asymptote  (<i>χ</i>)", shapes = shape_vec))
	if slides
		relayout!(p1, font_family = "Lato", font_size = 16, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", width = 800, height = 450)
	else
		relayout!(p1, font_family="Linux Libertine", font_size = 17, width=900, height=450)
	end

	return p1
end

function plot_announcements(;slides::Bool=true, darkslides::Bool=false, exts::Vector=[], cond::Bool=false, add_opt::Bool=false, cond_t::Bool=false)
	xvec = 0:0.25:10

	darkslides ? slides = true : nothing

	cond_t ? cond = true : nothing

	slides ? colorpal = ColorSchemes.davos : colorpal = ColorSchemes.southwest
	slides ? correction = 0.8 : correction = 1
	darkslides ? fl = 0.2 : fl = 0

	line_opt = scatter(;x=xvec, y=((1.750133)-(0.784)) * exp.(-0.4951.*(4.0.*xvec)).+(0.784), showlegend=false, marker_color="#d62728", line_width=3, line_dash="dash")

	lines = [scatter(;x=xvec, y=(a0-χ) * exp.(-ω.*(xvec)).+χ, showlegend=false, marker_color=get(colorpal, fl + correction*χ/2)) for a0 in range(0,2, length=5) for ω in range(0,0.8,length=3) for (jχ, χ) in enumerate(range(2,0,length=5)) if ω != 0.0]

	plotname = "announcements"
	annotations = []
	if cond
		lines = [lines[43]]
		plotname *= "_cond"
		te = 9*4+1
		xe = lines[1][:x][te]
		ye = lines[1][:y][te]
		col_line = lines[1][:marker][:color]
		push!(annotations, attr(; x=xe, y=ye+0.05, text="<i>c", font_color=col_line, showarrow=false))
	end

	if add_opt
		push!(lines, line_opt)
		plotname *= "_w_opt"
	end

	darkslides ? shape_col = get(ColorSchemes.davos, 0.9) : shape_col = get(ColorSchemes.darkrainbow, 0.12)
	shapes = []
	if cond_t
		tt = 11
		x0 = lines[1][:x][tt]
		y0 = lines[1][:y][tt]
		shapes = [vline(x0, line_dash = "dash", line_color=shape_col); attr(;x0=x0-1*0.03, x1 = x0+1*0.03, y0 = y0-1*0.01, y1=y0+1*0.01, line_color=shape_col, fillcolor=shape_col, type="circle")]
		push!(annotations,attr(; x=x0 + 0.05, y=y0 + 0.01, text="<i>a<sub>t</sub><sup>c</sup>", ax=35, font_color = shape_col, font_size=24, font_family="Lato"))
		plotname *="_t"
	end

	if slides
		style = slides_def
		if darkslides
			style = slides_dark
		end
		plotname *= "_slides"
		width = slides_dark.layout.fields[:width] * 0.9
	else
		style = paper
		# font = "Linux Libertine"
		# bg = ""
		# plotname *= "_paper"
		# width = 400
	end

	p1 = plot(lines, style=style, Layout(;xaxis_title="<i>Quarters", yaxis_range=[-0.1;2.1], yaxis_title="%", title="Inflation announcements", shapes = shapes, annotations=annotations)
		)

	if darkslides
		relayout!(p1, font_color="white")
	end
	if length(exts) > 0
		for (jj, ext) in enumerate(exts)
			savefig(p1, pwd()*"/../Graphs/"*plotname*"."*ext)
		end
		return nothing
	end

	return p1
end


function plot_bayes(; center=1.5, dist=0.5, σ=0.5, p=0.25, distplot=4*sqrt(dist))

	a = center-dist
	g = center+dist

	ϵ_vec = range(center-distplot, center+distplot, length=1001)

	fa(x) = pdf.(Normal(0,σ), x .- a)
	fg(x) = pdf.(Normal(0,σ), x .- g)

	Bayes(p,x) = p .+ p.*(1.0.-p) .* (fa(x) .- fg(x)) ./ (p.*fa(x) .+ (1.0.-p).*fg(x))

	_, jj = findmin((fa(ϵ_vec) - fg(ϵ_vec)).^2)
	ϵstar = ϵ_vec[jj]

	shapes = [vline(ϵstar, fa(ϵstar), Bayes(p,ϵstar), line_dash="dashdot")]
	annotations = [attr(x=ϵstar, xanchor="left", y=(Bayes(p,ϵstar)+fa(ϵstar))/2, yanchor="center", text="π*", showarrow=false)]

	p1 = plot([
		scatter(;x=ϵ_vec, y=fa(ϵ_vec), name="<i>f(a-π)")
		scatter(;x=ϵ_vec, y=fg(ϵ_vec), name="<i>f(g-π)")
		scatter(;x=ϵ_vec, y=Bayes(p,ϵ_vec), name="<i>B(p,π,a,g)", line_width=2)
		],
		Layout(;yaxis_range=[-0.05;1.05], xaxis_range=[0,3], shapes=shapes, annotations=annotations)
		)

	relayout!(p1, font_family = "Lato", font_size = 16, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", width=1200, height=500)


	return p1
end

plot_plans_p(mt::MultiType; decay::Bool=true, slides::Bool=false) = plot_plans_p(mt.ct, mt.L_mat, mt.ωgrid, mt.χgrid, decay=decay, slides=slides)
function plot_plans_p(ct::CrazyType, L_mat, ωgrid, χgrid; decay::Bool=true, slides::Bool=false)

	ωvec = zeros(ct.Np)
	avec = zeros(ct.Np)
	χvec = zeros(ct.Np)

	data = zeros(ct.Np,3)
	for jp in 1:ct.Np
		_, jj = findmin(L_mat[:,:,jp,:])

		if decay
			data[jp, 1] = perc_rate(ωgrid[jj[1]])
		else
			data[jp, 1] = ωgrid[jj[1]]
		end
		data[jp, 2] = annualized.(ct.agrid[jj[3]])
		data[jp, 3] = annualized.(χgrid[jj[2]])
	end


	datanames = ["ω", "a<sub>0", "χ"]
	cols = [get(ColorSchemes.southwest, jj) for jj in [0, 0.5, 1]]
	ls = Vector{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, 0)

	yax = ["y2", "y1", "y1"]
	for jj in 1:3
		col = cols[jj]
		push!(ls, scatter(;x=ct.pgrid[3:end], y=data[3:end, jj], line_width = 2.5, yaxis=yax[jj], marker_color=col, name="<i>"*datanames[jj]*"</i>"))
	end

	layout = Layout(
		yaxis = attr(domain=[0, 0.45], zeroline=false, title="<i>%"),
		yaxis2 = attr(domain=[0.55, 1], zeroline=false, title="<i>%"),
		xaxis = attr(zeroline=false, title="<i>p<sub>0"),
		legend = attr(orientation="h", x=0.05),
		width = 900, height = 450,
		font_size=17, font_family="Linux Libertine", title="Preferred plans"
		)

	p1 = plot(ls, layout)

	# pω = plot(scatter(;x=ct.pgrid[3:end], y=ωvec[3:end], line_width=2.5, name="<i>ω", marker_color=get(ColorSchemes.southwest, 0.0)));
	# pχa= plot([
	# 	scatter(;x=ct.pgrid[3:end], y=annualized.(avec[3:end]), line_width=2.5, name="<i>a", marker_color=get(ColorSchemes.southwest, 0.5))
	# 	scatter(;x=ct.pgrid[3:end], y=annualized.(χvec[3:end]), line_width=2.5, name="<i>χ", marker_color=get(ColorSchemes.southwest, 0.99))
	# 	], Layout(;yaxis_title="%", xaxis_title="<i>p</i>"));

	# relayout!(pω,  xaxis_zeroline=false, yaxis_zeroline=false)
	# relayout!(pχa, xaxis_zeroline=false, yaxis_zeroline=false)

	# p1 = [pω; pχa]
	if slides
		relayout!(p1, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa")
		relayout!(p1, height=600, width=900, font_family="Lato")
	end

	return p1
end

function make_colorbar(ct::CrazyType; slides::Bool=true)
	agr = annualized.(ct.agrid)
	Na = length(ct.agrid)

	if slides
		layout = Layout(font_family="Lato", font_size=18, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", width=600, height=1000)
	else
		layout = Layout(font_family = "Linux Libertine", font_size=18)
	end

	p1 = plot(contour(;x=range(0,1,length=2), y=range(0,1,length=Na), z=[jy for jx in 1:2, jy in agr], colorscale=[[vv, get(ColorSchemes.lapaz, vv)] for vv in range(0,1,length=Na)], colorbar=attr(title="<i>a"), contours_coloring="heatmap")
		, layout)

	return p1

end

function plot_mimic_z(mt::MultiType, N=50; style::Style=slides_def, decay::Bool=true, CIs::Bool=false)

	data, datanames, zgrid = mimic_z(mt, N, decay=decay, annualize=true)

	cols = [get(ColorSchemes.southwest, jj) for jj in [0, 0.5, 1]]
	ls = Vector{PlotlyBase.GenericTrace{Dict{Symbol,Any}}}(undef, 0)

	yax = ["y2", "y1", "y1"]
	for jj in 1:3
		col = cols[jj]
		if CIs
			push!(ls, scatter(;x = zgrid, y = data[:,jj]+data[:,jj+3], yaxis=yax[jj], marker_color=col, mode="lines", opacity = 0.5, showlegend=false, line_width=0.01, hoverinfo="skip"))
			push!(ls, scatter(;x = zgrid, y = data[:,jj]-data[:,jj+3], yaxis=yax[jj], marker_color=col, mode="lines", opacity = 0.5, fill="tonexty", showlegend=false, line_width=0.01, hoverinfo="skip"))
		end
		push!(ls, scatter(;x=zgrid, y=data[:, jj], yaxis=yax[jj], marker_color=col, name="𝔼[<i>"*datanames[jj]*"</i>]"))
	end

	layout = Layout(
		yaxis = attr(domain=[0, 0.45], title="<i>%"),
		yaxis2 = attr(domain=[0.55, 1], title="<i>%"),
		xaxis = attr(title="<i>z"),
		legend = attr(orientation="h", x=0.05),
		title="Average plans",
		width = 700, height = 300, font_size = 16,
		)

	p1 = plot(ls, layout, style=style)

	return p1
end

function save_plot_mimic_z(mt::MultiType, N=50; slides::Bool=true, CIs::Bool=false)
	slides ? sty = slides_def : sty = paper
	p1 = plot_mimic_z(mt, N; style=sty, CIs=CIs)
	savejson(p1, "../Graphs/tests/mimics$(ifelse(CIs, "_CI", ""))$(ifelse(slides, "_slides", "_paper")).json")
	# if !CIs
	# 	p1, p2 = strategy_μ(mt, style=sty)
	# 	savejson(p1, "../Graphs/tests/marg_achi$(ifelse(slides, "_slides", "_paper")).json")
	# 	savejson(p2, "../Graphs/tests/marg_omegachi$(ifelse(slides, "_slides", "_paper")).json")		
	# end
	nothing
end

function strategy_μ(mt::MultiType; style, save_stats=false)

	c1, c2 = strategy_μ(mt, save_stats)
	
	l1 = Layout(
		title = "lim<sub>z→0</sub>∫<i>μ<sub>z</sub></i> (<i>ω, χ, a<sub>0</sub></i>) d<i>ω</i>",
		font_size=16,
		xaxis = attr(title="Initial inflation (<i>a<sub>0</sub></i>)"),
		yaxis = attr(title="Asymptote (<i>χ</i>)"),
	)
	l2 = Layout(
		title = "lim<sub>z→0</sub>∫<i>μ<sub>z</sub></i> (<i>ω, χ, a<sub>0</sub></i>) d<i>a<sub>0</sub></i>",
		xaxis = attr(title="Decay rate (<i>%</i>)"),
		yaxis = attr(title=""),
		)

	p1 = plot(c1, l1, style=style)
	p2 = plot(c2, l2, style=style)

	return p1, p2
end

function strategy_μ(mt::MultiType, save_stats = false)

	χgrid = mt.χgrid
	ωgrid = mt.ωgrid
	agrid = mt.ct.agrid

	marg_aχ = [sum(mt.μ[:, jχ, ja]) for jχ in 1:size(mt.μ,2), ja in 1:size(mt.μ,3)];

	marg_ωχ = [sum(mt.μ[jω, jχ, :]) for jω in 1:size(mt.μ, 1), jχ in 1:size(mt.μ,2)]

	# min_z1, max_z1 = extrema(marg_aχ)
	# min_z2, max_z2 = extrema(marg_ωχ)
	# min_z = min(min_z1, min_z2)
	# max_z = max(max_z1, max_z2)

	p1 = contour(y=annualized.(χgrid), x=annualized.(agrid), z=marg_aχ', colorscale=[[jj, get(ColorSchemes.lapaz, jj)] for jj in range(0,1,length=50)])

	p2 = contour(x=perc_rate.(ωgrid), y=annualized.(χgrid), z=marg_ωχ, colorscale=[[jj, get(ColorSchemes.lapaz, jj)] for jj in range(0,1,length=50)])
	# relayout!(p2, xaxis_title="Decay rate (<i>%</i>)", title="lim<sub>z→0</sub>∫<i>μ<sub>z</sub></i> (<i>ω, χ, a<sub>0</sub></i>) d<i>a<sub>0</sub></i>")


	P = sum([sum(mt.μ[:,jχ,ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if agrid[ja]>χgrid[jχ]])
	print("P(a_0 > χ) = $(@sprintf("%0.3g",100P))%")
	save_stats && write("../pa_chi.txt", "$(@sprintf("%0.3g",100P))\\%.")

	P = sum([sum(mt.μ[:,jχ,ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if agrid[ja]>5χgrid[jχ]])
	save_stats && write("../pa_chi5.txt", "$(@sprintf("%0.3g",100P))\\%.")

	P = sum([sum(mt.μ[:,jχ,ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if χgrid[jχ]==0])
	save_stats && write("../pa_chi0.txt", "$(@sprintf("%0.3g",100P))\\%.")

	P = sum([sum(mt.μ[jω,jχ,ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid), jω in 1:length(ωgrid) if perc_rate(ωgrid[jω])<=10])
	save_stats && write("../pa_omega0.txt", "$(@sprintf("%0.3g",100P))\\%.")

	return p1, p2
end

function comp_plot_planner(mt::MultiType; makeplots::Bool=false, slides::Bool=false)
	ct = mt.ct
	rp = Ramsey(ct)

	vfi!(rp)
	πR, θR = simul_plan(rp)

	# tvec = 1:length(πR)
	tvec = 1:11

	mean_ω, mean_a, mean_χ, sd_ω, sd_a, sd_χ = find_plan_μ(mt; annualize=false, decay=false)

	πC = (mean_a - mean_χ) * exp.(-mean_ω * (tvec.-1)) .+ mean_χ

	minL, jj = findmin(mt.L_mat[:, :, 3, :])
	ωK = mt.ωgrid[jj[1]]
	χK = mt.χgrid[jj[2]]
	aK = mt.ct.agrid[jj[3]]

	πK = (aK - χK) * exp.(-ωK * (tvec.-1)) .+ χK
	p1 = plot()

	iter_slides = [true, false]
	if slides
		iter_slides = [true]
	end
	for slides in iter_slides
		if slides
			ff = "Lato"
			bgcol = "#fafafa"
			wi = 800
		else
			ff = "Linux Libertine"
			bgcol = ""
			wi = 1000
		end

		layout = Layout(title="Plans", yaxis_title="%", xaxis_title="<i>Quarters", font_size=18, font_family = ff, width = wi, height=350, paper_bgcolor=bgcol, plot_bgcolor=bgcol, xaxis_zeroline=false, yaxis_zeroline=false, legend=attr(orientation="h", x=0.05))
		p1 = plot([
			scatter(x=tvec.-1, y=annualized.(πR)[tvec], marker_color=get(ColorSchemes.southwest, 0.01), name="<i>Ramsey")
			scatter(x=tvec.-1, y=annualized.(πC)[tvec], marker_color=get(ColorSchemes.southwest, 0.99), name="<i>Average eq'm")
			scatter(x=tvec.-1, y=annualized.(πK)[tvec], marker_color=get(ColorSchemes.southwest, 0.5), name="<i>Kambe eq'm")
			], layout)

		if makeplots
			savejson(p1, pwd()*"/../Graphs/tests/ramsey$(ifelse(slides, "_slides", "_paper")).json")
		end
	end
	return p1
end


function make_sustainable_plots(mt::MultiType, K; pc::DataType=Fwd_strategy, makeplots::Bool=false)

	ct = mt.ct
	rp = Ramsey(ct)

	vfi!(rp, verbose=false)
	πR, θR = simul_plan(rp)

	# tvec = 1:length(πR)
	tvec = 1:11

	if pc == Fwd_GP
		mult = range(0.0,1.75,length=K)
	elseif pc == Fwd_strategy
		mult = range(0.0,0.75,length=K)
	end

	π_sust = zeros(length(tvec), K)
	a_sust = zeros(length(tvec), K)
	show_vec = Vector{Bool}(undef, K)

	if pc == Fwd_strategy
		plotname = "sustainable"
	elseif pc == Fwd_GP
		plotname = "GP"
	end

	sp = Sustainable(ct, pc = pc, ξ = 0.0);
	vfi!(sp)
	old_g = sp.g
	old_v = sp.v
	for (jj, jv) in enumerate(mult)
		
		sp = Sustainable(ct, ξ = jv*Nash(ct), pc = pc)
		sp.g = old_g
		sp.v = old_v
		println("$jj")
		flag = vfi!(sp, verbose=true, maxiter = 750)
		if !flag
			flag2 = false
			subiter = 0
			while flag2 == false && subiter < 2
				subiter += 1
				flag2 = vfi!(sp, verbose=true, maxiter = 750)
			end
			if flag2
				flag = true
			end
		end
		old_g = sp.g
		old_v = sp.v

		πv, θv = simul_plan(sp)

		π_sust[:, jj] = πv[tvec]
		a_sust[:, jj] = θv[tvec]
		show_vec[jj] = flag
	end
	p1 = plot()
	p2 = plot()
	for slides in [true, false]
		if slides
			ff = "Lato"
			bgcol = "#fafafa"
			wi = 800
		else
			ff = "Linux Libertine"
			bgcol = ""
			wi = 1000
		end

		layout = Layout(title="Plans", yaxis_title="%", xaxis_title="<i>Quarters", font_size=18, font_family = ff, width = wi, height=350, paper_bgcolor=bgcol, plot_bgcolor=bgcol, xaxis_zeroline=false, yaxis_zeroline=false, legend=attr(orientation="h", x=0.05))

		xs = [mean(tvec.-1) for jj in 1:K if show_vec[jj]]
		ys = [mean(annualized.(π_sust[tvec, jj])) for jj in 1:K if show_vec[jj]]
		cols = range(0,1,length=length(xs))
		colscale = [[vv, get(ColorSchemes.davos, 0.8*vv)] for vv in range(0,1,length=100)]
		colnames = round.(range(extrema(mult)..., length=6), digits = 2)

		p1 = plot([
			scatter(mode = "markers", x=xs, y=ys, showlegend = false, 
				marker = attr(color=cols, opacity = 0, colorscale=colscale, colorbar = attr(tickvals=range(0,1,length=length(colnames)), title="&nbsp;<i>ξ/π<sup>N", ticktext=colnames))
				)
			[scatter(x=tvec.-1, y=annualized.(π_sust[tvec, jj]), mode="lines", opacity=0.9, line_width=5, marker_color=get(ColorSchemes.davos, 0.8(jj)/K), name="$(mult[jj])", showlegend=false) for jj in 1:K if show_vec[jj]]
			scatter(x=tvec.-1, y=annualized.(πR[tvec]), line_dash="dash", line_width = 2, marker_color=get(ColorSchemes.lajolla, 0.6), name="<i>Ramsey")
			], layout)
		p2 = plot([
			[scatter(x=tvec.-1, y=annualized.(π_sust[tvec, jj]), mode="lines", opacity=0.9, line_width=5, line_dash="dash", marker_color=get(ColorSchemes.davos, 0.8(jj)/K), name="$(mult[jj])", showlegend=false) for jj in 1:K if show_vec[jj]]
			[scatter(x=tvec.-1, y=annualized.(a_sust[tvec, jj]), mode="lines", opacity=0.9, line_width=3, marker_color=get(ColorSchemes.lajolla, 1-0.8(jj)/K), name="$(mult[jj])", showlegend=false) for jj in 1:K if show_vec[jj]]
			# scatter(x=tvec.-1, y=annualized.(πR[tvec]), line_dash="dash", marker_color=get(ColorSchemes.lajolla, 0.6), name="<i>Ramsey")
			], layout)

		if makeplots
			savejson(p1, pwd()*"/../Graphs/tests/$(plotname)_$(ifelse(slides, "slides", "paper")).json")
		end
	end

	return p1, π_sust, a_sust, p2

end

function plots_recursive(dk::DovisKirpalani, ct=CrazyType(dk); makeplots::Bool=false)

	rp = Ramsey(ct)
	vR = vfi!(rp)
	πR, θR = simul_plan(rp)

	tvec = 1:11
	T = length(tvec)
	K = length(dk.pgrid)
	amat = zeros(T, K)

	for jp in 1:K
		p_vec, a_vec, π_vec, y_vec, g_vec, L_vec = simul(dk; jp0=jp, T=T, noshocks=true)
		amat[:,jp] = annualized.(a_vec)
	end

	p1 = plot()
	for slides in [true,false]
		if slides
			ff = "Lato"
			bgcol = "#fafafa"
			wi = 800
		else
			ff = "Linux Libertine"
			bgcol = ""
			wi = 1000
		end

		layout = Layout(title="Plans", yaxis_title="%", xaxis_title="<i>Quarters", font_size=18, font_family = ff, width = wi, height=350, paper_bgcolor=bgcol, plot_bgcolor=bgcol, xaxis_zeroline=false, yaxis_zeroline=false, legend=attr(orientation="h", x=0.05))

		xs = [mean(tvec.-1) for (jp, pv) in enumerate(dk.pgrid) if jp >= 3]
		ys = [mean(amat[:, jp]) for (jp, pv) in enumerate(dk.pgrid) if jp >= 3]
		cols = range(0,1,length=length(xs))
		colscale = [[vv, get(ColorSchemes.davos, 0.8*vv)] for vv in range(0,1,length=100)]
		colnames = round.(range(extrema(dk.pgrid)..., length=6), digits = 2)


		p1 = plot([
			scatter(x=xs, y=ys, mode="markers", showlegend=false,
				marker = attr(color=cols, opacity = 0, colorscale=colscale, colorbar = attr(tickvals=range(0,1,length=length(colnames)), title="&nbsp;&nbsp;<i>p<sub>0", ticktext=colnames))
				)
			[scatter(x=tvec.-1, y=amat[:, jp], mode="lines+markers", opacity=0.9, line_width=2, marker_color=get(ColorSchemes.davos, 0.05+0.8*jp/K), name="<i>p₀ = $(pv)", showlegend=false) for (jp, pv) in enumerate(dk.pgrid) if jp >= 3]
			scatter(x=tvec.-1, y=annualized.(πR[tvec]), line_dash="dash", line_width = 3, marker_color=get(ColorSchemes.lajolla, 0.6), name="<i>Ramsey")
			], layout)

		# p2 = plot([
		# 	scatter(x=dk.pgrid[2:end], y = [amat[1,jp] for jp in 2:length(dk.pgrid)], name="<i>a<sub>0")
		# 	scatter(x=dk.pgrid[2:end], y = [amat[end,jp] for jp in 2:length(dk.pgrid)], name="<i>χ")
		# 	], Layout(yaxis_title="%", xaxis_title="<i>Quarters", font_size=18, font_family = ff, width = wi, height=350, paper_bgcolor=bgcol, plot_bgcolor=bgcol, xaxis_zeroline=false, yaxis_zeroline=false, legend=attr(orientation="h", x=0.05)))

		# p1 = [p1; p2]

		if makeplots
			savejson(p1, pwd()*"/../Graphs/tests/all_recursive_$(T)_$(ifelse(slides, "slides", "paper")).json")
		end
	end
	return p1
end




function comp_plot_recursive(dk::DovisKirpalani{T}, mt::MultiType; makeplots::Bool=false) where T<:PhillipsCurve

	ct = mt.ct
	rp = Ramsey(ct)

	vR = vfi!(rp)
	πR, θR = simul_plan(rp)

	# tvec = 1:length(πR)
	tvec = 1:11

	mean_ω, mean_a, mean_χ, sd_ω, sd_a, sd_χ = find_plan_μ(mt; annualize=false, decay=false)

	πC = (mean_a - mean_χ) * exp.(-mean_ω * (tvec.-1)) .+ mean_χ

	minL, jj = findmin(mt.L_mat[:, :, 3, :])
	ωK = mt.ωgrid[jj[1]]
	χK = mt.χgrid[jj[2]]
	aK = mt.ct.agrid[jj[3]]

	πK = (aK - χK) * exp.(-ωK * (tvec.-1)) .+ χK
	vK = findmin(ct.L[3,:])[1]

    p_vec, a_vec, π_vec, y_vec, g_vec, L_vec = simul(dk; jp0=3, T=4*length(tvec), noshocks=true)


	p1 = plot()

	print("\nRamsey initial: $(round(annualized(πR[1]),digits=3))")
	print("\nRamsey final: $(round(annualized(πR[end]),digits=3))")
	print("\nRamsey value: $(round(show_value(rp),digits=4))")
	print("\nKambe initial: $(round(annualized(πK[1]),digits=3))")
	print("\nKambe final: $(round(annualized(χK),digits=3))")
	print("\nKambe value: $(round(vK,digits=4))")
	print("\nRecursive initial: $(round(annualized(a_vec[1]),digits=3))")
	print("\nRecursive final: $(round(annualized(a_vec[end]),digits=3))")
	print("\nRecursive value: $(round(L_vec[1],digits=4))")

	for slides in [true, false]
		if slides
			ff = "Lato"
			bgcol = "#fafafa"
			wi = 800
		else
			ff = "Linux Libertine"
			bgcol = ""
			wi = 1000
		end

		layout = Layout(title="Plans", yaxis_title="%", xaxis_title="<i>Quarters", font_size=18, font_family = ff, width = wi, height=350, paper_bgcolor=bgcol, plot_bgcolor=bgcol, xaxis_zeroline=false, yaxis_zeroline=false, legend=attr(orientation="h", x=0.05))
		p1 = plot([
			scatter(x=tvec.-1, y=annualized.(πR[tvec]), marker_color=get(ColorSchemes.southwest, 0.01), name="<i>Ramsey")
			scatter(x=tvec.-1, y=annualized.(πC)[tvec], marker_color=get(ColorSchemes.southwest, 0.99), name="<i>Average eq'm")
			scatter(x=tvec.-1, y=annualized.(πK)[tvec], marker_color=get(ColorSchemes.southwest, 0.5), name="<i>Kambe eq'm")
			scatter(x=tvec.-1, y=annualized.(a_vec)[tvec], marker_color=get(ColorSchemes.lajolla, 0.6), name="<i>Recursive")
			], layout)

		if makeplots
			savejson(p1, pwd()*"/../Graphs/tests/recursive_$(T)_$(ifelse(slides, "slides", "paper")).json")
		end
	end
	return p1
end

dist_plans(a1::Vector, a2::Vector, pp::Plan) = dist_plans(a1,a2,pp.β)
function dist_plans(a1::Vector, a2::Vector, β::Float64)
	if length(a1) == length(a2)
		T = length(a1)
		return sqrt(sum([β^tt * (a1[tt]-a2[tt])^2 for tt in 1:T]) / T )
	else
		throw(error("Both plans should have the same length"))
	end
end

function eval_close_ct(dk, ω, χ, a0, avec)
	T = length(avec)
	tvec = 1:T
	πvec = (a0 - χ) * exp.(-ω * (tvec.-1)) .+ χ
	return dist_plans(annualized.(πvec), annualized.(avec), dk)
end

function find_closest_ct(dk::DovisKirpalani, mt::MultiType; cont::Bool=false)
	p_vec, a_vec, π_vec, y_vec, g_vec, L_vec = simul(dk; jp0=3, T=40, noshocks=true)

	guess = [0.5; Nash(dk); Nash(dk)]

	if cont
		res = Optim.optimize(
			x -> eval_close_ct(dk, x[1], x[2], x[3], a_vec), guess
			)
		ωP, χP, aP = res.minimizer
	else
		χgrid, ωgrid = mt.χgrid, mt.ωgrid
		agrid = mt.ct.agrid
		curr_min = 1e8

		jωP, jχP, jaP = 0, 0, 0
		for (jω, ωv) in enumerate(ωgrid), (jχ, χv) in enumerate(χgrid), (ja, av) in enumerate(agrid)
			vvv = eval_close_ct(dk, ωv, χv, av, a_vec)

			if vvv < curr_min
				ωP, χP, aP = ωv, χv, av
				jωP, jχP, jaP = jω, jχ, ja
				curr_min = vvv
			end
		end
		Lct = mt.L_mat[jωP, jχP, 3, jaP]
		print("\nRecursive value: $(round(L_vec[1],digits=4))")
		print("\nProjection value: $(round(Lct,digits=4))")
	end


	tvec = 1:11
	πP = (aP - χP) * exp.(-ωP * (tvec.-1)) .+ χP
	plot([
		scatter(x=tvec.-1, y = annualized.(a_vec[tvec]), name="Recursive")
		scatter(x=tvec.-1, y = annualized.(πP), name="Projection")
		])
end
