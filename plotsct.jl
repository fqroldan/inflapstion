using PlotlyJS, ColorSchemes

""" Define styles """

sand() = "#F5F3F1"
darkbgd() = "#272929"
lightgrid() = "#353535"
darkgrid() = "#e2e2e2"
gridcol(dark=false) = ifelse(dark, lightgrid(), darkgrid())

q_axis(dark) = attr(showgrid=true, gridcolor=gridcol(dark), gridwidth=0.5, zeroline=false)
bgcol(slides, dark) = ifelse(slides, ifelse(dark, darkbgd(), sand()), "white")
qleg() = attr(orientation="h", x=0.05, xanchor="left")

qwidth(slides) = 864
qheight(slides) = ceil(Int, qwidth(slides) * ifelse(slides, 10 / 16, 7 / 16))

function qtemplate(; dark=false, slides=!dark)
	axis = q_axis(dark)
	width = 864 #1920 * 0.45
	l = Layout(
		xaxis=axis, yaxis=axis,
		width=width,
		height=width * ifelse(slides, 10 / 16, 7 / 16),
		font=attr(
			family=ifelse(slides, "Lato", "Linux Libertine"),
			size=16, color=ifelse(dark, sand(), darkbgd())
		),
		paper_bgcolor=bgcol(slides, dark), plot_bgcolor=bgcol(slides, dark),
		legend=qleg(),
	)
	return Template(layout=l)
end

function hplot(ct::CrazyType; kwargs...)
    y = [ct.gπ[jp, ja] - av for jp in eachindex(ct.pgrid), (ja, av) in enumerate(ct.agrid)]

    ctplot(ct, annualized.(y), title = "<i>g<sup>⋆</sup> - a"; kwargs...)
end


function Eplot(ct::CrazyType; kwargs...)
    y = [ct.Ep[jp, ja] .- pv for (jp, pv) in enumerate(ct.pgrid), ja in eachindex(ct.agrid)]

    ctplot(ct, y, title = "𝔼[<i>p'-p</i>]"; kwargs...)
end

Cplot(cc::CrazyType; kwargs...) = ctplot(ct, ct.C; kwargs...)
gplot(ct::CrazyType; kwargs...) = ctplot(ct, annualized.(ct.gπ); kwargs...)
Lplot(ct::CrazyType; kwargs...) = ctplot(ct, ct.L, title = "𝓛"; kwargs...)
function ctplot(ct::CrazyType, y::Array; slides=true, dark=false, mod_a = 1, kwargs...)

    vv = [ja / length(ct.agrid) for ja in eachindex(ct.agrid)]
    if dark
        col = get(ColorSchemes.lapaz, 1 .- 0.9 * vv, :clamp)
    else
        col = get(ColorSchemes.lapaz, 0.9 * vv, :clamp)
    end

    min_a, max_a = annualized.(extrema(ct.agrid)) ./ annualized.(Nash(ct))

    jp = round(Int, length(ct.pgrid)/2)
	xs = [ct.pgrid[jp] for _ in axes(y, 2)]
    ys = [y[jp, ja] for ja in axes(y, 2)]
    cols = range(0,1,length=length(xs))
	colscale = [[(jp-1)/(length(col)-1), col[jp]] for jp in eachindex(col)]
    colnames = round.(range(min_a, max_a, length=6), digits=2)

    data = [
        scatter(mode = "markers", marker_opacity = 0,
				x = xs, y = ys, showlegend=false,
				marker = attr(color=cols, reversescale=false, colorscale=colscale, colorbar = attr(tickvals=range(min_a, max_a,length=length(colnames)), title="&nbsp;<i>a/π<sup>N", ticktext=colnames))
				)
        [scatter(x = ct.pgrid, y = y[:, ja], line_width = 2, marker_color = col[ja], name = "a = $(round(annualized(av), sigdigits=2))") for (ja, av) in enumerate(ct.agrid) if ja % mod_a == 0]
    ]

    template = qtemplate(dark = dark, slides = slides)

    layout = Layout(template = template, xaxis_title = "<i>p", hovermode ="x unified", showlegend=false; kwargs...)

    plot(data[end:-1:1], layout)
end

function Cplot(mt::MultiType; jp = 2, kwargs...)

    Nω, Nχ, Np, Na = size(mt.L_mat)
    C = zeros(Na, Nχ)
    for ja in axes(C, 1), jχ in axes(C, 2)
        _, jω = findmin(mt.L_mat[:, jχ, jp, ja])
        C[ja, jχ] = mt.C_mat[jω, jχ, jp, ja]
    end

    ctωplot(mt, C, title="lim<sub><i>p→0</i></sub> 𝓛(<i>p,a,ω*,χ</i>)"; kwargs...)
end


function Lplot(mt::MultiType; jp = 2, kwargs...)
    L = [minimum(mt.L_mat[jω, jχ, jp, :]) for jω in axes(mt.L_mat,1), jχ in axes(mt.L_mat, 2)]

    jj = findfirst(L .== minimum(L))
	xmin = perc_rate(mt.ωgrid[jj[1]])
	ymin = annualized(mt.χgrid[jj[2]])

    title = "lim<sub><i>p→0</i></sub> min<sub><i>a</i></sub> 𝓛(<i>p,a,ω,χ</i>)"
    shape_vec = [attr(;x0=xmin-0.001, x1 = xmin+0.001, y0 = ymin-0.002, y1=ymin+0.002, line_color = "#08282e", fillcolor="#08282e", type = "circle")]

    ctplot(mt, L; title = title, shapes = shape_vec, kwargs...)
end

function Lωplot(mt::MultiType; jp=2, slides=true, dark=false, kwargs...)
    L = [minimum(mt.L_mat[:, jχ, jp, ja]) for ja in axes(mt.L_mat, 4), jχ in axes(mt.L_mat, 2)]

    jj = findfirst(L .== minimum(L))
    xmin = annualized(mt.ct.agrid[jj[1]])
    ymin = annualized(mt.χgrid[jj[2]])

    shape_vec = [
        attr(; x0=xmin - 0.001, x1=xmin + 0.001, y0=ymin - 0.002, y1=ymin + 0.002, line_color="darkred", fillcolor="darkred", type="circle")
        attr(; x0=0, x1 = annualized(mt.χgrid[end]), y0=0, y1 = annualized(mt.χgrid[end]), type = "line", line_width = 1, line_dash="dash", line_color = "darkred")
    ]
    
    ctωplot(mt, L, shapes = shape_vec, 
        title = "lim<sub><i>p→0</i></sub> min<sub><i>ω</i></sub> 𝓛(<i>p,a,ω,χ</i>)", 
        slides=slides, dark=dark; kwargs...)
end

function ctωplot(mt::MultiType, y::Array; title = "", slides = true, dark = false, kwargs...)
    Na, Nχ = length(mt.ct.agrid), length(mt.χgrid)
    @assert size(y) == (Na, Nχ)

    colpal = ColorSchemes.lapaz
    xt = "Initial Inflation (<i>a<sub>0</sub></i>)"
    yt = "Asymptote (χ)"

    data = contour(
        z = y', x = annualized.(mt.ct.agrid), y = annualized.(mt.χgrid),
        colorscale = [[jj, get(colpal, 1-jj, :clamp)] for jj in range(0,1,length=50)]

    )
    layout = Layout(
        xaxis_title = xt, yaxis_title = yt, title = title,
        template = qtemplate(slides=slides, dark=dark);
        kwargs...
    )
    
    plot(data, layout)
end


function ctplot(mt::MultiType, y::Array; slides = true, dark = false, kwargs...)
    
	colpal = ColorSchemes.lapaz
    xt = "Decay rate (%)"
    yt = "Asymptote (χ)"

    data = contour(
        z = y', x = perc_rate.(mt.ωgrid), y = annualized.(mt.χgrid),
        colorscale = [[jj, get(colpal, 1-jj, :clamp)] for jj in range(0,1,length=50)]

    )
    layout = Layout(
        xaxis_title = xt, yaxis_title = yt, 
        template = qtemplate(slides=slides, dark=dark);
        kwargs...
    )
    
    plot(data, layout)
end

function plansp(mt::MultiType; slides = true, dark = false)
    Np = mt.ct.Np

	data = zeros(Np, 3)
    for jp in axes(data, 1)
        
        _, jj = findmin(mt.L_mat[:, :, jp, :])
        
        data[jp, 1] = perc_rate(mt.ωgrid[jj[1]])
        data[jp, 2] = annualized(mt.ct.agrid[jj[3]])
        data[jp, 3] = annualized(mt.χgrid[jj[2]])
    end
    
    datanames = ["ω", "a<sub>0", "χ"]
    yax = ["y2", "y1", "y1"]
    
    cols = [get(ColorSchemes.southwest, jj, :clamp) for jj in [0, 0.5, 1]]
    lines = [
        scatter(;x=mt.ct.pgrid[2:end], y=data[2:end, jj], line_width = 2.5, yaxis=yax[jj], marker_color=cols[jj], name="<i>"*datanames[jj]*"</i>") for jj in eachindex(datanames)
    ]

	layout = Layout(
        template = qtemplate(dark = dark, slides = slides),
        hovermode = "x",
		yaxis = attr(domain=[0, 0.45], zeroline=false, title="<i>%"),
		yaxis2 = attr(domain=[0.55, 1], zeroline=false, title="<i>%"),
		xaxis = attr(zeroline=false, title="<i>p<sub>0"),
		legend = attr(orientation="h", x=0.05),
		title="Preferred plans",
    )

    plot(lines, layout)
end

function avgplans(mt::MultiType, N = 50; decay=true, CIs=false, slides = true, dark = false)
    data, datanames, zgrid = mimic_z(mt, N, decay=decay, annualize=true)

    cols = [get(ColorSchemes.southwest, jj, :clamp) for jj in (0, 0.5, 1)]
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
        template = qtemplate(slides=slides, dark=dark),
		yaxis = attr(domain=[0, 0.45], title="<i>%"),
		yaxis2 = attr(domain=[0.55, 1], title="<i>%"),
		xaxis = attr(title="<i>z"),
		legend = attr(orientation="h", x=0.05),
		title="Average plans",
		)

	plot(ls, layout)
end

function strategy_μ(mt::MultiType; save_stats=false, slides = true, dark = false)

    χgrid = mt.χgrid
    ωgrid = mt.ωgrid
    agrid = mt.ct.agrid

    marg_aχ = [sum(mt.μ[:, jχ, ja]) for jχ in axes(mt.μ, 2), ja in axes(mt.μ, 3)]

    marg_ωχ = [sum(mt.μ[jω, jχ, :]) for jω in axes(mt.μ, 1), jχ in axes(mt.μ, 2)]

    c1 = contour(x=annualized.(agrid), y=annualized.(χgrid), z=marg_aχ, colorscale=[[jj, get(ColorSchemes.lapaz, jj, :clamp)] for jj in range(0, 1, length=50)])

    c2 = contour(x=perc_rate.(ωgrid), y=annualized.(χgrid), z=marg_ωχ', colorscale=[[jj, get(ColorSchemes.lapaz, jj, :clamp)] for jj in range(0, 1, length=50)])

    P = sum([sum(mt.μ[:, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if agrid[ja] > χgrid[jχ]])
    print("P(a_0 > χ) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write("pa_chi.txt", "$(@sprintf("%0.3g",100P))\\%.")

    P = sum([sum(mt.μ[:, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if agrid[ja] > 5χgrid[jχ]])
    print("P(a_0 > 5χ) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write("pa_chi5.txt", "$(@sprintf("%0.3g",100P))\\%.")

    P = sum([sum(mt.μ[:, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if χgrid[jχ] == 0])
    print("P(a_0 > 0) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write("pa_chi0.txt", "$(@sprintf("%0.3g",100P))\\%.")

    P = sum([sum(mt.μ[jω, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid), jω in 1:length(ωgrid) if perc_rate(ωgrid[jω]) <= 10])
    print("P(ω ≤ 10%) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write("pa_omega0.txt", "$(@sprintf("%0.3g",100P))\\%.")

    l1 = Layout(
        template = qtemplate(slides=slides, dark=dark),
        title="lim<sub>z→0</sub>∫<i>μ<sub>z</sub></i> (<i>ω, χ, a<sub>0</sub></i>) d<i>ω</i>",
        font_size=16,
        xaxis=attr(title="Initial inflation (<i>a<sub>0</sub></i>)"),
        yaxis=attr(title="Asymptote (<i>χ</i>)"),
        shapes = [attr(x0=0, x1 = annualized(χgrid[end]), y0=0, y1=annualized(χgrid[end]), type = "line", line_color="darkred", line_dash="dash", line_width = 1)]
    )
    l2 = Layout(
        template = qtemplate(slides=slides, dark=dark),
        title="lim<sub>z→0</sub>∫<i>μ<sub>z</sub></i> (<i>ω, χ, a<sub>0</sub></i>) d<i>a<sub>0</sub></i>",
        xaxis=attr(title="Decay rate (<i>%</i>)"),
        yaxis=attr(title=""),
    )

    p1 = plot(c1, l1)
    p2 = plot(c2, l2)

    return p1, p2
end