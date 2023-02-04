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

gplot(ct::CrazyType; kwargs...) = ctplot(ct, annualized.(ct.gπ); kwargs...)
Lplot(ct::CrazyType; kwargs...) = ctplot(ct, ct.L, title = "𝓛"; kwargs...)
function ctplot(ct::CrazyType, y::Array; slides=true, dark=false, kwargs...)

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
        [scatter(x = ct.pgrid, y = y[:, ja], line_width = 2, marker_color = col[ja], name = "a = $(round(annualized(av), sigdigits=2))") for (ja, av) in enumerate(ct.agrid) if ja % 1 == 0]
    ]

    template = qtemplate(dark = dark, slides = slides)

    layout = Layout(template = template, xaxis_title = "<i>p", hovermode ="x unified", showlegend=false; kwargs...)

    plot(data[end:-1:1], layout)
end

function Lplot(mt::MultiType; jp = 2, kwargs...)
    L = [minimum(mt.L_mat[jω, jχ, jp, :]) for jω in axes(mt.L_mat,1), jχ in axes(mt.L_mat, 2)]

    jj = findfirst(L .== minimum(L))
	xmin = perc_rate(mt.ωgrid[jj[1]])
	ymin = annualized(mt.χgrid[jj[2]])

    title = "lim<sub><i>p→0</i></sub> min<sub><i>a</i></sub> 𝓛(<i>p,a,ω,χ</i>)"
    shape_vec = [attr(;x0=xmin-0.001, x1 = xmin+0.001, y0 = ymin-0.002, y1=ymin+0.002, line_color = "#08282e", fillcolor="#08282e", type = "circle")]

    xt = "Decay rate (%)"
    yt = "Asymptote (χ)"

    ctplot(mt, L; title = title, xaxis_title = xt, yaxis_title = yt, shapes = shape_vec, kwargs...)
end

function ctplot(mt::MultiType, y::Array; slides = true, dark = false, kwargs...)
    
	colpal = ColorSchemes.lapaz

    data = contour(
        z = y', x = perc_rate.(mt.ωgrid), y = annualized.(mt.χgrid),
        colorscale = [[jj, get(colpal, 1-jj, :clamp)] for jj in range(0,1,length=50)]

    )
    layout = Layout(
        template = qtemplate(slides=slides, dark=dark);
        kwargs...
    )
    
    plot(data, layout)
end

