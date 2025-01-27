using PlotlyJS, ColorSchemes

""" Define styles """

OwlBlue() = "#0098E9"
OwlRed() = "#FF5CA8"
OwlYellow() = "#F29318"
OwlOrange() = "#F97860"
OwlGreen() = "#5AA800"
OwlPurple() = "#807AC9"

defblue() = "#1f77b4"

Owls() = [OwlBlue(), OwlOrange(), OwlGreen(), OwlPurple(), OwlYellow(), OwlRed()]
Owls(k::Int64) = cycle(k - 1, Owls())

BlRd() = [
    "#001621",
    "#004164",
    "#006DA6",
    "#0098E9",
    "#3FB1ED",
    "#7DC9F2",
    "#BCE2F6",
    "#F4CEC3",
    "#EB8B71",
    "#E2481E"
]
BlRd(k::Int64) = cycle(k, BlRd())

cycle(k::Int64, v::Vector) = v[k%length(v)+1]

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


function plot_announcements(; slides = true, dark = false, add_opt=false, cond_t=false, cond=cond_t)
	xvec = 0:0.25:10

	slides ? colorpal = ColorSchemes.davos : colorpal = ColorSchemes.southwest
	slides ? correction = 0.8 : correction = 1
	dark ? fl = 0.2 : fl = 0
    
    colorpal = ColorSchemes.lapaz
    correction = 0.9

	line_opt = scatter(;x=xvec, y=((1.750133)-(0.784)) * exp.(-0.4951.*(4.0.*xvec)).+(0.784), showlegend=false, marker_color="#d62728", line_width=3, line_dash="dash")

	lines = [scatter(;x=xvec, y=(a0-χ) * exp.(-ω.*(xvec)).+χ, showlegend=false, marker_color=get(colorpal, fl + correction*χ/2, :clamp)) for a0 in range(0,2, length=5) for ω in range(0,0.8,length=3) for (jχ, χ) in enumerate(range(2,0,length=5)) if ω != 0.0]

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

	dark ? shape_col = get(ColorSchemes.davos, 0.9, :clamp) : shape_col = get(ColorSchemes.darkrainbow, 0.12, :clamp)
	shapes = []
	if cond_t
		tt = 11
		x0 = lines[1][:x][tt]
		y0 = lines[1][:y][tt]
		shapes = [vline(x0, line_dash = "dash", line_color=shape_col); attr(;x0=x0-1*0.03, x1 = x0+1*0.03, y0 = y0-1*0.01, y1=y0+1*0.01, line_color=shape_col, fillcolor=shape_col, type="circle")]
		push!(annotations,attr(; x=x0 + 0.05, y=y0 + 0.01, text="<i>a<sub>t</sub><sup>c</sup>", ax=35, font_color = shape_col, font_size=24, font_family="Lato"))
		plotname *="_t"
	end

    layout = Layout(
        template = qtemplate(slides=slides, dark=dark),
        xaxis_title="<i>Quarters", yaxis_range=[-0.1;2.1], yaxis_title="%",
        title="Inflation announcements", shapes = shapes, 
        annotations=annotations,
        )

	plot(lines, layout)
end


function hplot(ct::CrazyType; share = false, kwargs...)
    y = [ct.gπ[jp, ja] - av for jp in eachindex(ct.gr[:p]), (ja, av) in enumerate(ct.gr[:a])]
    y = annualized.(y)

    title = "<i>g<sup>⋆</sup> - a"
    yaxis_title = "%"
    if share
        y *= 1/annualized(Nash(ct))
        yaxis_title = "Share of Nash inflation"
    end

    ctplot(ct, y; title, yaxis_title, kwargs...)
end


function Eplot(ct::CrazyType; kwargs...)
    y = [ct.Ep[jp, ja] .- pv for (jp, pv) in enumerate(ct.gr[:p]), ja in eachindex(ct.gr[:a])]

    ctplot(ct, y, title = "𝔼[<i>p'-p</i>]"; kwargs...)
end

Cplot(ct::CrazyType; kwargs...) = ctplot(ct, ct.C; kwargs...)
gplot(ct::CrazyType; kwargs...) = ctplot(ct, annualized.(ct.gπ); kwargs...)
Lplot(ct::CrazyType; kwargs...) = ctplot(ct, ct.L, title = "𝓛"; kwargs...)
function ctplot(ct::CrazyType, y::Array; slides=true, dark=false, mod_a = 1, kwargs...)

    vv = [ja / length(ct.gr[:a]) for ja in eachindex(ct.gr[:a])]
    if dark
        col = get(ColorSchemes.lapaz, 1 .- 0.9 * vv, :clamp)
    else
        col = get(ColorSchemes.lapaz, 0.95 * vv, :clamp)
    end

    min_a, max_a = annualized.(extrema(ct.gr[:a])) ./ annualized.(Nash(ct))

    jp = round(Int, length(ct.gr[:p])/2)
	xs = [ct.gr[:p][jp] for _ in axes(y, 2)]
    ys = [y[jp, ja] for ja in axes(y, 2)]
    cols = range(0,1,length=length(xs))
	colscale = [[(jp-1)/(length(col)-1), col[jp]] for jp in eachindex(col)]
    colnames = round.(range(min_a, max_a, length=6), digits=2)

    data = [
        scatter(mode = "markers", marker_opacity = 0,
				x = xs, y = ys, showlegend=false,
				marker = attr(color=cols, reversescale=false, colorscale=colscale, colorbar = attr(tickvals=range(min_a, max_a,length=length(colnames)), title="&nbsp;<i>a/π<sup>N", ticktext=colnames))
				)
        [scatter(x = ct.gr[:p], y = y[:, ja], line_width = 1.75, marker_color = col[ja], name = "a = $(round(annualized(av), sigdigits=2))") for (ja, av) in enumerate(ct.gr[:a]) if ja % mod_a == 0]
    ]

    template = qtemplate(dark = dark, slides = slides)

    layout = Layout(template = template, xaxis_title = "<i>p", hovermode ="x", showlegend=false; kwargs...)

    plot(data[end:-1:1], layout)
end

function astar(ct::CrazyType; dark = false, slides = true, kwargs...)

    avec = [annualized(ct.gr[:a][findmin(ct.L[jp,:])[2]]) for jp in eachindex(ct.gr[:p])]

    sc = scatter(x = ct.gr[:p], y = avec)

    template = qtemplate(; dark, slides)

    layout = Layout(; template, xaxis_title = "<i>p", yaxis_title = "<i>%", hovermode = "x", kwargs...)

    plot(sc, layout)
end



function Cplot(mt::MultiType; jp = 2, kwargs...)

    Nω, Nχ, Np, Na = size(mt.L_mat)
    C = zeros(Na, Nχ)
    for ja in axes(C, 1), jχ in axes(C, 2)
        _, jω = findmin(mt.L_mat[:, jχ, jp, ja])
        C[ja, jχ] = mt.C_mat[jω, jχ, jp, ja]
    end

    ctωplot(mt, C, title="lim<sub><i>p→0</i></sub> C(<i>p,a,ω*,χ</i>)"; kwargs...)
end

function Lplot_fixed_ω(mt::MultiType; jp = 2, kwargs...)
    _, jj = findmin(mt.L_mat[:,:,2,:])

    jχ = jj[2]

    L = mt.L_mat[:, jχ, jp, :]
    title = "lim<sub><i>p→0</i></sub> 𝓛(<i>p,a,ω,χ*</i>)"

    ctaplot(mt, L, title = title; kwargs...)
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

function Lωplot(m0::Prequel; jp = 2, jA = 1, slides=true, dark=false, kwargs...)
    L = [minimum(m0.L[:, jχ, ja, jp, jA]) for ja in axes(m0.L, 3), jχ in axes(m0.L, 2)]

    jj = findfirst(L .== minimum(L))
    xmin = annualized(m0.agrid[jj[1]])
    ymin = annualized(m0.χgrid[jj[2]])

    shape_vec = [
        attr(; x0=xmin - 0.001, x1=xmin + 0.001, y0=ymin - 0.002, y1=ymin + 0.002, line_color="darkred", fillcolor="darkred", type="circle")
        attr(; x0=0, x1 = annualized(m0.χgrid[end]), y0=0, y1 = annualized(m0.χgrid[end]), type = "line", line_width = 1, line_dash="dash", line_color = "darkred")
    ]
    
    ctωplot(m0, L, shapes = shape_vec, 
        title = "lim<sub><i>p→0</i></sub> min<sub><i>ω</i></sub> 𝓛(<i>p,a,ω,χ</i>)", 
        slides=slides, dark=dark; kwargs...
    )
end

function Lωplot(mt::MultiType, L_mat = mt.L_mat; jp=2, slides=true, dark=false, kwargs...)
    L = [minimum(L_mat[:, jχ, jp, ja]) for ja in axes(L_mat, 4), jχ in axes(L_mat, 2)]

    jj = findfirst(L .== minimum(L))
    xmin = annualized(mt.ct.gr[:a][jj[1]])
    ymin = annualized(mt.χgrid[jj[2]])

    shape_vec = [
        attr(; x0=xmin - 0.001, x1=xmin + 0.001, y0=ymin - 0.002, y1=ymin + 0.002, line_color="darkred", fillcolor="darkred", type="circle")
        attr(; x0=0, x1 = annualized(mt.χgrid[end]), y0=0, y1 = annualized(mt.χgrid[end]), type = "line", line_width = 1, line_dash="dash", line_color = "darkred")
    ]
    
    ctωplot(mt, L, shapes = shape_vec, 
        title = "lim<sub><i>p→0</i></sub> min<sub><i>ω</i></sub> 𝓛(<i>p,a,ω,χ</i>)", 
        slides=slides, dark=dark; kwargs...)
end

ctωplot(m0::Prequel, y::Array; title = "", slides = true, dark = false, kwargs...) = ctωplot(m0, y, annualized.(m0.agrid), title = title, slides = slides, dark = dark; kwargs...)

function ctωplot(mt::Union{MultiType, Prequel}, y::Array, xgrid = annualized.(mt.ct.gr[:a]), ygrid = annualized.(mt.χgrid); title = "", slides = true, dark = false, kwargs...)
    Na, Nχ = length(xgrid), length(ygrid)
    @assert size(y) == (Na, Nχ)

    colpal = ColorSchemes.lapaz
    xt = "Initial inflation (<i>a<sub>0</sub></i>)"
    yt = "Asymptote (χ)"

    data = contour(
        z = y', x = xgrid, y = ygrid,
        colorscale = [[jj, get(colpal, 1-jj, :clamp)] for jj in range(0,1,length=50)]

    )
    layout = Layout(
        xaxis_title = xt, yaxis_title = yt, title = title,
        template = qtemplate(slides=slides, dark=dark);
        kwargs...
    )
    
    plot(data, layout)
end

function ctaplot(mt::MultiType, y::Array; slides = true, dark = false, kwargs...)

    colpal = ColorSchemes.lapaz
    xt = "Decay rate (%)"
    yt = "Initial inflation (<i>a<sub>0</sub></i>)"

    data = contour(
        z = y', x = perc_rate.(mt.ωgrid), y = annualized.(mt.ct.gr[:a]),
        colorscale = [[jj, get(colpal, 1-jj, :clamp)] for jj in range(0,1,length=50)]

    )
    layout = Layout(
        xaxis_title = xt, yaxis_title = yt, 
        template = qtemplate(slides=slides, dark=dark);
        kwargs...
    )
    
    plot(data, layout)
end 


function ctplot(mt::Union{MultiType, Prequel}, y::Array; slides = true, dark = false, kwargs...)
    
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
    Np = length(mt.ct.gr[:p])

	data = zeros(Np, 3)
    for jp in axes(data, 1)
        
        _, jj = findmin(mt.L_mat[:, :, jp, :])
        
        data[jp, 1] = perc_rate(mt.ωgrid[jj[1]])
        data[jp, 2] = annualized(mt.ct.gr[:a][jj[3]])
        data[jp, 3] = annualized(mt.χgrid[jj[2]])
    end
    
    datanames = ["ω", "a<sub>0", "χ"]
    yax = ["y2", "y1", "y1"]
    
    cols = [get(ColorSchemes.southwest, jj, :clamp) for jj in [0, 0.5, 1]]
    lines = [
        scatter(;x=mt.ct.gr[:p][2:end], y=data[2:end, jj], line_width = 2.5, yaxis=yax[jj], marker_color=cols[jj], name="<i>"*datanames[jj]*"</i>") for jj in eachindex(datanames)
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

function scatscol(z::AbstractMatrix, x::AbstractVector, y::AbstractVector; name = "", sc::ColorScheme=ColorSchemes.batlow, scmax = 0.95, kwargs...)
    @assert size(z) == (length(x), length(y))

    vv = eachindex(y) / length(y)
    m, M = extrema(y)

    col = get(sc, scmax * vv, :clamp)

    jx = round(Int, length(x)/2)
    xs = [x[jx] for _ in axes(z,2)]
    ys = [z[jx, jy] for jy in axes(z,2)]
    cols = range(m,M, length=length(xs))
    colscale = [[(jx-1) / (length(col) - 1), col[jx]] for jx in eachindex(col)]
    colnames = round.(range(m, M, length = 6), digits = 2)

    sc = [scatter(
        mode = "markers", marker_opacity = 0,
        x = xs, y = ys, showlegend=false,
        marker=attr(color=cols, reversescale=false, colorscale=colscale, colorbar=attr(tickvals=range(m, M, length=length(colnames)), title="&nbsp;<i>$name", ticktext=colnames))
    )]

    push!(sc, [
        scatter(x = x, y = z[:, jy], showlegend = false, name = "$name = $(@sprintf("%.3g", y[jy]))", marker_color = col[jy]; kwargs...) for jy in eachindex(y)
    ]...)

    return sc
end

function planspwide(mt::MultiType, dk = nothing; DKcrit = !isnothing(dk), slides = true, dark = false, T = 11, share = true)

    data = zeros(T, length(mt.ct.gr[:p]))

    for jp in axes(data, 2)
        if DKcrit
            _, jj = findmin(dk[:, :, jp, :])
        else
            _, jj = findmin(mt.L_mat[:, :, jp, :])
        end

        a = annualized(mt.ct.gr[:a][jj[3]])
        χ = annualized(mt.χgrid[jj[2]])
        ω = mt.ωgrid[jj[1]]
        
        for tt in 1:T
            data[tt, jp] = χ + exp(-ω*(tt-1)) * (a - χ)
        end
    end

    yaxis_title = "%"
    if share
        data .*= 1/annualized(Nash(mt))
        yaxis_title = "Share of Nash inflation"
    end

    layout = Layout(;
        template = qtemplate(;slides, dark), xaxis_title = "<i>Quarters", yaxis_title
    )

    plot(scatscol(data[:,2:end], 0:T-1, mt.ct.gr[:p][2:end], name="p"), layout)
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

function strategy_μ(mt::MultiType; save_stats=false, slides = true, dark = false, folder = "")

    χgrid = mt.χgrid
    ωgrid = mt.ωgrid
    agrid = mt.ct.gr[:a]

    marg_aχ = [sum(mt.μ[:, jχ, ja]) for jχ in axes(mt.μ, 2), ja in axes(mt.μ, 3)]

    marg_ωχ = [sum(mt.μ[jω, jχ, :]) for jω in axes(mt.μ, 1), jχ in axes(mt.μ, 2)]

    c1 = contour(x=annualized.(agrid), y=annualized.(χgrid), z=marg_aχ, colorscale=[[jj, get(ColorSchemes.lapaz, jj, :clamp)] for jj in range(0, 1, length=50)])

    c2 = contour(x=perc_rate.(ωgrid), y=annualized.(χgrid), z=marg_ωχ', colorscale=[[jj, get(ColorSchemes.lapaz, jj, :clamp)] for jj in range(0, 1, length=50)])

    P = sum([sum(mt.μ[:, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if agrid[ja] > χgrid[jχ]])
    print("P(a_0 > χ) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write(folder * "pa_chi.txt", "$(@sprintf("%0.3g",100P))\\%.")

    P = sum([sum(mt.μ[:, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if agrid[ja] > 5χgrid[jχ]])
    print("P(a_0 > 5χ) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write(folder * "pa_chi5.txt", "$(@sprintf("%0.3g",100P))\\%.")

    P = sum([sum(mt.μ[:, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid) if χgrid[jχ] == 0])
    print("P(a_0 > 0) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write(folder * "pa_chi0.txt", "$(@sprintf("%0.3g",100P))\\%.")

    P = sum([sum(mt.μ[jω, jχ, ja]) for jχ in 1:length(χgrid), ja in 1:length(agrid), jω in 1:length(ωgrid) if perc_rate(ωgrid[jω]) <= 10])
    print("P(ω ≤ 10%) = $(@sprintf("%0.3g",100P))%\n")
    save_stats && write(folder * "pa_omega0.txt", "$(@sprintf("%0.3g",100P))\\%.")

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

function get_Kambe(mt::MultiType; jp = 2)
	minL, jj = findmin(mt.L_mat[:, :, jp, :])
	ωK = mt.ωgrid[jj[1]]
	χK = mt.χgrid[jj[2]]
	aK = mt.ct.gr[:a][jj[3]]

    return ωK, χK, aK
end

function prep_plot_planner(mt::MultiType; k = 10)
	ct = mt.ct
	rp = Ramsey(ct)

	vfi!(rp, verbose = false)
	πR, θR = simul_plan(rp)

	# tvec = 1:length(πR)
    tvec = (0:k) .+ 1

	mean_ω, mean_a, mean_χ, sd_ω, sd_a, sd_χ = find_plan_μ(mt; annualize=false, decay=false)

	πC = (mean_a - mean_χ) * exp.(-mean_ω * (tvec.-1)) .+ mean_χ

    ωK, χK, aK = get_Kambe(mt)

	πK = (aK - χK) * exp.(-ωK * (tvec.-1)) .+ χK

    return tvec, πR, πC, πK
end

function get_Kambe(m0::Prequel, jA; jp = 2)
	minL, jj = findmin(m0.L[:, :, :, jp, jA])
	ωK = m0.ωgrid[jj[1]]
	χK = m0.χgrid[jj[2]]
	aK = m0.agrid[jj[3]]

    return ωK, χK, aK
end

function prep_plot_prequel(m0::Prequel; k = 10, jA = 1)
    tvec = (0:k) .+ 1
    
    ωK, χK, aK = get_Kambe(m0, jA)
	πK = (aK - χK) * exp.(-ωK * (tvec.-1)) .+ χK

    πK = vcat([m0.Agrid[jA]], πK)
    return 0:11, πK
end

function comp_plot_planner(mt::MultiType; k = 10, slides::Bool=true, dark = false, Kambe = true, Average = true, share = true, kwargs...)

    yaxis_title = "%"

    tvec, πR, πC, πK = prep_plot_planner(mt; k)



    if share
        yaxis_title = "Share of Nash inflation"
        πN = Nash(mt)
        πR *= 1/πN
        πC *= 1/πN
        πK *= 1/πN
    else
        πR = annualized.(πR)
        πC = annualized.(πC)
        πK = annualized.(πK)
    end

    layout = Layout(;
        title="Plans", xaxis_title="<i>Quarters", yaxis_title, legend=attr(orientation="h", x=0.05), template = qtemplate(slides = slides, dark = dark),
        kwargs...
    )
    data = [
        scatter(x=tvec.-1, y=πR[tvec], marker_color=get(ColorSchemes.southwest, 0.01, :clamp), name="<i>Ramsey", showlegend=true)
        ]
    if Kambe
        push!(data, scatter(x=tvec.-1, y=πK[tvec], marker_color=get(ColorSchemes.southwest, 0.5, :clamp), name="<i>Optimal eq'm"))
    end
    if Average
        push!(data, scatter(x=tvec.-1, y=πC[tvec], marker_color=get(ColorSchemes.southwest, 0.99, :clamp), name="<i>Average eq'm"))
    end

    plot(data, layout)
end

function show_μ(mt::MultiType; slides = true, dark = false)
	pgrid, agrid = mt.ct.gr[:p], mt.ct.gr[:a]
	ωgrid, χgrid = mt.ωgrid, mt.χgrid

    πN = Nash(mt.ct)
	tvec = 1:11

    L = mt.L_mat[:, :, 2, :]
    μ = (-L .+ maximum(L)) / (maximum(L) - minimum(L))

    marker = attr(color = get(ColorSchemes.oslo, 0.4, :clamp))

    sc = [
        scatter(x=tvec .- 1, y=100*((av - χv) * exp.(-ωv * tvec .- 1) .+ χv) / πN; showlegend=false, marker = get(ColorSchemes.oslo, 0.8 * (1-μ[jω, jχ, ja]), :clamp), hoverinfo="skip", line_width = 0.1, opacity=μ[jω, jχ, ja]^13, mode="lines") for (ja, av) in enumerate(agrid), (jω, ωv) in enumerate(ωgrid), (jχ, χv) in enumerate(χgrid)
    ] |> vec

    layout = Layout(template = qtemplate(; slides, dark),
        xaxis_title = "Quarters", yaxis_title = "Inflation (% of Nash)",
    )

    plot(sc, layout)
end

function comp_plot_planner(m0::Prequel, mt::MultiType; jA = 1, slides=true, dark=false)
    tvec, πR, πC, πK = prep_plot_planner(mt)

    tvec2, πK0 = prep_plot_prequel(m0; jA)

    layout = Layout(
        title="Plans", yaxis_title="%", xaxis_title="<i>Quarters",
        legend=attr(orientation="h", x=0.05),
        template = qtemplate(slides = slides, dark = dark)
    )
    data = [
        scatter(x=tvec.-1, y=annualized.(πR)[tvec], marker_color=get(ColorSchemes.southwest, 0.01, :clamp), name="<i>Ramsey")
        scatter(x=tvec.-1, y=annualized.(πC)[tvec], marker_color=get(ColorSchemes.southwest, 0.99, :clamp), name="<i>Average eq'm")
        scatter(x=tvec.-1, y=annualized.(πK)[tvec], marker_color=get(ColorSchemes.southwest, 0.5, :clamp), name="<i>Kambe eq'm")
        scatter(x=tvec2.-1, y=annualized.(πK0), name ="<i>Kambe eq'm with initial period", line_dash="dash")
    ]

    plot(data, layout)
end

function Lωplot(mt::Multiψ; jp = 2, slides = true, dark = false)

    L = [minimum(mt.L[:, :, jψ, jp, ja]) for jψ in axes(mt.L,3), ja in axes(mt.L, 5)]

    layout = Layout(
        template = qtemplate(slides=slides, dark=dark),
        xaxis_title = "Initial inflation (<i>a<sub>0</sub></i>)",
        yaxis_title = "ψ",
    )
    colpal = ColorSchemes.lapaz
    data = contour(
        z = L, y = mt.ψgrid, x = annualized.(mt.ct.gr[:a]),
        colorscale = [[jj, get(colpal, 1-jj, :clamp)] for jj in range(0,1,length=50)],
    )

    plot(data, layout)
end

function twolines(mt::Multiψ; show="L", jp = 2, slides = true, dark = false)
    cols = [get(ColorSchemes.southwest, x, :clamp) for x in (0.99, 0.5)]
    # L_reopt = [minimum(mt.L[:,:,jψ, jp,:]) for jψ in eachindex(mt.ψgrid)]
    L_reopt = zeros(length(mt.ψgrid))
    C_reopt = similar(L_reopt)
    for jψ in eachindex(mt.ψgrid)
        L, jj = findmin(mt.L[:,:,jψ, jp, :])
        L_reopt[jψ] = L
        C_reopt[jψ] = mt.C[jj[1], jj[2], jψ, jp, jj[3]]
    end

    _, jj = findmin(mt.L[:,:,1,jp,:])
    jω = jj[1]
    jχ = jj[2]
    ja = jj[3]

    L_og = mt.L[jω, jχ, :, jp, ja]
    C_og = mt.C[jω, jχ, :, jp, ja]

    data = [
        scatter(x=mt.ψgrid, y = L_reopt, marker_color=cols[1], name = "<b>a</b>*(ψ)")
        scatter(x=mt.ψgrid, y = L_og, marker_color=cols[2], name = "<b>a</b>*")
    ]
    if show == "C"
        data = [
            scatter(x=mt.ψgrid, y = C_reopt, marker_color=cols[1], name = "<b>a</b>*(ψ)")
            scatter(x=mt.ψgrid, y = C_og, marker_color=cols[2], name = "<b>a</b>*")
        ]
    end

    layout = Layout(
        template = qtemplate(slides=slides, dark=dark),
        xaxis_title = "<i>ψ"
    )

    plot(data, layout)
end

function implied_plan(mt::Multiψ; jp = 2, slides = true, dark = false)

    ψvec = mt.ψgrid

    ωvec = similar(ψvec)
    χvec = similar(ψvec)
    avec = similar(ψvec)

    for jψ in eachindex(ψvec)
        _, jj = findmin(mt.L[:,:,jψ,jp,:])
        ωvec[jψ] = mt.ωgrid[jj[1]]
        χvec[jψ] = mt.χgrid[jj[2]]
        avec[jψ] = mt.ct.gr[:a][jj[3]]
    end

    ωvec = perc_rate.(ωvec)
    χvec = annualized.(χvec)
    avec = annualized.(avec)

    cols = [get(ColorSchemes.southwest, jj, :clamp) for jj in [0, 0.5, 1]]

    lines = [
        scatter(x = ψvec, y = ωvec, line_width = 2.5, yaxis = "y2", marker_color = cols[1], name = "<i>ω*(ψ)")
        scatter(x = ψvec, y = avec, line_width = 2.5, yaxis = "y1", marker_color = cols[2], name = "<i>a<sub>0</sub>*(ψ)")
        scatter(x = ψvec, y = χvec, line_width = 2.5, yaxis = "y1", marker_color = cols[3], name = "<i>χ*(ψ)")
    ]

	layout = Layout(
        template = qtemplate(;slides, dark),
        hovermode = "x",
		yaxis = attr(domain=[0, 0.45], zeroline=false, title="<i>%"),
		yaxis2 = attr(domain=[0.55, 1], zeroline=false, title="<i>%"),
		xaxis = attr(zeroline=false, title="<i>ψ"),
		legend = attr(orientation="h", x=0.05),
		title="Preferred plans",
    )

    plot(lines, layout)
end

function implied_plan_wide(mt::Multiψ, dk = nothing; DKcrit=!isnothing(dk), jp = 2, slides = true, dark = false, share = true, T = 11)

    ψvec = mt.ψgrid

    data = zeros(T, length(ψvec))
    for jψ in axes(data, 2)
        if DKcrit
            _, jj = findmin(dk[:, :, jψ, jp, :])
        else
            _, jj = findmin(mt.L[:, :, jψ, jp, :])
        end

        a = annualized(mt.ct.gr[:a][jj[3]])
        χ = annualized(mt.χgrid[jj[2]])
        ω = mt.ωgrid[jj[1]]

        for tt in 1:T
            data[tt, jψ] = χ + exp(-ω * (tt - 1)) * (a - χ)
        end
    end

    yaxis_title = "%"
    if share
        data .*= 1 / annualized(Nash(mt))
        yaxis_title = "Share of Nash inflation"
    end

    vv = [jψ / length(ψvec) for jψ in eachindex(ψvec)]
    min_z, max_z = extrema(ψvec)
    col = get(ColorSchemes.batlow, 0.85 * vv, :clamp)

    jT = round(Int, T / 2)
    xs = [jT for _ in axes(data, 2)]
    ys = [data[jT, ja] for ja in axes(data, 2)]
    cols = range(0, 1, length=length(xs))
    colscale = [[(jT - 1) / (length(col) - 1), col[jT]] for jT in eachindex(col)]
    colnames = round.(range(min_z, max_z, length=6), digits=2)


    lines = [
        scatter()
        scatter(mode="markers", marker_opacity=0,
            x=xs, y=ys, showlegend=false,
            marker=attr(color=cols, reversescale=false, colorscale=colscale, colorbar=attr(tickvals=range(0, 1, length=length(colnames)), title="&nbsp;<i>ψ", ticktext=colnames))
        )
        [scatter(x=(1:T) .- 1, y=data[:, jψ], marker_color=col[jψ], showlegend=false) for jψ in axes(data, 2) if jψ > 1]
    ]
    layout = Layout(;
        template=qtemplate(; slides, dark), xaxis_title="<i>Quarters", yaxis_title, yaxis_range = [-0.01, maximum(data)*1.05]
    )

    # lines = scatscol(data[:,2:end], 0:T-1, ψvec[2:end], name=string(k))

    plot(lines, layout)
end




function allplans(mt::Multiψ; jp = 2, T = 11, slides = true, dark = false)

    ψvec = mt.ψgrid

    ωvec = similar(ψvec)
    χvec = similar(ψvec)
    avec = similar(ψvec)

    cvec = zeros(1:T, length(ψvec))

    for jψ in eachindex(ψvec)
        _, jj = findmin(mt.L[:,:,jψ,jp,:])
        ωvec[jψ] = mt.ωgrid[jj[1]]
        χvec[jψ] = mt.χgrid[jj[2]]
        avec[jψ] = mt.ct.gr[:a][jj[3]]

        a = annualized(avec[jψ])
        χ = annualized(χvec[jψ])
        Ω = exp(-ωvec[jψ])
        
        for tt in 1:T
            cvec[tt, jψ] = χ + Ω^(tt-1) * (a - χ)
        end
    end

    lines = [scatter(x = 0:T-1, y = cvec[:, jψ], name = "ψ = $(round(ψv, sigdigits=2))") for (jψ, ψv) in enumerate(ψvec)]

    layout = Layout(
        template = qtemplate(;slides, dark)
    )

    plot(lines, layout)
end

function reformat_x(xvec, k::Symbol)
    if k == :σ
        xvec = xvec * 400
    elseif k == :β
        xvec = xvec.^(-4)
    elseif k == :κ
    else
        throw(error("Wrong parameter name"))
    end
    xvec
end

function plot_compstats(k::Symbol, T = 11; share = false, slides = true, dark = false, showLmat = false)

    s = sk(k)

    ωvec, χvec, avec, xvec, Lmat, pvec = load("Output/JET/compstats_$s.jld2", "ωvec", "χvec", "avec", "$(s)vec", "Lmat", "pvec")

    K = length(xvec)
    plans = zeros(T, K)
    for jk in eachindex(xvec)

        a = annualized(avec[jk])
        χ = annualized(χvec[jk])
        Ω = exp(-ωvec[jk])
        
        for tt in 1:T
            plans[tt, jk] = χ + Ω^(tt-1) * (a - χ)
        end
    end

    xvec = reformat_x(xvec, k)

    yaxis_title = "%"
    if share
        plans .*= 1 / annualized(Nash(mt))
        yaxis_title = "Share of Nash inflation"
    end

    if showLmat
        return plot(
            contour(z = Lmat, x = pvec, y = xvec), 
            Layout(;template = qtemplate(;slides, dark), xaxis_title = "<i>p", yaxis_title = "<i>$(string(k))")
        )
    end

    plot(scatscol(plans, 0:T-1, xvec, name =string(k)), Layout(;
            template=qtemplate(; slides, dark), xaxis_title="<i>Quarters", yaxis_title, yaxis_range=[-0.01, maximum(plans) * 1.05]
        )
    )
end


findflex(y::Array{T, 2}; jp) where T = findmin(y[jp, :])[2]
findflex(y::Array{T, 3}; jp) where T = findmin(y[jp, :, :])[2]
findflex(y::Array{T, 4}; jp) where T = findmin(y[jp, :, :, :])[2]
findflex(y::Array{T, 5}; jp) where T = findmin(y[jp, :, :, :, :])[2]
findflex(y::Array{T, 6}; jp) where T = findmin(y[jp, :, :, :, :, :])[2]

function bestflex(z::Zero, mt=nothing; slides = true, dark = false, share = true)
    T = periods(z)

    a = zeros(T+4, T+4) * NaN

    y = z.L
    for (jt, tt) in enumerate(T:-1:2)

        jv = findflex(y, jp = 2)
        a_inv = annualized.(z.gr[:a][i] for i in Tuple(jv))

        a[periods(z)-tt+1:periods(z), jt] = a_inv[end:-1:1]
        
        if periods(z) < size(a,1)
            a[periods(z)+1:end, jt] .= a[periods(z), jt]
        end

        y = same_last_two(y)
    end

    yaxis_title = "%"
    if share
        yaxis_title = "Share of Nash inflation"
        a .*= 1/annualized(Nash(z))
    end

    sc = [
        [scatter(x=T-1:length(a), y=a[T+1:end, jj], name="<i>K</i> = $(periods(z)-jj+1)", marker_color=Owls(jj), legendgroup = jj, showlegend = false, line_dash="dot") for jj in 1:T-1]
        [scatter(x=axes(a, 1) .- 1, y=a[1:T, jj], name="<i>K</i> = $(periods(z)-jj)", marker_color=Owls(jj), legendgroup=jj) for jj in 1:T-1]
    ]
    My = maximum(a[.!(isnan.(a))])*1.05

    plans = zeros(T+3)
    if !isnothing(mt)
        _, jj = findmin(mt.L_mat[:, :, 2, :])
        a = annualized(mt.ct.gr[:a][jj[3]])
        χ = annualized(mt.χgrid[jj[2]])
        ω = mt.ωgrid[jj[1]]
        
        for tt in eachindex(plans)
            plans[tt] = χ + exp(-ω*(tt-1)) * (a - χ)
        end
        if share
            plans .*= 1/annualized(Nash(mt))
        end

        push!(sc,
            scatter(x = eachindex(plans).-1, y = plans, line_dash = "dash", marker_color = Owls(T), name = "Opt. announcement")
        )
    end
    Mp = maximum(plans * 1.05)

    My = max(My, Mp)

    plot(sc, Layout(; template=qtemplate(; slides, dark), xaxis_title="<i>Quarters", yaxis_title, yaxis_range = [-0.05, My]))
end
