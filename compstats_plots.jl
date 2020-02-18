using CSV, DataFrames, PlotlyJS, ColorSchemes

function makeplot_compstats(param::String; slides::Bool=true, temp::Bool=false)

	if temp 
		df = CSV.read("../HPC_output/output_compstats.csv");
	else
		df = CSV.read("../HPC_output/compstats_"*param*".csv");
	end
	df = sort(df, Symbol(param));

	varnames = ["<i>σ</i>"; "<i>ω</i>"; "<i>a<sub>0</sub></i>"; "<i>χ</i>"; "𝓛"]

	palette = ColorSchemes.southwest

	if param == "sigma"
		xax = "σ"
	elseif param == "beta"
		xax = "β"
	elseif param == "kappa"
		xax = "κ"
	end

	yax = [""; "y2"; "y1"; "y1"]

	ps = [scatter(x=df[!,1], y=df[!,jj], mode="lines", marker_color=get(palette, (jj-1)/(length(yax)-1)), name="𝔼["*varnames[jj]*"]", xaxis="x", yaxis=yax[jj]) for jj in 2:length(yax)]
	# pω = scatter(x=df[!,1]*4, y=df[!,3], name=varnames[3], xaxis="x", yaxis="y2")

	layout = Layout(
		title="Average plans",
		xaxis=attr(domain=[0,1],anchor="y", title="<i>"*xax),
		yaxis1=attr(domain=[0,0.45], title="%"),
		yaxis2=attr(domain=[0.55,1], anchor="x"),
		legend=attr(orientation="h", x=0.05)
		)

	p1 = plot(ps, layout)

	p2 = plot(scatter(x=df[!,1], y=df[!,end], mode="lines", fill="tozeroy"), Layout(xaxis = attr(title="<i>"*xax, zeroline=false), yaxis_zeroline=false))
	
	if slides
		relayout!(p1, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", font_family="Lato", font_size=16, width=900, height=600)
		relayout!(p2, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", font_family="Lato", font_size=16, width=900, height=600)
	else
		relayout!(p1, font_family="Linux Libertine", font_size=16, width=700, height=350)
		relayout!(p2, font_family="Linux Libertine", font_size=16, width=700, height=350)
	end

	return p1, p2
end