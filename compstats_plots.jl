using CSV, DataFrames, PlotlyJS, ColorSchemes

function makeplot_compstats(param::String; slides::Bool=true)

	df = CSV.read("../HPC_output/compstats_"*param*".txt");
	df = sort(df, Symbol(param));

	varnames = ["σ"; "a<sub>0"; "ω"; "χ"; "L"]

	palette = ColorSchemes.southwest

	yax = [""; "y1"; "y2"; "y1"; "y3"]

	ps = [scatter(x=df[!,1]*4, y=df[!,jj], marker_color=get(palette, (jj-1)/(length(yax)-1)), name=varnames[jj], xaxis="x", yaxis=yax[jj]) for jj in 2:5]
	# pω = scatter(x=df[!,1]*4, y=df[!,3], name=varnames[3], xaxis="x", yaxis="y2")

	layout = Layout(
		xaxis=attr(domain=[0,1],anchor="y", title="<i>$(varnames[1])"),
		yaxis1=attr(domain=[0,0.33]),
		yaxis2=attr(domain=[0.33,0.66]),
		yaxis3=attr(domain=[0.67,1], anchor="x"),
		legend=attr(orientation="h", x=0.05)
		)

	p1 = plot(ps[[1;3;2;4]], layout)

	if slides
		relayout!(p1, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", font_family="Lato", fontsize=16, width=700, height=400)
	else
		relayout!(p1, font_family="Linux Libertine", fontsize=16, width=1200, height=900)
	end

	return p1
end
