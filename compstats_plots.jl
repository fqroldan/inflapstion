using CSV, DataFrames, PlotlyJS, ColorSchemes

function makeplot_compstats(param::String; slides::Bool=true, temp::Bool=false)

	if temp 
		df = CSV.read("../HPC_output/output_compstats.csv");
	else
		df = CSV.read("../HPC_output/compstats_"*param*".csv");
	end
	df = sort(df, Symbol(param));

	varnames = ["<i>σ"; "<i>ω"; "<i>a<sub>0"; "<i>χ"]


	palette = ColorSchemes.southwest

	yax = [""; "y2"; "y1"; "y1"]

	ps = [scatter(x=df[!,1]*400, y=df[!,jj], marker_color=get(palette, (jj-1)/(length(yax)-1)), name=varnames[jj], xaxis="x", yaxis=yax[jj]) for jj in 2:length(yax)]
	# pω = scatter(x=df[!,1]*4, y=df[!,3], name=varnames[3], xaxis="x", yaxis="y2")

	layout = Layout(
		xaxis=attr(domain=[0,1],anchor="y", title="<i>$(varnames[1])"),
		yaxis1=attr(domain=[0,0.45], title="%"),
		yaxis2=attr(domain=[0.55,1], anchor="x"),
		legend=attr(orientation="h", x=0.05)
		)

	p1 = plot(ps, layout)

	if slides
		relayout!(p1, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa", font_family="Lato", font_size=16, width=900, height=600)
	else
		relayout!(p1, font_family="Linux Libertine", font_size=16, width=700, height=400)
	end

	return p1
end
