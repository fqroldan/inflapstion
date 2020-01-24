using CSV, DataFrames, PlotlyJS

function makeplot_compstats(param::String)

	df = CSV.read("../HPC_output/compstats_"*param*".txt");
	df = sort(df, Symbol(param))

	names = ["σ"; "a<sub>0"; "ω"; "χ"]

	paχ = [scatter(x=df[!,1]*4, y=df[!,jj], name=names[jj], xaxis="x", yaxis="y1") for jj in [2;4]]
	pω = scatter(x=df[!,1]*4, y=df[!,3], name=names[3], xaxis="x", yaxis="y2")

	layout = Layout(
		xaxis=attr(domain=[0,1],anchor="y", title="<i>$(names[1])"),
		yaxis1=attr(domain=[0,0.45]),
		yaxis2=attr(domain=[0.55,1], anchor="x"),
		legend=attr(orientation="h", x=0.05)
		)

	p1 = plot([pω; paχ], layout)





	pv = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 2*Nv)

	for jj in 1:Nv
		weight = (jj-1)/(Nv-1)
		pv[2*jj-1] = scatter(;x=(1:T)/1, y=X_vec[jj,:], line_dash="solid", name="true " * names[jj], marker_color=get(colorpal, weight), xaxis="x", yaxis="y$(jj)")
		pv[2*jj] = scatter(;x=(1:T)/1, y=x_vec[jj,:], line_dash="dashdot", name="estimated " * names[jj], marker_color=get(colorpal, weight), xaxis="x", yaxis="y$(jj)")
	end
	layout=Layout(
		xaxis=attr(domain=[0,1], anchor="y", title="<i>Years"),
		yaxis_domain=[0,0.33],
		yaxis2_domain=[0.33,0.67],
		yaxis3=attr(domain=[0.67,1], anchor="x")
		)
