using PlotlyJS, JSON


function retrieve_plot(fn::String)

	pl = JSON.parsefile(PlotlyJS.Plot, fn)

	return plot(pl)
end
