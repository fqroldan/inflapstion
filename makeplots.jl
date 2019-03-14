using PlotlyJS, JSON


function retrieve_plot(fn::String)

	pl = JSON.parsefile(PlotlyJS.Plot, fn)

	return plot(pl)
end

function save_summ_plots(N=10, ext="png")

	for jj in 1:N
		p1 = retrieve_plot(pwd()*"/../Graphs/tests/summary_jom_$(jj).json")
		relayout!(p1, height=800, width=1200)
		savefig(p1, pwd()*"/../Graphs/tests/summary_$(jj)."*ext)
	end

	nothing
end

