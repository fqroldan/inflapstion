using PlotlyJS, JSON


function retrieve_plot(fn::String)

	pl = JSON.parsefile(PlotlyJS.Plot, fn)

	return plot(pl)
end

function save_summ_plots(N=10, ext="png")

	for jj in 1:N
		p1 = retrieve_plot(pwd()*"/../Graphs/tests/summary_jom_$(jj).json")
		relayout!(p1, height=800, width=1200, font_size=18)
		savefig(p1, pwd()*"/../Graphs/tests/summary_$(jj)."*ext)
		p2 = retrieve_plot(pwd()*"/../Graphs/tests/simul_jom_$(jj).json")
		relayout!(p2, height=800, width=1200, font_family = "Fira Sans Light", font_size = 18, plot_bgcolor="rgba(250, 250, 250, 1.0)", paper_bgcolor="rgba(250, 250, 250, 1.0)")
		savefig(p2, pwd()*"/../Graphs/tests/simul_$(jj)."*ext)
	end

	nothing
end

function save_all_plots()

	plnames = ["first_L", "first_g", "first_p", "plans", "mimics", "mimics_CI", "contour", "Ccontour"]

	for slides in [true, false]
		for pln in plnames
			plname = pln*"$(ifelse(slides, "_slides", "_paper"))"
			p1 = retrieve_plot("../HPC_output/Graphs/"*plname*".json")
			savefig(p1, "../Graphs/new/"*plname*".pdf", format="pdf")
		end
	end
	nothing
end
		