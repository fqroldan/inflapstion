using Distributed, JLD2, FileIO

include("ct.jl")

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function create_or_load(T::DataType)
	ct = CrazyType(T, ω = 0.2, χ = 0.0);
	try
		print_save("Loading first file of previous run: ")
		ctt = load("../Output/ct_1.jld2", "ct")
		if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
			ct.gπ=ctt.gπ
		end
		print_save("✓\n")
	catch
		print_save("failed.")
		try
			print_save("Loading optimum of previous run: ")
			ctt = load("../Output/ct_opt.jld2", "ct");
			if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
				ct.gπ=ctt.gπ
			end
			print_save("✓\n")
		catch
			print_save("failed. Using new type as guess.\n")
		end
	end
	return ct
end

function makeplots_mimics_marginals(mt::MultiType)
	find_equil!(mt)
	for slides in [true, false]
		slides ? sty = slides_def : sty = paper
		slides ? yh = 0.7 : yh = 1
		p1 = strategy_μ(mt, style=sty, yh = yh)
		savejson(p1, pwd()*"/../Graphs/tests/marg$(ifelse(slides, "_slides", "_paper")).json")
		savefig(p1, pwd()*"/../Graphs/tests/marg$(ifelse(slides, "_slides", "_paper")).pdf")
	end
	for slides in [true,false]
		for CI in [true, false]
			save_plot_mimic_z(mt, CIs=CI, slides=slides)
		end
	end
	nothing
end
function makeplots_planner(mt::MultiType)
	find_equil!(mt)
	comp_plot_planner(mt, makeplots=true)
	print_save("\nDone with plot vs planner")
	make_sustainable_plots(mt, 50, makeplots=true, pc = Fwd_strategy)
	print_save("\nDone with plot vs sustainable plans")
	make_sustainable_plots(mt, 50, makeplots=true, pc = Fwd_GP)
	print_save("\nDone with plot vs sust plans with reverting triggers")
	nothing
end


ct = create_or_load(Forward)
initial_report(ct)

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 20
Nχ = 20
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na);
a, ω, χ, mt = choose_ω!(L_mat, ct)

find_equil!(mt)
save("../Output/mt.jld2", "mt", mt)

makeplots_mimics_marginals(mt)

makeplots_planner(mt)
