using Distributed, JLD

include("ct.jl")

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function create_or_load(T::DataType)
	ct = CrazyType(T, ω = 0.2, χ = 0.0);
	try
		print_save("Loading first file of previous run: ")
		ctt = load("../../ct_1.jld", "ct")
		if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
			ct.gπ=ctt.gπ
		end
		print_save("✓\n")
	catch
		print_save("failed.")
		try
			print_save("Loading optimum of previous run: ")
			ctt = load("../../ct_opt.jld", "ct");
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

ct = create_or_load(Forward)

initial_report(ct)

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 60
Nχ = 30
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na)
a, ω, χ, mt = choose_ω!(L_mat, ct)

find_equil!(mt)

for slides in [true, false]
	for CI in [true, false]
		save_plot_mimic_z(mt, CIs=CI, slides=slides)
		save("../../mt.jld", "mt", mt)
	end
	p1, p2 = strategy_μ(mt, slides=slides)
	savejson(p1, pwd()*"/../Graphs/tests/marg_achi$(ifelse(slides, "_slides", "_paper")).json")
	savejson(p2, pwd()*"/../Graphs/tests/marg_omegachi$(ifelse(slides, "_slides", "_paper")).json")
end

nothing


# for (jp, pv) in enumerate(ct.pgrid)
# 	findmin(L_mat[:,:,jp,:])
# end