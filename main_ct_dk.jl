using Distributed, JLD

include("ct.jl")

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function create_or_load(T::DataType)
	ct = CrazyType(T, ω = 0.2, χ = 0.0, Np=20, Na=20);
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

ct = create_or_load(Simultaneous)

dk = DovisKirpalani(ct);

initial_report(dk)

Epfi!(dk, tol=1e-4, tempplots=true, upd_η = 0.1)
save("../../dk.jld", "dk", dk)

dk = switch_PC(dk, Simultaneous);

Epfi!(dk, tol=1e-4, tempplots=true, upd_η = 0.1)
save("../../dk_simultaneous.jld", "dk", dk)

nothing