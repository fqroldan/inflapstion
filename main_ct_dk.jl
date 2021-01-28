using Distributed, JLD2, FileIO

include("ct.jl")

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function create_or_load(T::DataType)
	ct = CrazyType(T, ω = 1.5, χ = 0.0);
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

# ct = CrazyType(Forward, ω = 1.5, χ = 0.0, Np=20, Na=20);
# update_ga!(ct)

ct = create_or_load(Forward)

dk = DovisKirpalani(ct);

initial_report(dk)

solve!(dk)
save("../Output/dk.jld2", "dk", dk)

dk = switch_PC(dk, Simultaneous);

solve!(dk, tol=4e-4)
save("../Output/dk_simultaneous.jld2", "dk", dk)

nothing