using Distributed, JLD

include("ct.jl")

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function establish_remote()
	machine_name = ""
	try
		machine_name = read("/name.txt", String)
	catch
	end
	return !(machine_name=="qlaptop")
end
machine_remote = establish_remote()

function create_or_load(T::DataType)
	ct = CrazyType(T, ω = 0.2, χ = 0.0);
	try
		print_save("Loading first file of previous run: ")
		ctt = load("../../ct_1.jld", "ct")
		if typeof(ctt) == typeof(ct)
			ct=ctt
		end
		print_save("✓\n")
	catch
		print_save("failed.")
		try
			print_save("Loading optimum of previous run :")
			ctt = load("../../ct_opt.jld", "ct");
			if typeof(ctt) == typeof(ct)
				ct=ctt
			end
			print_save("✓\n")
		catch
			print_save("failed. Using new type as guess.\n")
		end
	end
	return ct
end

ct = create_or_load(Backward)

initial_report(ct)

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 20
Nχ = 15 
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na)
ωmin = choose_ω!(L_mat, ct; remote = machine_remote)


# for (jp, pv) in enumerate(ct.pgrid)
# 	findmin(L_mat[:,:,jp,:])
# end