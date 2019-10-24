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

function create_or_load()
	ct = CrazyType(ω = 0.2, χ = 0.0);
	try
		print_save("trying for ct_1")
		ct = load("../../ct_1.jld", "ct")
		print_save("Loaded first file of previous run")
	catch
		try
			print_save("failed. ")
			ct = load("../../ct_opt.jld", "ct");
			print_save("Loaded optimum ω of previous run")
		catch
		end
	end
	return ct
end

ct = create_or_load()

initial_report(ct)

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 15
Nχ = 15 
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na)
ωmin = choose_ω!(L_mat, ct; remote = machine_remote)


# for (jp, pv) in enumerate(ct.pgrid)
# 	findmin(L_mat[:,:,jp,:])
# end