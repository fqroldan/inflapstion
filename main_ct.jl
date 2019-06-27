using Distributed

@everywhere include("ct.jl")

write(pwd()*"/../output.txt", "")

function establish_remote()
	machine_name = ""
	try
		machine_name = read("/name.txt", String)
	catch
	end
	return !(machine_name=="qlaptop")
end
machine_remote = establish_remote()

ct = CrazyType(ω = 0.2, χ = 0.0);
initial_report(ct)
# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 22
Nχ = 20
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na)
ωmin = choose_ω!(L_mat, ct; remote = machine_remote)


# for (jp, pv) in enumerate(ct.pgrid)
# 	findmin(L_mat[:,:,jp,:])
# end