using Distributed, JLD

function qload(s)
	if length(s) == 2
		s1 = parse(Int, s[1])
		s2 = parse(Int, s[2])
	else
		throw(error("ARGS is too long (length = $(length(s))"))
	end

	return s1, s2
end
function prepare_results(run_number, Nruns, smin, smax)
	if run_number == 1
		write("../../comments_compstats.txt", "Nruns = $Nruns. σ between $smin and $smax")
		write("../../output_compstats.txt", "sigma,a,omega,chi,Lmin\n")
	end
	nothing
end

run_number, Nruns = qload(ARGS)

include("ct.jl")
smin, smax = 0.0075, 0.02
σvec = range(smin/4, smax/4, length=Nruns)
σs = σvec[run_number]

prepare_results(run_number, Nruns, smin, smax)

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function create_or_load(T::DataType)
	ct = CrazyType(T, σ = σs, ω = 0.2, χ = 0.0);
	try
		print_save("Loading first file of previous run: ")
		ctt = load("../ct_1_temp.jld", "ct")
		if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
			ctt.σ = σs
			ct.gπ = ctt.gπ
		end
		print_save("✓\n")
	catch
		print_save("failed. ")
		try
			print_save("Loading first file of generic run: ")
			ctt = load("../../ct_1.jld", "ct");
			if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
				ctt.σ = σs
				ct.gπ = ctt.gπ
			end
			print_save("✓\n")
		catch
			print_save("failed. .\n")
			try
				print_save("Loading optimum of generic run: ")
				ctt = load("../../ct_opt.jld", "ct");
				if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
					ctt.σ = σs
					ct.gπ = ctt.gπ
				end
				print_save("✓\n")
			catch
				print_save("failed. Using new type as guess.\n")
			end
		end
	end
	return ct
end

ct = create_or_load(Forward)

initial_report(ct)

function fill_in_results(σ, a, ω, χ, Lmin)
	s = read("../../compstats.txt", String)
	write("../../compstats.txt", s*"$(σ), $(a), $(ω), $(χ), $(Lmin)\n")
	nothing
end

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 30
Nχ = 15 
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na)

a, ω, χ = choose_ω!(L_mat, ct)
Lmin = minimum(L_mat[:,:,3,:])
fill_in_results(σs, a, ω, χ, Lmin)

nothing
