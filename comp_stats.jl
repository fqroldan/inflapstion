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
function prepare_results(run_number, Nruns, smin, smax, param)
	if run_number == 1
		write("../../comments_compstats.txt", "Nruns = $Nruns. σ between $smin and $smax")
		write("../../output_compstats.txt", param*",omega,a,chi,Lmin\n")
	end
	nothing
end

run_number, Nruns = qload(ARGS)

include("ct.jl")

param = "sigma"
if param == "sigma"
	smin, smax = 0.0075, 0.02
elseif param == "beta"
	smin, smax = 0.99, 0.999
elseif param == "kappa"
	smin, smax = 0.1, 0.25
end
σvec = range(smin/4, smax/4, length=Nruns)
σs = σvec[run_number]

prepare_results(run_number, Nruns, smin, smax, param)

write(pwd()*"/../output.txt", "")
write(pwd()*"/../temp.txt", "")

function create_or_load(T::DataType, param)
	if param == "sigma"
		ct = CrazyType(T, σ = σs, ω = 0.2, χ = 0.0);
	elseif param == "beta"
		ct = CrazyType(T, β = σs, ω = 0.2, χ = 0.0);
	elseif param == "kappa"
		ct = CrazyType(T, κ = σs, ω = 0.2, χ = 0.0);
	end
	try
		print_save("Loading first file of previous run: ")
		ctt = load("../ct_1_temp.jld", "ct")
		if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
			ct.gπ = ctt.gπ
		end
		print_save("✓\n")
	catch
		print_save("failed. ")
		try
			print_save("Loading first file of generic run: ")
			ctt = load("../../ct_1.jld", "ct");
			if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
				ct.gπ = ctt.gπ
			end
			print_save("✓\n")
		catch
			print_save("failed. .\n")
			try
				print_save("Loading optimum of generic run: ")
				ctt = load("../../ct_opt.jld", "ct");
				if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
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

ct = create_or_load(Forward, param)

initial_report(ct)

function fill_in_results(par, ω, a, χ, Lmin)
	s = read("../../output_compstats.txt", String)
	write("../../output_compstats.txt", s*"$(par), $(ω), $(a), $(χ), $(Lmin)\n")
	nothing
end

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)
Nω = 30
Nχ = 15 
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na)

a, ω, χ = choose_ω!(L_mat, ct)
Lmin = minimum(L_mat[:,:,3,:])
fill_in_results(σs, ω, a, χ, Lmin)

nothing
