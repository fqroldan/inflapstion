using JLD2

function qload(s)
	if length(s) == 3
		s1 = parse(Int, s[1])
		s2 = parse(Int, s[2])
		s3 = s[3]
	else
		throw(error("ARGS is too long (length = $(length(s))"))
	end

	return s1, s2, s3
end
function prepare_results(run_number, Nruns, smin, smax, param)
	if run_number == 1
		write("../comments_compstats.txt", "Nruns = $Nruns. $param between $smin and $smax")
		write("../output_compstats.csv", param*",omega,a,chi,Lmin\n")
	end
	nothing
end

run_number, Nruns, param = qload(ARGS)

run_number == 1 && println("Starting $Nruns runs for parameter $param")

include("../../../../Codes/ct.jl")

if param == "sigma"
	smin, smax = 0.75/400, 1.25/400
elseif param == "beta"
	smin, smax = 1.01^(-0.25), 1.05^(-0.25)
elseif param == "kappa"
	smin, smax = 0.1, 0.25
end
σvec = range(smin, smax, length=Nruns)
σs = σvec[run_number]

if param == "sigma"
	show_σs = σs * 400
elseif param == "beta"
	show_σs = 100 * (σs^-4 - 1)
else
	show_σs = σs
end

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
		ctt = load("../../../ct_1.jld2", "ct")
		if typeof(ctt) == typeof(ct) && ct.Np == ctt.Np && ct.Na == ctt.Na
			ct.gπ=ctt.gπ
		end
		print_save("✓\n")
	catch
		print_save("failed.")
		try
			print_save("Loading optimum of previous run: ")
			ctt = load("../../../ct_opt.jld2", "ct");
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

ct = create_or_load(Forward, param)

initial_report(ct)

function fill_in_results(par, ω, a, χ, Lmin)
	s = read("../output_compstats.csv", String)
	write("../output_compstats.csv", s*"$(par), $(ω), $(a), $(χ), $(Lmin)\n")
	nothing
end

# Epfi!(ct, tol=1e-4, tempplots=true, upd_η = 0.1)

Nω = 20
Nχ = 20
print_save("\nNω, Nχ = $Nω, $Nχ")
L_mat = zeros(Nω, Nχ, ct.Np, ct.Na);
a, ω, χ, mt = choose_ω!(L_mat, ct, verbose = false)

# Lmin = minimum(L_mat[:,:,3,:])
# fill_in_results(show_σs, ω, a, χ, Lmin)

z0 = mt.ct.pgrid[3]
Lavg = find_equil!(mt, z0)
mean_ω, mean_χ, mean_a, sd_ω, sd_χ, sd_a = find_plan_μ(mt)
fill_in_results(show_σs, mean_ω, mean_a, mean_χ, Lavg)


nothing

