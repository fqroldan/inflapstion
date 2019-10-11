function print_save(s::String, temp=false)
	print(s)

	filename = pwd()*"/../output.txt"
	output = read(filename, String)

	if temp
		filename = pwd() * "/../temp.txt"
		output = read(filename, String)
	end

	write(filename, output * s)

	nothing
end

function time_print(tfloat::Float64)
	t = floor(Int, tfloat)
	if t < 1
		t_print = "no time"
	elseif t < 60
		t_print = "$(t) second"
		t == 1 ? nothing : t_print = t_print * "s"
	elseif t < 3600
		t_print = "$(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? nothing : t_print = t_print * "s"
		if t % 60 != 0
			t_print = t_print * " and $(t%60) second"
			t % 60 == 1 ? nothing : t_print = t_print * "s"
		end
	else
		t_print = "$(floor(Int,t/3600)) hour"
		floor(Int, t/3600) == 1 ? nothing : t_print = t_print * "s"
		t = t % 3600
		t_print = t_print * " and $(floor(Int,t/60)) minute"
		floor(Int,t/60) == 1 ? nothing : t_print = t_print * "s"
		# if t % 60 != 0
		# 	t_print = t_print * " and $(t%60) second"
		# 	t%60 == 1 ? nothing : t_print = t_print * "s"
		# end
	end
	return t_print
end

function initial_report(ct::CrazyType)

	π_Nash = ct.κ / (1.0 - ct.β + ct.κ^2*ct.γ) * ct.ystar
	π_Nash = (1+π_Nash)^4 - 1
	real_rate = (1/ct.β^4 - 1)*100
	print_save("Credibility Dynamics and Disinflation Plans\n")
	print_save("\nStarting run on $(Threads.nthreads()) threads at $(Dates.format(now(),"HH:MM"))")
	print_save("\nNash inflation is $(@sprintf("%.3g",100*π_Nash))%, real rate is $(@sprintf("%.3g",real_rate))%")
	print_save("\nGrid for 𝑎 goes up to $(@sprintf("%.3g",maximum(ct.agrid))) ($(@sprintf("%.3g",annualized(maximum(ct.agrid))))% annual)")

	nothing
end