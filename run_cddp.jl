include("ctRoot.jl")
include("prequel.jl")


mt = load("Output/mt.jld2", "mt")

z = makeZero(mt, Na = 30)

solve!(z)

extend(z) # T = 3
extend(z) # T = 4

save("z.jld2", "z", z)