using PumpProbeModels
# using PyPlot
using BenchmarkTools
# pygui(false)

wavenumbers = 2800:3000
delaytimes = 0:100

model = PumpProbeModel(wavenumbers, delaytimes, 2)
model.parameters.A = [3.0, 3.0]
model.parameters.a = [0.4, 0.4]
model2 = PumpProbeModel(wavenumbers, delaytimes, 2; pumpmode=:triangle)
model2.parameters = model.parameters
model2.pumptimes = -1:0.1:1

F = Array{Float64,2}(undef, length(delaytimes), length(wavenumbers))
w,d,G = PumpProbeModels.convert1D(wavenumbers,delaytimes,F)
H = deepcopy(G)

# b = @benchmark PumpProbeModels.evaluate1!($F, $model)
#
# @assert F[40,20] == -0.000318328733471701

# c = @benchmark PumpProbeModels.evaluate!($G, $w, $d, $model)
# d = @benchmark PumpProbeModels.evaluate!($H, $w, $d, $model2)

@btime PumpProbeModels.evaluate!(G, w, d, model)
@btime PumpProbeModels.evaluate!(H, w, d, model2)

X = reshape(G, size(F))
Y = reshape(H, size(F))
@show X[40,20]
@show Y[40,20]
@assert X[40,20] == -0.002614753495806977
@assert Y[40,20] == -0.0026144358299938697

XN = X' ./ minimum(X)
YN = Y' ./ minimum(Y);

# figure()
# subplot(121)
# pcolormesh(YN); title("YN")
# subplot(122)
# pcolormesh(XN); title("XN")

nothing
