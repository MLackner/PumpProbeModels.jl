# using PumpProbeModels
using PyPlot
using BenchmarkTools
pygui(false)

wavenumbers = 2800:3000
delaytimes = 0:100

model = PumpProbeModel(wavenumbers, delaytimes, 2)
model2 = PumpProbeModel(wavenumbers, delaytimes, 3; pumpmode=:triangle)
model2.parameters.A = [0.8, 1.0, 0.9]
model2.parameters.τ = [[200.0, 10, 90]]
model2.parameters.Δω = [0.0, 0.0, -20]
model2.parameters.σ = [[0.0, 0.0, 20.0]]

F = Array{Float64,2}(undef, length(delaytimes), length(wavenumbers))
w,d,G = PumpProbeModels.convert1D(wavenumbers,delaytimes,F)
H = deepcopy(G)

# b = @benchmark PumpProbeModels.evaluate1!($F, $model)
#
# @assert F[40,20] == -0.000318328733471701

# c = @benchmark PumpProbeModels.evaluate!($G, $w, $d, $model)
# d = @benchmark PumpProbeModels.evaluate!($H, $w, $d, $model2)

PumpProbeModels.evaluate!(G, w, d, model)
PumpProbeModels.evaluate!(H, w, d, model2)

X = reshape(G, size(F))
Y = reshape(H, size(F))
@show X[40,20]
@assert X[40,20] == -0.0002233944636494317

XN = X' ./ minimum(X)
YN = Y' ./ minimum(Y)

figure()
subplot(121)
pcolormesh(YN); title("YN")
subplot(122)
pcolormesh(XN); title("XN")

c, d
