module PumpProbeModels

export PumpProbeModel, evaluate!

mutable struct ParamContainer
    A::Array{Real,1}
    a::Array{Real,1}
    τ::Array{Array{Real,1}}
    ω::Array{Real,1}
    Δω::Array{Real,1}
    σ::Array{Array{Real,1}}
    Γ::Array{Real,1}
    ΔΓ::Array{Real,1}
    γ::Array{Array{Real,1}}
    t0::Real
end

mutable struct PumpProbeModel
    wavenumbers::AbstractArray{Real,1}
    delaytimes::AbstractArray{Real,1}
    parameters::ParamContainer
    decaymode::Symbol
    pumpmode::Symbol
    pumpfunction::Function
    pumptimes::AbstractArray{Real,1}
end

"""
"""
function PumpProbeModel(wavenumbers, delaytimes, N::Int; decaymode=:exp, pumpmode=:δ)
    A = fill(1.0, N)
    a = fill(0.3, N)
    # Distribute the resonance positions evenly across the wavenumber span as a default
    wnspan = maximum(wavenumbers) - minimum(wavenumbers)
    ω = fill(0.0, N)
    for i = 1:N
        ω[i] = i*wnspan / (N+1) + minimum(wavenumbers)
    end
    Δω = fill(0.0, N)
    σ = [fill(Inf, N)]
    Γ = fill(6.0, N)
    ΔΓ = fill(0.0, N)
    γ = [fill(Inf, N)]

    if decaymode == :exp
        # For τ put in the delay span for every peak as a default
        delayspan = maximum(delaytimes) - minimum(delaytimes)
        τ = [fill(delayspan/2, N)]
    else
        error("No decaymode $decaymode available")
    end

    if pumpmode == :δ
        # just put one value in t0
        t0 = 0.0
        pumpfunction = t -> (t == t0) ? 1.0 : 0.0
        pumptimes = [0.0]
    elseif pumpmode == :triangle
        t0 = 30.0
        pumpfunction = t -> (t <= 0) ? (t+20)/20 : -(t-10)/10
        pumptimes = -20:1:10
    else
        error("No pumpmode $pumpmode available")
    end

    parameters = ParamContainer(A, a, τ, ω, Δω, σ, Γ, ΔΓ, γ, t0)
    model = PumpProbeModel(wavenumbers, delaytimes, parameters, decaymode, pumpmode, pumpfunction, pumptimes)
end


function evaluate!(F::T, x, t, m::PumpProbeModel) where T

    F .= zero(eltype(T))
    r = fill(0.0+0.0im, size(F))
    s = fill(0.0+0.0im, size(F))

    function compute_spectra!(r, s, m::PumpProbeModel, x, t, t_step)

        function myexp(x::T) where T
            if x > 0
                return zero(T)
            else
                return exp(x)
            end
        end

        A = m.parameters.A
        ω = m.parameters.ω
        Γ = m.parameters.Γ
        a = m.parameters.a
        τ = m.parameters.τ
        σ = m.parameters.σ
        γ = m.parameters.γ
        t0 = m.parameters.t0 + t_step
        Δω = m.parameters.Δω
        ΔΓ = m.parameters.ΔΓ

        for i = 1:length(A)
            @. r += A[i] / (x - ω[i] - 1im * Γ[i])
            # Checking for zero values prevents computing unneccesary terms
            if iszero(Δω[i]) && iszero(ΔΓ[i])
                # println("both zero")
                @. s += A[i] * (1 - a[i] * myexp(-(t-t0) / τ[1][i])) /
                        (x - ω[i] - 1im * Γ[i])
            elseif iszero(ΔΓ[i])
                # println("gamma zero")
                @. s += A[i] * (1 - a[i] * myexp(-(t-t0) / τ[1][i])) /
                        (x - ω[i] + Δω[i] * myexp(-(t-t0) / σ[1][i]) -
                        1im * Γ[i])
            elseif iszero(Δω[i])
                # println("omega zero")
                @. s += A[i] * (1 - a[i] * myexp(-(t-t0) / τ[1][i])) /
                        (x - ω[i] -
                        1im * Γ[i] + ΔΓ[i] * myexp(-(t-t0) / γ[1][i]))
            elseif !iszero(Δω[i]) && !iszero(ΔΓ[i])
                # println("none zero")
                @. s += A[i] * (1 - a[i]   * myexp(-(t-t0) / τ[1][i])) /
                        (x - ω[i]  + Δω[i] * myexp(-(t-t0) / σ[1][i]) -
                        1im * Γ[i] + ΔΓ[i] * myexp(-(t-t0) / γ[1][i]))
            else
                error("🍣")
            end
        end
        return r, s
    end

    for t_step in m.pumptimes

        # We are in principle computing a bunch of delta functions here
        # and adding a weight to them
        pumpintensity = m.pumpfunction(t_step)
        r .= 0.0
        s .= 0.0
        compute_spectra!(r, s, m, x, t, t_step)
        @. F += (abs2(s) - abs2(r)) * pumpintensity
    end
    # Integrate Intesity and normalize
    @show intensity_total = m.pumpfunction.(m.pumptimes) |> sum
    F ./= intensity_total  #Renormalize
    return F
end

function convert1D(x, y, z)
    xarray = Float64[]
    yarray = Float64[]
    zarray = Float64[]
    for _x in x, _y in y
        push!(xarray, _x)
        push!(yarray, _y)
    end
    for i = 1:size(z,2)
        push!(zarray, z[:,i]...)
    end
    xarray, yarray, zarray
end

end # module
