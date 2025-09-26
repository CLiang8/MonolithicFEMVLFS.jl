using Revise
using Gridap
using Printf
using Plots
using DrWatson
using WaveSpec
using .Constants
using .WaveTimeSeries
using .Jonswap
using WriteVTK
using DataFrames:DataFrame
using DataFrames:Matrix


name::String = "data/sims_202509/trial_dsn"
order::Int = 2
vtk_output::Bool = true
filename = name*"/mem"
mkpath(filename)


H0 = 10 #m #still-water depth


# Wave parameters 
ω, S, η₀ = jonswap(0.4, 2.5; 
    plotflag=true, plotloc=filename, nω=200)
println(ω[1], "\t", ω[2], "\t", ω[end])
ω = ω[2:end]
S = S[2:end]
η₀ = η₀[2:end]
α = randomPhase(ω; seed=100)
# ω = [2*π/2.0, 2*π/2.5, 2*π/4]
# η₀ = [0.3, 0.4, 0.1]
# @show α = randomPhase(ω; seed=100)

# 最主要频率以及波高（对应S最大）
# === [新增] 打印 S 最大值对应的 ω 和 η₀ ===
i_peak = argmax(S)
println("Peak of S(ω): ω = ", ω[i_peak], " rad/s, η₀ = ", η₀[i_peak], " m")

# 累计能量（支持不等距网格）
function cumtrapz(ω, y)
    n = length(ω)
    E = zeros(Float64, n)
    for i in 2:n
        Δ = ω[i] - ω[i-1]
        E[i] = E[i-1] + 0.5*(y[i] + y[i-1]) * Δ
    end
    return E
end

E   = cumtrapz(ω, S)
m0  = E[end]
cov = 0.95                     # 想更紧就 0.95；更宽就 0.99
El  = (1 - cov)/2 * m0         # 左尾能量
Er  = (1 + cov)/2 * m0         # 右侧累计到这里
iL  = findfirst(>=(El), E)
iR  = findfirst(>=(Er), E)

ω_sel  = ω[iL:iR]
S_sel  = S[iL:iR]
η₀_sel = η₀[iL:iR]


# Wave parameters
k = dispersionRelAng.(H0, ω; msg=false)
@show λ = 2π/k
@show T = 2π/ω
ηᵢₙ(x,t) = sum( η₀ .* cos.(k*x[1] - ω*t + α) )
ϕᵢₙ(x,t) = (η₀.*ω./k) .* (cosh.(k*(H0 + x[2])) ./ 
  sinh.(k*H0)) .* sin.(k*x[1]-ω*t + α)
vᵢₙ(x,t) = sum( -(η₀.*ω) .* (cosh.(k*(H0 + x[2])) ./ 
  sinh.(k*H0)) .* cos.(k*x[1]-ω*t + α) )
vzᵢₙ(x,t) = sum( ω .* η₀ .* sin.(k*x[1]-ω*t + α) )
ηᵢₙ(t::Real) = x -> ηᵢₙ(x,t)
ϕᵢₙ(t::Real) = x -> ϕᵢₙ(x,t)
vᵢₙ(t::Real) = x -> vᵢₙ(x,t)
vzᵢₙ(t::Real) = x -> vzᵢₙ(x,t)