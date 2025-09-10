module run_lrmm

using DrWatson
using Parameters
using WaveSpec
using Gridap
using .Constants

@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("lrmm","mem_freq_damp_lrmm_fnc.jl"))
include("plot.jl")

resDir::String = "data/sims_mem_freq_lrmm"

# # Warm-up run
# params = Memb_undamped_2D.Memb_params_warmup(name = resDir*"/Coefficients")
# Memb_undamped_2D.main(params)

# single rω run for coefficients and contour
params = Memb_undamped_2D.Memb_params(name = resDir*"/Coefficients")
Memb_undamped_2D.main(params)

include("plot.jl")
using .plot_response_contour
using .plot_energy_coefficients

name = "data/sims_mem_freq_lrmm/Coefficients/"
if isfile(joinpath(name, "mem_data.jld2"))
    plot_contour(name)
    plot_coefficients(name)
else
    @warn "Missing data file at $name"
end

params = Memb_undamped_2D.Memb_params(name = resDir*"/Memcase")
Memb_undamped_2D.main(params)

include("plot.jl")
using .plot_response_contour
using .plot_energy_coefficients

name = "data/sims_mem_freq_lrmm/Memcase/"
if isfile(joinpath(name, "mem_data.jld2"))
    plot_contour(name)
    plot_coefficients(name)
else
    @warn "Missing data file at $name"
end



# production run
@with_kw struct run_params
  name = resDir*"/Damping_3hz"
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [2*π/2.53079486745378, 2*π/2.0]
  # η₀ = [0.25, 0.25]
  ω = 0.7:0.05:5.0 
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)
  # k = dispersionRelAng.(H0, ω; msg=false)

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff

  #oscillator parameters
  rM = 1.0e3 #Kg
  rω = 3.0 #rad/s
  rK = rM*rω^2 #N/m
  ζ = 0.01 # 0.05 #damping ratio
  rC = 2*ζ*sqrt(rK*rM) #N*s/m

  #DiracDelta
  # pts = [Point(95.0, 0.0)]
  xr::Float64 = 90.0

  # Domain 
  nx = 1650
  ny = 20
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  Ld = 15*H0 #damping zone length
  LΩ = 18*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 8*H0
  xm₁ = xm₀ + Lm

  # Probes
  prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
          52.7, 53.7, 55, 60.0, 80.0, 
          85.0, 90.0, 95.0, 100.0, 120.0, 
          125.0, 140.0, 160.0, 180.0 ]
  prbPowx=[ 55.0, 125.0 ]

end
params = run_params()
# Memb_undamped_2D.main(params)

# plot_contour(params.name)
# plot_coefficients(params.name)

# ------------ rω sweep (natural frequency of spring) -----------------
rω_list = 3.65:0.05:5.0  # rad/s  运行时中间3.65卡了一次
# rω_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]  # rad/s
for rω in rω_list
    rM = 1.0e3  # 固定质量
    rK = rM * rω^2
    ζ = 0.0
    rC = 2 * ζ * sqrt(rK * rM)

    name = joinpath(resDir, "Frequency_scan", "rω_$(round(rω; digits=2))")
    mkpath(name)

    # 构造一个参数结构体，覆写 rM, rK, rC, name
    params = run_params(;
        name = name,
        rM = rM,
        rω = rω,
        rK = rK,
        ζ = ζ,
        rC = rC
    )

    @show params
    println("▶ Running for rω = $rω rad/s ...")
    Memb_undamped_2D.main(params)
end

end

# ---------------------- Post Processing Plot ----------------------
include("plot.jl")
using .plot_response_contour
using .plot_energy_coefficients

rω_list = 0.7:0.05:5.0  # rad/s

for rω in rω_list
    name = "data/sims_mem_freq_lrmm/Frequency_scan/rω_$(round(rω; digits=2))"
    if isfile(joinpath(name, "mem_data.jld2"))
        plot_contour(name)
        plot_coefficients(name)
    else
        @warn "Missing data file at $name"
    end
end
