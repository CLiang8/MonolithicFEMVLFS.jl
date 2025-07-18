module run_lrmm

using DrWatson
using Parameters
using WaveSpec
using .Constants
@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("lrmm","mem_freq_damp_free_fnc.jl"))

resDir::String = "data/sims_mem_freq_lrmm/Coefficients"

# Warm-up run
params = Memb_undamped_2D.Memb_params_warmup(name = resDir)
Memb_undamped_2D.main(params)

# Production run
@with_kw struct run_params
  name = resDir
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
  ω = 0.75:0.05:5
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
  rK = 5.9e3 #N/m
  ζ = 0 # 0.05 #damping ratio
  rC = 2*ζ*sqrt(rK*rM) #N*s/m

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
Memb_undamped_2D.main(params)
end


# ---------------------- Plot Energy Coefficients ----------------------
module plot_energy_coefficients

using JLD2
using Plots
using DataFrames

println("Loading and plotting energy coefficients...")

name::String = "data/sims_mem_freq_lrmm/Coefficients"
file = name * "/mem_data.jld2"
data = load(file)
prbPow = Matrix(data["prbPow"])
ω = data["ω"]

Pin = prbPow[:, 1]
Prf = prbPow[:, 2]
Ptr = prbPow[:, 3]
Pd_total  = prbPow[:, 4]

Kr = Prf ./ Pin
Kt = Ptr ./ Pin
Ka = Pd_total  ./ Pin
Err = 1 .- (Kr .+ Kt .+ Ka)

# 反射系数 Kr
plot(ω, Kr, label="Kr (Reflected)", lw=2)
xlabel!("ω (rad/s)")
ylabel!("Kr")
title!("Reflection Coefficient vs Frequency")
savefig(joinpath(name, "Kr_vs_omega.png"))
println("✅ Saved: Kr_vs_omega.png")

# 透射系数 Kt
plot(ω, Kt, label="Kt (Transmitted)", lw=2)
xlabel!("ω (rad/s)")
ylabel!("Kt")
title!("Transmission Coefficient vs Frequency")
savefig(joinpath(name, "Kt_vs_omega.png"))
println("✅ Saved: Kt_vs_omega.png")

# 吸收系数 Ka
plot(ω, Ka, label="Ka (Absorbed)", lw=2)
xlabel!("ω (rad/s)")
ylabel!("Ka")
title!("Absorption Coefficient vs Frequency")
savefig(joinpath(name, "Ka_vs_omega.png"))
println("✅ Saved: Ka_vs_omega.png")

# 4 in 1 plot
plot(ω, Kr, label="Kr (Reflected)", lw=2)
plot!(ω, Kt, label="Kt (Transmitted)", lw=2)
plot!(ω, Ka, label="Ka (Absorbed)", lw=2)
plot!(ω, Err, label="Energy Error", lw=1, linestyle=:dash)

xlabel!("ω (rad/s)")
ylabel!("Energy Coefficient")
title!("Energy Coefficients vs Frequency")
savefig( name * "/energy_coeffs.png")
println("Plot saved to: $(name)/energy_coeffs.png")


# # 探针点响应图（随频率）
# # --------------------------
# ω = data["ω"]
# k = data["k"]
# H0 = 10.0  # 你模拟中用的 H0，用于计算 kh
# prbxy = data["prbxy"]
# prbDa = data["prbDa"]
# prbDa_x = data["prbDa_x"]
# prbPow = data["prbPow"]

#   for lprb in 1:length(prbxy)
#     xloc = prbxy[lprb][1]

#     plt1 = plot(k .* H0, abs.(prbDa[:, lprb]), lw=2,
#         xlabel = "kh", ylabel = "|A| (m)", title = "Amplitude")

#     plt2 = plot(k .* H0, abs.(prbDa_x[:, lprb]), lw=2,
#         xlabel = "kh", ylabel = "|dA/dx|", title = "Slope Magnitude")

#     plt3 = plot(k .* H0, angle.(prbDa[:, lprb]), lw=2,
#         xlabel = "kh", ylabel = "∠A (rad)", title = "Phase")

#     plt4 = plot(k .* H0, angle.(prbDa_x[:, lprb]), lw=2,
#         xlabel = "kh", ylabel = "∠dA/dx (rad)", title = "Slope Phase")

#     plt_all = plot(plt1, plt2, plt3, plt4, layout=(4,1), size=(600, 800),
#         plot_title = "Probe at x = $xloc")

#     savefig(plt_all, name * "/_dxPrb_$lprb.png")
#   end
end # module plot energy coefficients
