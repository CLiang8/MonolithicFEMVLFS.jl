module run_spatial_scan
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
  ζ = 0.05 # 0.05 #damping ratio
  rC = 2*ζ*sqrt(rK*rM) #N*s/m

  #DiracDelta
  #pts = [Point(95.0, 0.0)]
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

ζ_list = [0.0, 0.02, 0.05, 0.1, 0.5]

#--------------------------- Run damping scan ---------------------------
for ζ in ζ_list
    rω = 2.41 #3 # rad/s  
    rM = 1.0e3
    rK = rM * rω^2
    rC = 2 * ζ * sqrt(rK * rM)
    # xr = 90.0 

    name = joinpath(resDir, "Damping_scan", "ζ_$(round(ζ; digits=3))")
    mkpath(name)

    params = run_params(;
        name = name,
        rM = rM,
        rω = rω,
        rK = rK,
        ζ = ζ,
        rC = rC,
        # xr = xr  
    )

    @show params.ζ
    println("▶ Running for ζ = $ζ ...")
    Memb_undamped_2D.main(params)
end

end

# -------------------------- Post Processing Plot -------------------------
include("plot.jl")
using .plot_response_contour
using .plot_energy_coefficients

# damping ratio list
ζ_list = [0.0, 0.02, 0.05, 0.1, 0.5]

for ζ in ζ_list
    name = "data/sims_mem_freq_lrmm/Damping_scan/ζ_$(round(ζ; digits=3))"
    if isfile(joinpath(name, "mem_data.jld2"))
        plot_contour(name)          # 会自动从 name 中识别 xr 并绘制箭头
        plot_coefficients(name)
    else
        @warn "Missing data file at $name"
    end
end

# --------------------------------Plot damping scan heatmap----------------------------------
module plot_damping_scan
using JLD2, Plots, DataFrames, LaTeXStrings

# parameters
resDir = "data/sims_mem_freq_lrmm/Damping_scan"
ζ_list = [0.0, 0.02, 0.05, 0.1, 0.5]
Kr_list, Kt_list, Ka_list, Err_list = [], [], [], []

# extract ω from the first file
first_file = joinpath(resDir, "ζ_$(round(ζ_list[1]; digits=3))", "mem_data.jld2")
ω_ref = load(first_file)["ω"]

# load data for each ζ
for ζ in ζ_list
    folder = joinpath(resDir, "ζ_$(round(ζ; digits=3))")
    file = joinpath(folder, "mem_data.jld2")
    if !isfile(file)
        @warn "Missing data at $folder"
        continue
    end
    data = load(file)
    ω = data["ω"]
    prbPow = Matrix(data["prbPow"])

    Pin = prbPow[:, 1]
    Prf = prbPow[:, 2]
    Ptr = prbPow[:, 3]
    Pabs = prbPow[:, 4]

    Kr = Prf ./ Pin
    Kt = Ptr ./ Pin
    Ka = Pabs ./ Pin
    Err = 1 .- (Kr .+ Kt .+ Ka)

    push!(Kr_list, Kr)
    push!(Kt_list, Kt)
    push!(Ka_list, Ka)
    push!(Err_list, Err)
end

# reshape to Matrix（n_ζ * n_ω）
Kr_mat = hcat(Kr_list...)'
Kt_mat = hcat(Kt_list...)'
Ka_mat = hcat(Ka_list...)'
Err_mat = hcat(Err_list...)'

out_dir = joinpath(resDir, "merged")
mkpath(out_dir)

# plot function
function plot_energy_map(Z, ω, ζ_list; title="", savepath="")
    heatmap(
        ω, ζ_list, Z,
        xlabel = L"wave frequency $\omega\ (rad/s)$",
        ylabel = L"Damping ratio $\zeta$",
        xlims = (minimum(ω), maximum(ω)),
        ylims = (minimum(ζ_list), maximum(ζ_list)),
        c = cgrad([RGB(0.95,0.95,0.95), RGB(0.2,0.2,0.2)]),
        colorbar_title = "Value",
        clims = (0, 1.0),
        dpi = 600,
        title = title
    )
    # add dash line
    ω_ticks = range(minimum(ω), stop=maximum(ω), length=5)
    vline!(ω_ticks, lw=1, linestyle=:dot, color=:white, label="")
    hline!([0.02, 0.05, 0.1, 0.5], lw=1, linestyle=:dot, color=:white, label="")

    if savepath != ""
        savefig(savepath)
        println("✅ Saved: $savepath")
    end
end

# plot output
plot_energy_map(Kr_mat, ω_ref, ζ_list, title=L"K_r", savepath=joinpath(out_dir, "Kr_map.png"))
plot_energy_map(Kt_mat, ω_ref, ζ_list, title=L"K_t", savepath=joinpath(out_dir, "Kt_map.png"))
plot_energy_map(Ka_mat, ω_ref, ζ_list, title=L"K_a", savepath=joinpath(out_dir, "Ka_map.png"))
plot_energy_map(Err_mat, ω_ref, ζ_list, title="Energy Balance Error", savepath=joinpath(out_dir, "Error_map.png"))

end # module plot_damping_scan
