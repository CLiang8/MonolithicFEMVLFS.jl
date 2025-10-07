module run_damping_scan
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
  ω = 0.95:0.05:5.0

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
  
# rMᵨ_list = [0.01, 0.05, 0.1, 0.5, 1.0]
rMᵨ_list = [0.01]


#--------------------------- Run damping scan ---------------------------
for rMᵨ in rMᵨ_list
    ρw = 1025
    rM = rMᵨ * ρw
    rω = 2.41 #3 # rad/s  
    rK = rM * rω^2
    rC = 0
    # xr = 90.0 

    name = joinpath(resDir, "MassRatio_scan", "rMᵨ_$(round(rMᵨ; digits=2))")
    mkpath(name)

    params = run_params(;
        name = name,
        rM = rM,
        rω = rω,
        rK = rK,
        rC = rC,
        # xr = xr  
    )

    @show params.rM
    println("▶ Running for rMᵨ = $rMᵨ ...")
    Memb_undamped_2D.main(params)
end

end

# -------------------------- Post Processing Plot -------------------------
include("plot.jl")
using .plot_response_contour
using .plot_energy_coefficients

# mass ratio list
rMᵨ_list = [0.01, 0.05, 0.1, 0.5, 1.0]

for rMᵨ in rMᵨ_list
    name = "data/sims_mem_freq_lrmm/MassRatio_scan/rMᵨ_$(round(rMᵨ; digits=2))"
    if isfile(joinpath(name, "mem_data.jld2"))
        plot_contour(name)          # 会自动从 name 中识别 xr 并绘制箭头
        plot_coefficients(name)
    else
        @warn "Missing data file at $name"
    end
end

# --------------------------------Plot damping scan 2D----------------------------------

module plot_damping_scan_2D
using JLD2, Plots, DataFrames, LaTeXStrings, Dierckx, QuadGK

# parameters
resDir = "data/sims_mem_freq_lrmm/MassRatio_scan"
rMᵨ_list = [1.0, 0.5, 0.1, 0.05, 0.01]
# rMᵨ_list = [1.0, 0.1, 0.01]
Kr_list, Kt_list, Ka_list, Err_list = [], [], [], []

# extract ω from the first file
first_file = joinpath(resDir, "rMᵨ_$(round(rMᵨ_list[1]; digits=2))", "mem_data.jld2")
ω_ref = load(first_file)["ω"]

# load data for each rMᵨ
for rMᵨ in rMᵨ_list
    folder = joinpath(resDir, "rMᵨ_$(round(rMᵨ; digits=2))")
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
    Err = abs.(1 .- (Kr .+ Kt .+ Ka))

    push!(Kr_list, Kr)
    push!(Kt_list, Kt)
    push!(Ka_list, Ka)
    push!(Err_list, Err)
end

# ----------------- 绘图函数 -----------------
function plot_coeff_comparison(ω, coeff_list, rMᵨ_list, ylabel_str, title_str, savepath)
    plt = plot(
        xlabel=L"\omega\ (rad/s)",
        ylabel=ylabel_str,
        title=title_str,
        xlims = (1, maximum(ω)),
        dpi=600,
        legend=:left,
        grid=true,
        size=(800, 500)
    )

    # 手动设定颜色和线型
    colors = [:red, :blue, :red, :blue, :red]
    styles = [:solid, :solid, :dot, :dot, :dash]

    for (i, coeff) in enumerate(coeff_list)
        plot!(
            ω, coeff,
            lw=1.5,
            ylims=(-0.05,1.05),
            color=colors[i],
            linestyle=styles[i],
            label=L"rM_{\rho} = "*string(rMᵨ_list[i])
        ) 
    end

    savefig(plt, savepath)
    println("✅ Saved: $savepath")
end

# ----------------- 绘制各图 -----------------
plot_coeff_comparison(ω_ref, Kr_list, rMᵨ_list, L"K_r", "Reflection Coefficient",
    joinpath(resDir, "merged/Kr_vs_omega.png"))

plot_coeff_comparison(ω_ref, Ka_list, rMᵨ_list, L"K_a", "Absorption Coefficient",
    joinpath(resDir, "merged/Ka_vs_omega.png"))

plot_coeff_comparison(ω_ref, Kt_list, rMᵨ_list, L"K_t", "Transmission Coefficient",
    joinpath(resDir, "merged/Kt_vs_omega.png"))

plot_coeff_comparison(ω_ref, Err_list, rMᵨ_list, L"Error", "Computational Error",
    joinpath(resDir, "merged/Err_vs_omega.png"))


# ----------------- Kr 积分计算 -----------------
function compute_integrated_Kr(ω, Kr_list, rMᵨ_list)
    println("📊 Integrating fitted Kr curves...")
    areas = []

    for (i, Kr) in enumerate(Kr_list)
        # 使用三次样条插值拟合
        spline = Spline1D(ω, Kr)

        # 积分范围
        ω_min, ω_max = 1, 5

        # 使用 QuadGK 对拟合曲线积分
        area, err = quadgk(spline, ω_min, ω_max)
        push!(areas, area)

        println("rMᵨ = $(rMᵨ_list[i]): ∫Kr(ω)dω = $(round(area; digits=4)) ± $(round(err; digits=2))")
    end

    return areas
end

# 调用函数计算并返回所有积分值
Kr_areas = compute_integrated_Kr(ω_ref, Kr_list, rMᵨ_list)

# # 查看spline拟合情况
# function plot_Kr_with_spline1D(ω, Kr_list, rMᵨ_list, savepath)
#     plt = plot(
#         xlabel = L"\omega\ (rad/s)",
#         ylabel = L"K_r",
#         title = "Kr with Spline1D Fit",
#         dpi = 600,
#         legend = :bottomleft,
#         grid = true,
#         size = (900, 500)
#     )

#     colors = [:red, :blue, :green, :orange, :purple]

#     for (i, Kr) in enumerate(Kr_list)
#         plot!(ω, Kr,
#             lw = 1,
#             color = colors[i],
#             linestyle = :solid,
#             label = L"Raw,\ rM\_rho = "*string(rMᵨ_list[i])
#         )

#         # 构造 spline 拟合
#         spline = Spline1D(ω, Kr)  # 默认 cubic
#         ω_dense = range(minimum(ω), stop=maximum(ω), length=300)
#         Kr_fit = spline.(ω_dense)

#         plot!(ω_dense, Kr_fit,
#             lw = 1.5,
#             color = colors[i],
#             linestyle = :dashdot,
#             label = L"Spline,\ rM\_rho = "*string(rMᵨ_list[i])
#         )

#         # 可选：积分
#         area, err = quadgk(spline, minimum(ω), maximum(ω))
#         println("rMᵨ = $(rMᵨ_list[i]): ∫Kr(ω)dω = $(round(area; digits=4))")
#     end

#     savefig(plt, savepath)
#     println("✅ Saved: $savepath")
# end


# plot_Kr_with_spline1D(ω_ref, Kr_list, rMᵨ_list,
#     joinpath(resDir, "merged/Kr_vs_omega_spline1d.png"))

end