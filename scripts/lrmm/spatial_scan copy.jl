module run_spatial_scan
using DrWatson
using Parameters
using WaveSpec
using Gridap
using .Constants

@quickactivate "MonolithicFEMVLFS.jl"

# 这个文件是针对一个 rω 情况下的研究，（example case） 包含绘制每次运行（每个位置下）的结构响应contour，和 Kr，Kt vs ω 图
# 另外一个文件是 rω，xr双循环

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
  ω = 1.4:0.05:3.0

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
  ζ = 0 # 0.05 #damping ratio
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

resDir = "data/sims_mem_freq_lrmm"
xr_list = 80.1:0.1:90.1 # try increase precision，利用对称性只关注一半
# xr_list = [80.05, 82.5, 80+10/3, 85, 80+20/3, 87.5, 90, 92.5, 90+10/3, 95, 90+20/3, 97.5, 99.95]

#--------------------------- Run spatial scan ---------------------------
for xr in xr_list
    rω = 2.41 # rad/s  
    rM = 1.0e3
    rK = rM * rω^2
    ζ = 0.0
    rC = 2 * ζ * sqrt(rK * rM)

    name = joinpath(resDir, "Spatial_scan_2.41", "xr_$(round(xr; digits=1))")
    mkpath(name)

    params = run_params(;
        name = name,
        rM = rM,
        rω = rω,
        rK = rK,
        ζ = ζ,
        rC = rC,
        xr = xr  
    )

    @show params.xr
    println("▶ Running for x_r = $xr m ...")
    Memb_undamped_2D.main(params)
end

end

# -------------------------- Post Processing Plot -------------------------
include("plot.jl")
using .plot_response_contour
using .plot_energy_coefficients

# 振子位置列表（单位：米）
# xr_list = [80.05, 82.5, 80+10/3, 85, 80+20/3, 87.5, 90, 92.5, 90+10/3, 95, 90+20/3, 97.5, 99.95]
xr_list = 80.1:0.1:90.1

for xr in xr_list
    name = "data/sims_mem_freq_lrmm/Spatial_scan_2.41/xr_$(round(xr; digits=1))"
    if isfile(joinpath(name, "mem_data.jld2"))
        plot_contour(name)          # 会自动从 name 中识别 xr 并绘制箭头
        plot_coefficients(name)
    else
        @warn "Missing data file at $name"
    end
end

# --------------------------------Plot spatial scan heatmap----------------------------------
module plot_spatial_scan
using JLD2, Plots, DataFrames, LaTeXStrings

# 参数
resDir = "data/sims_mem_freq_lrmm/Spatial_scan_2.41"
# xr_list = [80.05, 82.5, 80+10/3, 85, 80+20/3, 87.5, 90, 92.5, 90+10/3, 95, 90+20/3, 97.5, 99.95]
xr_list = 80.1:0.1:90.1
Kr_list, Kt_list, Ka_list, Err_list = [], [], [], []

# 读取第一个文件的 ω 作为统一频率参考
# first_rω_dir = joinpath(baseDir, "Spatial_scan_rω_$(round(rω_list[1]; digits=2))")
first_file = joinpath(resDir, "xr_$(round(xr_list[1]; digits=1))", "mem_data.jld2")
ω_ref = load(first_file)["ω"]

for xr in xr_list
    folder = joinpath(resDir, "xr_$(round(xr; digits=1))")
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
size(Kr_list)
println(Kr_list)
# 转为矩阵（size: n_ω*n_x）
Kr_mat = hcat(Kr_list...)
Kt_mat = hcat(Kt_list...)
Ka_mat = hcat(Ka_list...)
Err_mat = hcat(Err_list...)
# Kr_mat = transpose(Kr_mat)'

out_dir = joinpath(resDir, "merged")
mkpath(out_dir)

# 绘制单个 ωr heatmap
function plot_energy_map(Z, ω, xr_list; title="", savepath="")
    heatmap(
        (xr_list .- 80) ./ 20.0, ω, Z,  # 👈 x轴为位置，y轴为频率
        xlabel = L"L_m",                    
        ylabel = L"\omega\ (rad/s)",
        xlims = (0, 0.5),
        ylims = (minimum(ω), maximum(ω)),
        c = cgrad([RGB(0.95,0.95,0.95),RGB(0.2,0.2,0.2)]),
        colorbar_title = "Value",
        clims = (0, 1.0),
        # aspect_ratio = 2/3,
        title = title,
        
        dpi = 600
    )
    # 添加横向分割线
    ω_ticks = range(1, stop=maximum(ω), length=5)
    hline!(ω_ticks, lw=1, linestyle=:dot, color=:white, label="")
    # add dash line
    vline!([0.25,0.5,0.75], lw=1, linestyle=:dot, color=:white, label="")

    if savepath != ""
        savefig(savepath)
        println("✅ Saved: $savepath")
    end
end

# 输出图像
plot_energy_map(Kr_mat, ω_ref, xr_list, title=L"K_r", savepath=joinpath(out_dir, "Kr_map.png"))
plot_energy_map(Kt_mat, ω_ref, xr_list, title=L"K_t", savepath=joinpath(out_dir, "Kt_map.png"))
plot_energy_map(Ka_mat, ω_ref, xr_list, title=L"K_a", savepath=joinpath(out_dir, "Ka_map.png"))
plot_energy_map(Err_mat, ω_ref, xr_list, title="Energy Balance Error", savepath=joinpath(out_dir, "Error_map.png"))

end # module plot_spatial_scan