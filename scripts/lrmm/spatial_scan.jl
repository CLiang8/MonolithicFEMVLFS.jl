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

# --------------------------- Run spatial-frequency scan over rω & xr ---------------------------
resDir = "data/sims_mem_freq_lrmm/Spatial_scan"
# mkdir(resDir)

xr_list = [80.05, 82.5, 80+10/3, 85, 80+20/3, 87.5, 90, 92.5, 90+10/3, 95, 90+20/3, 97.5, 99.95]
rω_list = [1.0, 1.55, 2.0, 2.41, 3.0, 3.47, 4.0, 4.64, 5.0 ]

# xr_list = [80.05, 90,99.95]
# rω_list = [1.0, 2.0]

for rω in rω_list
    rM = 1.0e3
    rK = rM * rω^2
    ζ = 0.0
    rC = 2 * ζ * sqrt(rK * rM)

    # 对每个 rω 建立自己的结果子目录
    subdir = "Spatial_scan_rω_$(round(rω; digits=2))"

    for xr in xr_list
        name = joinpath(resDir, subdir, "xr_$(round(xr; digits=1))")
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

        println("▶ Running for rω = $rω rad/s at x_r = $xr m ...")
        Memb_undamped_2D.main(params)
    end
end

end # module run_spatial_scan

# --------------------------------Plot spatial scan optimal points----------------------------------
module plot_spatial_scan
using JLD2, Plots, DataFrames, LaTeXStrings, Measures

# 设置参数
resDir = "data/sims_mem_freq_lrmm/Spatial_scan"
# rω_list = [1.0, 2.0]
# xr_list = [80.05, 90,99.95]
xr_list = [80.05, 82.5, 80+10/3, 85, 80+20/3, 87.5, 90, 92.5, 90+10/3, 95, 90+20/3, 97.5, 99.95]
rω_list = [1.0, 1.55, 2.0, 2.41, 3.0, 3.47, 4.0, 4.64, 5.0 ]

# 存储最优点结果
optimal_points_Kr = DataFrame(rω=Float64[], ω=Float64[], xr=Float64[])
optimal_points_Kt = DataFrame(rω=Float64[], ω=Float64[], xr=Float64[])
optimal_points_Ka = DataFrame(rω=Float64[], ω=Float64[], xr=Float64[])

for rω in rω_list
    Kr_list, Kt_list, Ka_list = [], [], []
    ω_ref = nothing  # 每个 rω 要统一 ω
    for xr in xr_list
        folder = joinpath(resDir, "Spatial_scan_rω_$(round(rω; digits=2))", "xr_$(round(xr; digits=1))")
        file = joinpath(folder, "mem_data.jld2")
        if !isfile(file)
            @warn "Missing data at $folder"
            continue
        end

        data = load(file)
        ω = data["ω"]
        ω_ref === nothing && (ω_ref = ω)  # 初始化统一参考 ω

        prbPow = Matrix(data["prbPow"])
        Pin, Prf, Ptr, Pabs = prbPow[:, 1], prbPow[:, 2], prbPow[:, 3], prbPow[:, 4]
        Kr, Kt, Ka = Prf ./ Pin, Ptr ./ Pin, Pabs ./ Pin

        push!(Kr_list, Kr)
        push!(Kt_list, Kt)
        push!(Ka_list, Ka)
    end

    # 构造响应矩阵（行：xr，列：ω）
    Kr_mat = hcat(Kr_list...)
    Kt_mat = hcat(Kt_list...)
    Ka_mat = hcat(Ka_list...)

    # find Kr max
    max_Kr = maximum(Kr_mat)
    for idx in findall(x -> isapprox(x, max_Kr; atol=5e-4), Kr_mat)
        row, col = Tuple(idx)
        push!(optimal_points_Kr, (rω=rω, ω=ω_ref[row], xr=xr_list[col]))
    end

    # find Kt min
    min_Kt = minimum(Kt_mat)
    for idx in findall(x -> isapprox(x, min_Kt; atol=1e-8), Kt_mat)
        row, col = Tuple(idx)
        push!(optimal_points_Kt, (rω=rω, ω=ω_ref[row], xr=xr_list[col]))
    end

    # find Ka max
    max_Ka = maximum(Ka_mat)
    for idx in findall(x -> isapprox(x, max_Ka; atol=1e-8), Ka_mat)
        row, col = Tuple(idx)
        push!(optimal_points_Ka, (rω=rω, ω=ω_ref[row], xr=xr_list[col]))
    end
end

# plot
function plot_optimal_points(df::DataFrame, title::String, savepath::String;
    marker_shape=:circle, color_palette=palette(:tab10))

    plt = plot(
    xlabel = L"L_m",
    ylabel = L"\omega\ (rad/s)",
    title = title,
    legend = :outerbottomright,
    # top_margin = 10mm,
    dpi = 600,
    xlims = (-0.05, 1.2),
    xticks = 0.0:0.1:1.0,
    ylims = (0.7, 5.2)
    )

    for (i, rω_val) in enumerate(sort(unique(df.rω)))
        subdf = filter(row -> row.rω == rω_val, df)
        scatter!(
            (subdf.xr.-80) ./ 20, subdf.ω,
            label="rω = $(rω_val)",
            marker=marker_shape,
            color=color_palette[i],
            markersize=6,
            markerstrokewidth=0
        )
    end

    # add dash line
    vline!([1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8], lw=1, linestyle=:dot, color=:blue, label="1/8")
    vline!([1/6, 1/3, 1/2, 2/3, 5/6], lw=1, linestyle=:dot, color=:red, label="1/6")

    savefig(plt, savepath)
    println("✅ Saved: $savepath")
end

outdir = joinpath(resDir, "Optimal_points")
mkpath(joinpath(resDir, "Optimal_points"))

plot_optimal_points(optimal_points_Kr, "Max Kr Optimal Points", joinpath(outdir, "optimal_Kr.png"))
plot_optimal_points(optimal_points_Kt, "Min Kt Optimal Points", joinpath(outdir, "optimal_Kt.png"); marker_shape=:utriangle)
plot_optimal_points(optimal_points_Ka, "Max Ka Optimal Points", joinpath(outdir, "optimal_Ka.png"); marker_shape=:diamond)
end # module plot_spatial_scan