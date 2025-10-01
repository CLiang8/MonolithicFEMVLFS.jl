module plot_spatial_scan_contour
using JLD2, Gridap, DataFrames, Plots, LaTeXStrings

# 这个脚本用于绘制每个频率文件夹下merge文件夹中的图片

# 参数设置
baseDir = "data/sims_mem_freq_lrmm/Spatial_scan"
# xr_list = [80.05, 90, 99.95]
# rω_list = [1.0, 2.0]
xr_list = [80.05, 82.5, 80+10/3, 85, 80+20/3, 87.5, 90, 92.5, 90+10/3, 95, 90+20/3, 97.5, 99.95]
rω_list = [1.0, 1.55, 2.0, 2.41, 3.0, 3.47, 4.0, 4.64, 5.0 ]

# ---------- 固定读取 ω_ref ----------
first_rω_dir = joinpath(baseDir, "Spatial_scan_rω_$(round(rω_list[1]; digits=2))")
first_file = joinpath(first_rω_dir, "xr_$(round(xr_list[1]; digits=1))", "mem_data.jld2")
ω_ref = load(first_file)["ω"]

# ---------- 定义绘图函数 ----------
function plot_energy_map(Z, ω, xr_list; title="", savepath="")
    heatmap(
        (xr_list .- 80) ./ 20.0, ω, Z,
        xlabel = L"x_r / L_m",
        ylabel = L"\omega\ (rad/s)",
        xlims = (0, 1.0),
        ylims = (minimum(ω), maximum(ω)),
        c = cgrad([RGB(0.95,0.95,0.95), RGB(0.2,0.2,0.2)]),
        colorbar_title = "Value",
        clims = (0, 1.0),
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

# ---------- 主循环：每个 rω ----------
for rω in rω_list
    resDir = joinpath(baseDir, "Spatial_scan_rω_$(round(rω; digits=2))")
    Kr_list, Kt_list, Ka_list, Err_list = [], [], [], []

    for xr in xr_list
        folder = joinpath(resDir, "xr_$(round(xr; digits=1))")
        file = joinpath(folder, "mem_data.jld2")
        if !isfile(file)
            @warn "Missing data at $folder"
            continue
        end
        data = load(file)
        prbPow = Matrix(data["prbPow"])
        Pin, Prf, Ptr, Pabs = prbPow[:, 1], prbPow[:, 2], prbPow[:, 3], prbPow[:, 4]
        Kr, Kt, Ka = Prf ./ Pin, Ptr ./ Pin, Pabs ./ Pin
        Err = 1 .- (Kr .+ Kt .+ Ka)

        push!(Kr_list, Kr)
        push!(Kt_list, Kt)
        push!(Ka_list, Ka)
        push!(Err_list, Err)
    end

    Kr_mat = hcat(Kr_list...)
    Kt_mat = hcat(Kt_list...)
    Ka_mat = hcat(Ka_list...)
    Err_mat = hcat(Err_list...)

    out_dir = joinpath(resDir, "merged")
    mkpath(out_dir)

    plot_energy_map(Kr_mat, ω_ref, xr_list, title=L"K_r", savepath=joinpath(out_dir, "Kr_map.png"))
    plot_energy_map(Kt_mat, ω_ref, xr_list, title=L"K_t", savepath=joinpath(out_dir, "Kt_map.png"))
    plot_energy_map(Ka_mat, ω_ref, xr_list, title=L"K_a", savepath=joinpath(out_dir, "Ka_map.png"))
    plot_energy_map(Err_mat, ω_ref, xr_list, title="Energy Balance Error", savepath=joinpath(out_dir, "Error_map.png"))
end

end # module plot_spatial_scan_contour
