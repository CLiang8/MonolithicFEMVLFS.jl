module plot_frequency_scan
using JLD2, Plots, Printf

# 参数
resDir = "data/sims_mem_freq_lrmm/Frequency_scan"
rω_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]  # 与仿真一致
Kr_list, Kt_list, Ka_list, Err_list = [], [], [], []
ω_ref = nothing  # 统一ω向量

for rω in rω_list
    folder = joinpath(resDir, "rω_$(round(rω; digits=2))")
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

    if ω_ref === nothing
        ω_ref = ω
    end
end

# 转为矩阵
Kr_mat = hcat(Kr_list...)'  # size: nrω × nω
Kt_mat = hcat(Kt_list...)'
Ka_mat = hcat(Ka_list...)'
Err_mat = hcat(Err_list...)'

out_dir = joinpath(resDir, "merged")
mkpath(out_dir)

# 绘图函数
function plot_energy_map(Z, ω, rω_list; title="", savepath="")
    graymap = cgrad([RGB(0.95,0.95,0.95), RGB(0.2,0.2,0.2)])
    heatmap(
        ω, rω_list, Z,
        xlabel = "Wave frequency ω (rad/s)",
        ylabel = "Resonator rω (rad/s)",
        colorbar_title = "Value",
        title = title,
        c = cgrad([RGB(0.95,0.95,0.95), RGB(0.2,0.2,0.2)]),
        clims = (0, 1.0),
        dpi = 600
    )
    if savepath != ""
        savefig(savepath)
        println("✅ Saved: $savepath")
    end
end

# 绘图并保存
plot_energy_map(Kr_mat, ω_ref, rω_list, title="Reflection Coefficient Kr", savepath=joinpath(out_dir, "Kr_map.png"))
plot_energy_map(Kt_mat, ω_ref, rω_list, title="Transmission Coefficient Kt", savepath=joinpath(out_dir, "Kt_map.png"))
plot_energy_map(Ka_mat, ω_ref, rω_list, title="Absorption Coefficient Ka", savepath=joinpath(out_dir, "Ka_map.png"))
plot_energy_map(Err_mat, ω_ref, rω_list, title="Energy Balance Error", savepath=joinpath(out_dir, "Error_map.png"))

end # module plot_frequency_scan