module plot_modes_vs_rM

using JLD2, Plots, LaTeXStrings, Printf

# === 用户参数设置 ===
rω = 2.40
rM_list = [0.3, 0.5, 1.0]  # ✅ 改为循环 rM
n_modes = 6
Lm = 20.0
base_dir = "data/sims_memmodes/lrmm_modes_wet/"
output_dir = joinpath(base_dir, "plot_rM_sweep")

mkpath(output_dir)  # 自动创建输出目录

# 存储每个模态的所有曲线
modal_data = [Dict("xp" => nothing, "shapes" => []) for _ in 1:n_modes]
labels = String[]

# === 加载数据 ===
for rM in rM_list
    folder = @sprintf("lrmm_modes_rω%.2f_rmass%.2f", rω, rM)
    path = joinpath(base_dir, folder, "mem_modesdata.jld2")

    if !isfile(path)
        @warn "文件缺失: $path"
        continue
    end

    data = load(path)
    xp = (data["xp"] .- Lm) ./ 10  # ✅ 单位化横坐标（可调整偏移）
    Vlist = data["V"][2:end]  # Skip the first mode (rigid body mode)

    for i in 1:n_modes
        V = Vlist[i]
        shape = real.(V) ./ maximum(abs.(V))
        push!(modal_data[i]["shapes"], shape)
    end

    modal_data[1]["xp"] = xp
    push!(labels, @sprintf("rM = %.2f", rM))
end

# === 绘图 ===
for i in 1:n_modes
    xp = modal_data[1]["xp"]
    shapes = modal_data[i]["shapes"]

    plt = plot(
        xlabel = L"x / L_m", ylabel = "Normalized shape",
        title = "Mode $i (rω = $(rω))",
        legend = :right, lw = 2,
        xlims = (-0.05, 1.0),
        ylims = (-1.0, 1.0)
    )

    for (j, shape) in enumerate(shapes)
        plot!(plt, xp, shape, label = labels[j])
    end

    savefile = joinpath(output_dir, "mode_$(i)_rω$(rω).png")
    savefig(savefile)
    println("✅ Saved: $savefile")
end

end # module