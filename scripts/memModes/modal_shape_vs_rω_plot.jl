module plot_modes_vs_rω
    
using JLD2, Plots, LaTeXStrings, Printf

# === 用户参数设置 ===
rM = 1000.0
rω_list = [1.5, 2.0, 2.5, 3.0]  # 用户自定义
n_modes = 6
Lm = 20.0
base_dir = "data/sims_memmodes/lrmm_modes_wet/"
output_dir = joinpath(base_dir, "plot_rω_sweep")

mkpath(output_dir)  # 自动创建输出目录

# 存储每个模态的所有曲线
modal_data = [Dict("xp" => nothing, "shapes" => []) for _ in 1:n_modes]
labels = String[]

# === 加载数据 ===
for rω in rω_list
    folder = @sprintf("lrmm_modes_rω%.2f_rmass%.2f", rω, rM)
    path = joinpath(base_dir, folder, "mem_modesdata.jld2")

    if !isfile(path)
        @warn "文件缺失: $path"
        continue
    end

    data = load(path)
    xp = data["xp"]
    Vlist = data["V"]  # Vlist[i] 是第 i 个模态的 ComplexF64[]

    for i in 1:n_modes
        V = Vlist[i]
        xp = (data["xp"] .-20) ./ 10  # 单位化横坐标
        shape = real.(V) ./ maximum(abs.(V))  # 实部归一化
        push!(modal_data[i]["shapes"], shape)
    end
    modal_data[i]["xp"] = xp 
    push!(labels, @sprintf("rω = %.2f", rω))
end

# === 绘图 ===
for i in 1:n_modes
    xp = modal_data[1]["xp"]
    shapes = modal_data[i]["shapes"]

    plt = plot(
        xlabel = L"x / L_m", ylabel = "Normalized shape",
        title = "Mode $i",
        legend = :right, lw = 2,
        xlims = (-0.05, 1.0),
        ylims = (-1.0, 1.0)  
    )

    for (j, shape) in enumerate(shapes)
        plot!(plt, xp, shape, label = labels[j])
    end

    savefile = joinpath(output_dir, "mode_$(i)_rM$(Int(rM)).png")
    savefig(savefile)
    println("✅ Saved: $savefile")
end

end # module