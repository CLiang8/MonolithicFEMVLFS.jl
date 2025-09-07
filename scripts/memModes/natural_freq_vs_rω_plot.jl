module plot_rω_vs_ωₙ

using JLD2
using Plots
using LaTeXStrings
using Printf

# ========= 用户配置区域 =========
base_dir = "data/sims_memmodes/lrmm_modes_wet/"
output_dir = joinpath(base_dir, "plot")  # 新增保存图片的目录
# rM_list = [50.0, 500.0, 1000.0]
rω_list = [1.01, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50]
rM_list = [1000.0]
# rω_list = 1.50:0.1:2.50
# rω_list = [1.01, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50 ]
n_modes = 6


# ========= 主绘图函数 =========
function plot_modes_for_rM(rM::Float64)

    println("📊 Plotting for rM = $rM ...")

    ωₙ_matrix = [] # 一维数组
    valid_rω = []

    for rω in rω_list
        # 拼接路径
        folder = @sprintf("lrmm_modes_rω%.2f_rmass%.2f", rω, rM)
        path = joinpath(base_dir, folder, "mem_modesdata.jld2")

        if isfile(path)
            data = load(path)
            ωₙ = vec(data["ωₙ"])
            @assert length(ωₙ) == n_modes "模态数量不一致 at $(path)"
            push!(ωₙ_matrix, ωₙ)
            push!(valid_rω, rω)
        else
            @warn "文件缺失: $path"
        end
    end

    if isempty(valid_rω)
        @warn "⚠ 没有找到 rM = $rM 的任何有效数据"
        return
    end

    # 转换为矩阵：每行一个 rω 对应的模态频率向量
    ωₙ_mat = hcat(ωₙ_matrix...)'  # size = (num_rω × n_modes)

    # ========= 绘图 =========
    plt = plot()
    for i in 1:n_modes
        plot!(valid_rω, ωₙ_mat[:, i],
              label = "mode $i",
              lw = 2,
              ylims = (0.0, 6.2 ),
              marker = :circle,
              markersize = 3)
    end

    xlabel!(L"\omega_r")
    ylabel!(L"\omega_n^{(i)}")
    title!("Wet natural frequencies vs resonator frequency\n(rM = $(Int(rM)))")
    plot!(legend = :topright, grid = true)

    # savefile = joinpath(output_dir, @sprintf("ωn_vs_rω_rM%d.png", Int(rM)))
    savefile = joinpath(output_dir, @sprintf("ωn_vs_rω_rM%d.png", Int(rM)))
    savefig(savefile)
    println("✅ Saved: $savefile")
end


# ========= 循环绘制所有 rM =========
for rM in rM_list
    plot_modes_for_rM(rM)
end

end