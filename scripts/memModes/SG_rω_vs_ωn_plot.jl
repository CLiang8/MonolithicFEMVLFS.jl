module plot_lrmm_shagun_vs_rω

using JLD2
using Plots
using Printf
using LaTeXStrings
using DataFrames
using MonolithicFEMVLFS.Resonator
using Main.Membrane_modes

# ===================== 配置参数 =====================
base_dir = "data/sims_202508/mem_modes_lrmm_SG/"
output_dir = joinpath(base_dir, "plot")
rM_list = [1.00]
# rω_list = [1.01, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50]
rω_list = [1.01, 1.50]
n_modes = 6

# ===================== 主函数 =====================
function plot_modes_for_rM(rM::Float64)
    println("📊 Plotting for rM = $rM ...")

    ωₙ_matrix = []  # 存储每个 rω 下的 ωₙ 向量
    valid_rω = []   # 实际成功加载的 rω

    for rω in rω_list
        folder = @sprintf("rω%.2f_rmass%.2f", 1.01, 1.00)
        path = joinpath(base_dir, folder, "mem_modesdata.jld2")
        
        # if isfile(path)
        #     data = load(path)

        #     # ✅ 支持 dfWet 结构
        #     if haskey(data, "dfWet")
        #         dfWet = data["dfWet"]

        #         if :ωn ∈ propertynames(dfWet)
        #             ωₙ = dfWet.ωn
        #             @assert length(ωₙ) == n_modes "⚠️ number of modes do not match at $(path)"
        #             push!(ωₙ_matrix, ωₙ)
        #             push!(valid_rω, rω)
        #         else
        #             @warn "⛔ 字段 :ωn 不存在 in dfWet. 实际字段有: $(propertynames(dfWet))"
        #         end
        #     else
        #     @warn "⛔ 文件中没有 dfWet in: $path"
        #     end
        # else
        #     @warn "⛔ 缺失文件: $path"
        # end

        # 直接加载数据
        data = load(path)
        dfWet = data["dfWet"]
        ωₙ = dfWet.ωn
        @assert length(ωₙ) == n_modes "⚠️ number of modes do not match at $(path)"
        push!(ωₙ_matrix, ωₙ)
        push!(valid_rω, rω)
    end

    if isempty(valid_rω)
        @warn "⚠ 没有找到 rM = $rM 的任何有效数据"
        return
    end

    # 转换为矩阵（每行是一个 rω 的模态频率）
    ωₙ_mat = hcat(ωₙ_matrix...)'  # size = (num_rω × n_modes)

    # ========= 绘图 =========
    plt = plot()
    for i in 1:n_modes
        plot!(valid_rω, real.(ωₙ_mat[:, i]),
              label = "mode $(i-1)",
              lw = 2,
              marker = :circle,
              markersize = 4)
    end

    xlabel!(L"\omega_r")
    ylabel!(L"\omega_n^{(i)}")
    title!("Wet natural frequencies vs resonator frequency\n(rM = $(Int(rM)))")
    plot!(legend = :right, grid = true)

    savefile = joinpath(output_dir, @sprintf("ωn_vs_rω_rM%d.png", Int(rM)))
    mkpath(output_dir)
    savefig(savefile)
    println("✅ 已保存图像: $savefile")
end

# ===================== 循环绘图 =====================
for rM in rM_list
    plot_modes_for_rM(rM)
end

end  # module
