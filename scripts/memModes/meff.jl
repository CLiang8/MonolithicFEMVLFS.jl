module MemModesMeff

using JLD2
using Printf
using Plots

# ==== 设置参数组 ====
# rM_list = [0.1, 0.3, 0.5, 1.0]
# rω_list = [1.5, 2.4]
rM_list = [1.0]
rω_list = [2.4]

# 输出目录
# output_dir = "data/sims_memmodes/lrmm_modes_wet/plot_meff"  # this is for transpose results
output_dir = "data/sims_memmodes/lrmm_modes_wet/plot_meff2" # this is for hermitian transpose results
# output_dir = "data/sims_memmodes/lrmm_modes_wet_realpart/plot_meff"   # this is for no damp (real part) results
mkpath(output_dir)

# 加载原始结构（无振子）结果
ref_data_path = "data/sims_memmodes/mem_modes_wet_free/mem_modesdata.jld2"
ref_data = load(ref_data_path)
# meff = ref_data["meff"]
meff_ref = real.(ref_data["meff"][2:end])

dry_datda_path = "data/sims_memmodes/mem_modes_dry_lrmm/mem_modesdata.jld2"
dry_data = load(dry_datda_path)
meff_dry = real.(dry_data["meff"][2:end])


# ==== 双层循环遍历所有组合 ====
for rM in rM_list
    for rω in rω_list

        # 构造文件名
        caseName = "rω"* @sprintf("%0.2f", rω) *"_rmass"* @sprintf("%0.2f", rM)
        # caseName = "rω"* @sprintf("%0.2f", 3.5) *"_rmass"* @sprintf("%0.2f", 0.1) # 测试单个
        # name = "data/sims_memmodes/lrmm_modes_wet/lrmm_modes_" * caseName  # this is for complex results with damp
        name = "data/sims_memmodes/lrmm_modes_wet_realpart/lrmm_modes_" * caseName  # this is for no dammp results (real part)
        filename = name * "/mem_modesdata.jld2"

        # println("✅ Loading: ", filename)
        data = load(filename)
        meff = data["meff"][2:end]  # skip the first mode (rigid body mode)
        meff_im = imag.(data["meff"])
        meff = real.(meff) 
        
        @show meff

        plt = plot(
            meff, seriestype = :scatter,
            xlabel = "Mode number",
            ylabel = "Effective mass",
            title = "Wet mode effective mass",
            ylims = (-1.0, 1.85),
            label = @sprintf("rω = %.2f, rM = %.2f", rω, rM),
            legend = :topright,
            ) 
            
        plot!(plt, meff_ref, linestyle = :dash, label = "Reference (no resonators)")
        plot!(plt, meff_dry, linestyle = :dash, color = :blue, label = "no water")

        savefile = joinpath(output_dir, @sprintf("meff_rω%.2f_rM%.2f.png", rω, rM))
        savefig(savefile)
        println("✅ Saved: $savefile")

    end
end

# # ==== 外层循环：按 rω 画图 ====
# for rω in rω_list
#     plt = plot(
#         xlabel = "Mode number",
#         ylabel = "Effective mass",
#         title = @sprintf("Wet mode effective mass at rω = %.2f", rω),
#         ylims = (-1.0, 1.85),
#         legend = :topright,
#         lw = 2,
#         marker = :auto,
#     )

#     for rM in rM_list
#         # 构造文件路径
#         caseName = "rω"* @sprintf("%0.2f", rω) *"_rmass"* @sprintf("%0.2f", rM)
#         name = "data/sims_memmodes/lrmm_modes_wet/lrmm_modes_" * caseName  # this is for complex results with damp
#         # name = "data/sims_memmodes/lrmm_modes_wet_realpart/lrmm_modes_" * caseName  # this is for no dammp results (real part)
#         filename = name * "/mem_modesdata.jld2"

#         # 加载数据
#         data = load(filename)
#         meff = real.(data["meff"][2:end])  # 跳过第一个模态

#         # 绘制当前质量下的曲线
#         plot!(plt, meff, label = @sprintf("rM = %.2f", rM))
#     end

#     # 添加参考线（无振子）
#     plot!(plt, meff_ref, linestyle = :dash, color = :black, label = "Reference")

#     # 保存图像
#     savefile = joinpath(output_dir, @sprintf("meff_rω%.2f.png", rω))
#     savefig(plt, savefile)
#     println("✅ Saved: $savefile")
# end


end