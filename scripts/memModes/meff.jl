module MemModesMeff

using JLD2
using Printf
using Plots

# ==== 设置参数组 ====
# rM_list = [0.1, 0.3, 0.5]
# rω_list = [1.5, 2.4, 3.5, 4.5]
rM_list = [0.1]
rω_list = [2.4]

# 输出目录
output_dir = "data/sims_memmodes/lrmm_modes_wet/plot_meff"
mkpath(output_dir)


# ==== 双层循环遍历所有组合 ====
for rM in rM_list
    for rω in rω_list

        # 构造文件名
        caseName = "rω"* @sprintf("%0.2f", rω) *"_rmass"* @sprintf("%0.2f", rM)
        # caseName = "rω"* @sprintf("%0.2f", 2.4) *"_rmass"* @sprintf("%0.2f", 0.1) # 测试单个
        name = "data/sims_memmodes/lrmm_modes_wet/lrmm_modes_" * caseName
        filename = name * "/mem_modesdata.jld2"

        # println("✅ Loading: ", filename)
        data = load(filename)
        meff = real.(data["meff"]) # very low imaginary part due to numerical errors
        @show meff

        plt = plot(
            meff, seriestype = :scatter,
            xlabel = "Mode number",
            ylabel = "Effective mass",
            title = "Wet mode effective mass",
            ylims = (-0.25, 1.85),
            label = @sprintf("rω = %.2f, rM = %.2f", rω, rM),
            legend = :topright,
            )   

        savefile = joinpath(output_dir, @sprintf("meff_rω%.2f_rM%.2f.png", rω, rM))
        savefig(savefile)
        println("✅ Saved: $savefile")

    end
end


end