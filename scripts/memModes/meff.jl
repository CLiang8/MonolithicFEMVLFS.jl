module MemModesMeff

using JLD2
using Printf
using Plots

# 指定要读取的参数
rM = 0.1
rω = 2.40

# 构造路径和文件名
caseName = "rω"* @sprintf("%0.2f", rω) *"_rmass"* @sprintf("%0.2f", rM)
name = "data/sims_memmodes/lrmm_modes_wet/lrmm_modes_" * caseName
filename = name * "/mem_modesdata.jld2"
output_dir = "data/sims_memmodes/lrmm_modes_wet/plot_meff"

mkpath(output_dir)

# 加载文件并提取 meff
println("Loading file: ", filename)
data = load(filename)

# 提取和打印 meff
meff = real.(data["meff"])
println("\nEffective mass (meff) values:")
println(meff)


plot(meff, seriestype=:scatter, 
    xlabel="Mode number", 
    ylabel="Effective mass", 
    title="Wet mode effective mass", 
    ylims = (-0.25, 1.85),
    legend= true

)

savefile = joinpath(output_dir, @sprintf("meff_rω%.2f_rM%.2f.png", rω, rM))
savefig(savefile)
println("✅ Saved: $savefile")



end