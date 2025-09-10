module plot_q_participation

using JLD2
using Plots
using LaTeXStrings
using Printf
# 输入参数
rω_list = [1.01, 3.0, 5.0]
rM = 1000.0
base_dir = "data/sims_memmodes/lrmm_modes_wet/plot_rω_sweep"
isdir(base_dir) || mkpath(base_dir)

# 循环绘图
for rω in rω_list
    case_name = @sprintf("rω%.2f_rmass%.2f", rω, rM)
    file_path = "data/sims_memmodes/lrmm_modes_wet/lrmm_modes_" * case_name * "/mem_modesdata.jld2"
    data = load(file_path)

    q_modes = data["q_modes"]
    q_abs = [real(q[1]) for q in q_modes]
    q_abs = abs.(q_abs)  # 展平为一维向量  # 可选用 real.(q_modes[:]) 画实部
    modes = 1:length(q_abs)

    plt = bar(
        modes, q_abs,
        xlabel = "Mode number",
        ylabel = L"|q|",
        legend = false,
        title = "Modal participation of resonator q",
        bar_width = 0.6,
        linewidth = 0,
        color = :dodgerblue,
        dpi = 300,
        size = (600, 400)
    )

    savefig(plt, base_dir * "/q_bar_rω_" * @sprintf("%.2f", rω) * ".png")
    println("✅ Saved: $savefile")
end

end # module

# # see V data complex 
# base_dir = "data/sims_memmodes/lrmm_modes_wet/plot_rω_sweep"
# case_name = @sprintf("rω%.2f_rmass%.2f", 1.01, 1000.0)
# file_path = "data/sims_memmodes/lrmm_modes_wet/lrmm_modes_" * case_name * "/mem_modesdata.jld2"
# println("🔍 Loading: ", file_path)
# data = load(file_path)
# v1 = data["V"][3]