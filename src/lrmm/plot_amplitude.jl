module PlotAmplitude

using JLD2
using Plots
using LaTeXStrings
using Measures

# === 参数设定 ===
Lm = 20.0
amp = 0.1    # 入射波幅值 η₀
setlw = 5    # 线宽

# === 读取数据 ===
filename = "data/sims_202507/mono_freq_lrmm_damp/lrmm_plotdata.jld2"
# filename = "data/sims_202507/mono_freq_lrmm_ndp/lrmm_plotdata.jld2"
# filename = "data/sims_202507/mono_freq_free/lrmm_plotdata.jld2"
data = load(filename)

x_η   = data["x_η"]
x_kr  = data["x_kr"]
x_kin = data["x_kin"]
x_kh  = data["x_kh"]

η_vals   = data["η_vals"]
kr_vals  = data["kr_vals"]
kin_vals = data["kin_vals"]
kh_vals  = data["kh_vals"]

# === 绘图 ===
plt1 = plot()

# Membrane deformation amplitude
plot!(plt1, x_η ./ Lm, abs.(η_vals) ./ amp,
    lc = :red, lw = setlw, label = false)

# Reflected wave amplitude
plot!(plt1, x_kr ./ Lm, abs.(kr_vals) ./ amp,
    lc = :blue, lw = setlw, label = false)

# Incoming wave amplitude
n1 = sum(x_kin .<= 80.0)
plot!(plt1, x_kin[1:n1] ./ Lm, abs.(kin_vals[1:n1]) ./ amp,
    lc = :blue, lw = setlw, ls = :dot, label = false)

plot!(plt1, x_kin[n1+1:end] ./ Lm, abs.(kin_vals[n1+1:end]) ./ amp,
    lc = :blue, lw = setlw, ls = :dot, label = false)

# free-surface wave amplitude
n2 = sum(x_kh .<= 80.0)
plot!(plt1, x_kh[1:n2] ./ Lm, abs.(kh_vals[1:n2]) ./ amp,
    lc = :black, lw = setlw, label = false)

plot!(plt1, x_kh[n2+1:end] ./ Lm, abs.(kh_vals[n2+1:end]) ./ amp,
    lc = :black, lw = setlw, label = false)

# label
lgndlw = 2
xlbl = 0:1:1
ylbl = xlbl .* 0

plot!(plt1, xlbl, ylbl, lc=:black, lw=lgndlw, label=L"|\kappa|")
plot!(plt1, xlbl, ylbl, lc=:red, lw=lgndlw, label=L"|\eta|")
plot!(plt1, xlbl, ylbl, lc=:blue, lw=lgndlw, label=L"|\kappa_r|")
plot!(plt1, xlbl, ylbl, lc=:blue, lw=lgndlw, ls=:dot, label=L"|\kappa_{in}|")
# plot!(plt1, xlbl, ylbl, lc=:green, lw=lgndlw, ls=:dash, label=L"|")  # 模态线保留空位

# === 坐标、标签、美化 ===
plot!(plt1,
    xlim = (2.5, 6.5),
    ylim = (0, 2),
    xticks = 2.5:0.5:6.5,
    yticks = 0.0:0.5:2.0,
    xlabel = L"x/L_m",
    ylabel = L"|\eta|/\kappa_0",
    title = L" rM = 1000 \, | \, rω = 2.4rad/s \, | \, ω = 2.4rad/s \,| \, ζ = 0.1",
    # title = L" Membrane-only",
    titlefontsize = 30,
    tickfontsize = 25,
    labelfontsize = 30,
    legendfontsize = 25,
    left_margin = 16mm, right_margin = 10mm,
    bottom_margin = 18mm, top_margin = 10mm,
    legend = :topright, legendcolumns = 4,
    dpi = 150, size = (2400, 700),
    grid = :true, gridcolor = :black,
    gridalpha = 0.5, gridlinestyle = :dot,
    gridlinewidth = 3
)

# === 辅助线 ===
vline!(plt1, [2.5], lw=3, color=:black, label=false)
hline!(plt1, [0], lw=3, color=:black, label=false)

# === 保存图像 ===
savefig(plt1, replace(filename, ".jld2" => "_ParaPlot.png"))
println("✅ Saved plot to ", replace(filename, ".jld2" => "_ParaPlot.png"))

end