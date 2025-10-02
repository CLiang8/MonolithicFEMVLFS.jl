module PlotDesign 

# 绘制 Prf vs ω 对比图，并计算积分； 还有一个power balance 对比图但仅作为参考

using JLD2, Plots, Gridap, DataFrames, MonolithicFEMVLFS.Resonator, LaTeXStrings, Measures, Dierckx, QuadGK

# constant parameters
ρw = 1025 
Lc = 10     # height

# === 你可以在这里自定义对比数据路径 ===
data_paths = Dict(
    # "No MEMMB or LRMM" => "data/sims_202509/empt/mem_data.jld2",
    "MEMB Without LRMM" => "data/sims_202509/onlymem/mem_data.jld2",
    "With LRMM ζ=0"     => "data/sims_202509/newcorefunc/mem_data.jld2",
    "With LRMM ζ=0.02"  => "data/sims_202509/mem_lrmm_ζ/ζ_0.02/mem_data.jld2",
    "With LRMM ζ=0.05"  => "data/sims_202509/mem_lrmm_ζ/ζ_0.05/mem_data.jld2",
    "With LRMM ζ=0.1"  => "data/sims_202509/mem_lrmm_ζ/ζ_0.1/mem_data.jld2",
    # "try monopile with lrmm" => "data/sims_202509/tryMonopile/mem_data.jld2",
)

# save_dir = raw"data\sims_202509\newcorefunc\mem_plots"
# save_dir = raw"data\sims_202509\onlymem\mem_plots"
# save_dir = raw"data\sims_202509\tryMonopile\mem_plots"
save_dir = raw"data\sims_202509\Fx_comparePlots"

# === 加载每个案例的数据 ===
results = Dict{String, Dict}()

for (label, path) in data_paths
    data = load(path)
    ω   = data["ω"]
    k   = data["k"]
    Fx  = abs.(data["prbForce"][:,2])   # x 方向力的实部
    Fx_dimless = Fx ./ (ρw.* Lc .* k .* ω.^2)  # 非维化表达式

    results[label] = Dict(
        :ω    => ω,
        :Prf  => data["prbPow"][:,2],
        :Pin  => data["prbPow"][:,1],
        :Ptr  => data["prbPow"][:,3],
        :Pd   => data["prbPow"][:,4],
        :Pd_total => data["prbPow"][:,4] + data["prbPow"][:,7],
        :Fcyl_x   => abs.(Fx),
        :Fx_dimless => Fx_dimless,
    )
end

# 自定义每条曲线的颜色,线型,顺序
plot_styles = Dict(
    # "No MEMMB or LRMM"   => (:black, :dot),
    "MEMB Without LRMM"  => (:black, :solid),
    "With LRMM ζ=0"      => (:blue, :solid),
    "With LRMM ζ=0.02"   => (:red, :solid),
    "With LRMM ζ=0.05"   => (:blue, :dot),
    "With LRMM ζ=0.1"    => (:red, :dot),
)

plot_order = [
    # "No MEMMB or LRMM",
    "MEMB Without LRMM",
    "With LRMM ζ=0",
    "With LRMM ζ=0.02",
    "With LRMM ζ=0.05",
    "With LRMM ζ=0.1",
]

# === 绘图：非维化 Fx vs ω 对比（自定义样式） ===
plt_fx = plot()

for label in plot_order
    r = results[label]
    color, style = plot_styles[label]
    plot!(plt_fx, r[:ω], r[:Fx_dimless],
        label = label,
        color = color,
        linestyle = style,
        dpi = 300,
        lw = 1.5)
end

xlabel!(plt_fx, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt_fx, L"F_x\ /(ρw D_c^2 L_c k ω^2)")
title!(plt_fx, "Non-dimensional Monopile Force Comparison")
savefig(plt_fx, joinpath(save_dir, "Fx_dimless_comparison.png"))
println("✅ Saved: Fx_dimless_comparison.png")

# 计算积分：∫ Prf(ω) dω
println("\n===== Force on Monopile Integral Summary =====")
for label in plot_order
    r = results[label]
    Fcyl_x = r[:Fcyl_x]
    ω = r[:ω]
    spline = Spline1D(ω, Fcyl_x)
    area, err = quadgk(spline, minimum(ω), maximum(ω))
    println("[$label] ∫ Fcyl_x(ω) dω ≈ ", round(area, digits=5), " W·rad/(m·s)  (±", round(err, digits=2), ")")
end


# === 绘图：K_A vs ω 对比 ===
plt_ka = plot()

for label in plot_order
    r = results[label]
    color, style = plot_styles[label]
    Ka = r[:Pd_total] ./ r[:Pin]   # K_A = Pd_total / Pin
    plot!(plt_ka, r[:ω], Ka, 
        label = label,
        color = color,
        linestyle = style,
        dpi = 300,
        lw = 1.5)
end

xlabel!(plt_ka, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt_ka, L"K_A")
title!(plt_ka, "Absorption Coefficient Comparison")

savefig(plt_ka, joinpath(save_dir, "Ka_comparison.png"))
println("✅ Saved: Ka_comparison.png")

# 计算积分：∫ Ka(ω) dω
println("\n===== Energy Absorbed by the System Summary =====")
for label in plot_order
    r = results[label]
    ω = r[:ω]
    Ka = r[:Pd_total] ./ r[:Pin]
    spline = Spline1D(ω, Ka)
    area, err = quadgk(spline, minimum(ω), maximum(ω))
    println("[$label] ∫ Kₐ(ω) dω ≈ ", round(area, digits=5), " W·rad/(m·s)  (±", round(err, digits=2), ")")
end


# === Power Balance subplot 对比 ===
plt1 = plot(); plt2 = plot(); plt3 = plot(); plt4 = plot()

for (label, r) in results
    ω = r[:ω]
    Pin, Prf, Ptr, Pd_total = r[:Pin], r[:Prf], r[:Ptr], r[:Pd_total]

    Kr = Prf ./ Pin
    Kt = Ptr ./ Pin
    Ka = Pd_total  ./ Pin
    Err = 1 .- (Kr .+ Kt .+ Ka)

    plot!(plt1, ω, Kr, label="", legend=:top, lw=2)
    plot!(plt2, ω, Kt, label="", lw=2)
    plot!(plt3, ω, Ka, label="", lw=2)
    plot!(plt4, ω, Err, label= label,legend=:topright, lw=2)
end

xlabel!(plt1, L"\omega\ \mathrm{(rad/s)}"); ylabel!(plt1, L"K_R"); title!(plt1, "Reflection Coefficient")
xlabel!(plt2, L"\omega\ \mathrm{(rad/s)}"); ylabel!(plt2, L"K_T"); title!(plt2, "Transmission Coefficient")
xlabel!(plt3, L"\omega\ \mathrm{(rad/s)}"); ylabel!(plt3, L"K_A"); title!(plt3, "Absorption Coefficient")
xlabel!(plt4, L"\omega\ \mathrm{(rad/s)}"); ylabel!(plt4, L"\mathrm{Error}"); title!(plt4, "Power Balance Error")
ylims!(plt4, (-0.05, 1.05))

fullplot = plot(plt1, plt2, plt3, plt4,
    layout=(2,2), size=(1500,900),
    plot_title="Power Balance Comparison",
    margin=3mm)

savefig(fullplot, joinpath(save_dir, "power_balance_comparison.png"))
println("✅ Saved: power_balance_comparison.png")

end