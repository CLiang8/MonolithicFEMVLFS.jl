
using JLD2, Plots, LaTeXStrings, Measures, Dierckx, QuadGK

# === 你可以在这里自定义对比数据路径 ===
data_paths = Dict(
    "With LRMM"    => "data/sims_202509/newcorefunc/mem_data.jld2",
    "Without LRMM" => "data/sims_202509/onlymem/mem_data.jld2"
)

save_dir = raw"data\sims_202509\newcorefunc\mem_plots"

# === 加载每个案例的数据 ===
results = Dict{String, Dict}()

for (label, path) in data_paths
    data = load(path)
    results[label] = Dict(
        :ω   => data["ω"],
        :Prf => data["prbPow"][:,2],
        :Pin => data["prbPow"][:,1],
        :Ptr => data["prbPow"][:,3],
        :Pd  => data["prbPow"][:,4]
    )
end

# === 绘图：Prf vs ω 对比 ===
plt_prf = plot()

for (label, r) in results
    plot!(plt_prf, r[:ω], r[:Prf], label=label, lw=2)
end

xlabel!(plt_prf, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt_prf, L"P_{rf}\ \mathrm{(W/m)}")
title!(plt_prf, "Reflected Power Comparison")
savefig(plt_prf, joinpath(save_dir, "Prf_comparison.png"))
println("✅ Saved: Prf_comparison.png")

# === 计算积分：∫ Prf(ω) dω ===
println("\n===== Reflected Power Integral Summary =====")
for (label, r) in results
    ω = r[:ω]
    Prf = r[:Prf]
    spline = Spline1D(ω, Prf)
    area, err = quadgk(spline, minimum(ω), maximum(ω))
    println("[$label] ∫ Prf(ω) dω ≈ ", round(area, digits=5), " W·rad/(m·s)  (±", round(err, digits=2), ")")
end

# === Power Balance subplot 对比 ===
plt1 = plot(); plt2 = plot(); plt3 = plot(); plt4 = plot()

for (label, r) in results
    ω = r[:ω]
    Pin, Prf, Ptr, Pd = r[:Pin], r[:Prf], r[:Ptr], r[:Pd]

    Kr = Prf ./ Pin
    Kt = Ptr ./ Pin
    Ka = Pd  ./ Pin
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
    layout=(2,2), size=(1000,700),
    plot_title="Power Balance Comparison",
    margin=3mm)

savefig(fullplot, joinpath(save_dir, "power_balance_comparison.png"))
println("✅ Saved: power_balance_comparison.png")

