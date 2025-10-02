module run_lrmmDesign

using Gridap
using DrWatson
using Parameters
using WaveSpec
using .Constants
using .Jonswap
using MonolithicFEMVLFS.Resonator
using TickTock


@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("lrmm","mem_freq_lrmm_dsn_fnc copy.jl"))
include(srcdir("empt","copy_clean_version.jl"))

resDir::String = "data/sims_202509/newcorefunc"
filename = resDir*"/mem"

function generate_wave_parameters(Hs::Float64, Tp::Float64, cov::Float64; 
    plotloc::String = resDir, 
    nω::Int = 200, 
    seed::Int = 100)

    # === Wave parameters ===
    ω, S, η₀ = jonswap(Hs, Tp; 
        plotflag = true, 
        plotloc = plotloc, 
        nω = nω)

    println(ω[1], "\t", ω[2], "\t", ω[end])

    ω = ω[2:end]
    S = S[2:end]
    η₀ = η₀[2:end]

    # === Calculate Energy ===
    function cumtrapz(ω, y)
        n = length(ω)
        E = zeros(Float64, n)
        for i in 2:n
            Δ = ω[i] - ω[i-1]
            E[i] = E[i-1] + 0.5*(y[i] + y[i-1]) * Δ
        end
        return E
    end

    E   = cumtrapz(ω, S)
    m0  = E[end]
    El  = (1 - cov)/2 * m0
    Er  = (1 + cov)/2 * m0
    iL  = findfirst(>=(El), E)
    iR  = findfirst(>=(Er), E)

    ω_sel  = ω[iL:iR]
    S_sel  = S[iL:iR]
    η₀_sel = η₀[iL:iR]
    α      = randomPhase(ω_sel; seed=seed)

    # === plot (ω,η₀) according to max S  ===
    i_peak = argmax(S_sel)
    println("Peak of S(ω): ω = ", ω_sel[i_peak], " rad/s, η₀ = ", η₀_sel[i_peak], " m")

    @show length(ω_sel)

    return ω_sel, η₀_sel, α
end

ω, η₀, α = generate_wave_parameters(4.0, 2.5, 0.95; plotloc=resDir*"/mem")

# # ===================== Run simulation with empty tank ============================
# params = Empty2D_design.Memb_params_warmup(;
#     name = "data/sims_202509/empt",
#     ω = ω,
#     η₀ = η₀,
#     α = α,
#     T = 2π ./ ω,
# )

# Empty2D_design.main(params)

# ================ Run simulation with Membrane only (no LRMM) ====================
params_onlymem = Memb2D_design.Memb_params_nolrmm(;
    name = "data/sims_202509/onlymem",
    ω = ω,
    η₀ = η₀,
    α = α,
    T = 2π ./ ω,
)

Memb2D_design.main(params_onlymem)


# # =================== Run simulation with Membrane + LRMM =========================
# params = Memb2D_design.Memb_params(;
#     name = "data/sims_202509/newcorefunc",
#     ω = ω,
#     η₀ = η₀,
#     α = α,
#     T = 2π ./ ω,
# )

# tick()
# Memb2D_design.main(params)
# tock()


# # ================= Run simulation with Membrane + LRMM + damping ==================
# ζ_list = [0.0, 0.02, 0.05, 0.1, 0.5]

# for ζ in ζ_list
#     # ζ = 0.1
#     name = joinpath("data/sims_202509/mem_lrmm_ζ", "ζ_$(round(ζ; digits=3))")
#     mkpath(name)
#     rM = 1000
#     rK = 5900
#     rC = 2 * ζ * sqrt(rK * rM)
#     xr = 90.0 

#     params = Memb2D_design.Memb_params(;
#        name = name,
#        ω = ω,
#        η₀ = η₀,
#        α = α,
#        T = 2π ./ ω,
#        rS = Resonator.Array1D(1, rM, rK, rC, [Point(xr, 0.0)]) 
#     )

#     # @show params.rSi.C
#     println("▶ Running for ζ = $ζ ...")
#     Memb2D_design.main(params)
# end


end # RUN SIMULATION MODULE



# ---------------------- Plot Energy Coefficients ----------------------
module plot_energy_coefficients

using JLD2
using Plots
using DataFrames
using LaTeXStrings
using Measures
using Dierckx
using QuadGK

println("Loading and plotting energy coefficients...")

# name::String = "data/sims_202509/newcorefunc"
# name::String = "data/sims_202509/onlymem"
name::String = "data/sims_202509/empt"
# name::String = "data/sims_202509/tryMonopile"

file = name * "/mem_data.jld2"
folder = name * "/mem_plots"
mkpath(folder)

data = load(file)
prbPow = Matrix(data["prbPow"])
prbForce = Matrix(data["prbForce"])
ω = data["ω"]
k = data["k"]

Pin = prbPow[:, 1]
Prf = prbPow[:, 2]
Ptr = prbPow[:, 3]
Pd_total  = prbPow[:, 4] + prbPow[:, 7]

Kr = Prf ./ Pin
Kt = Ptr ./ Pin
Ka = Pd_total  ./ Pin
Err = 1 .- (Kr .+ Kt .+ Ka)

# Fcyl_x = abs.(real.(prbForce[:, 2])) # this is wrong
Fcyl_x = abs.(prbForce[:, 2])
Fcyl_x_nd = Fcyl_x ./ (1025* 10 * k .* ω.^2)

# monopile force Fx
plt = plot(ω, Fcyl_x,
    label="Fcyl_x (Monopile Force)",
    lw=2,
    xlabel=L"\omega\ \mathrm{(rad/s)}",
    ylabel=L"F_{x}\ \mathrm{(N/m)}",
    title="Monopile Force vs Frequency",
    margin = 10mm,
    size=(1200, 700),
    dpi = 300
)

savefig(plt, folder * "/Fcylx_vs_omega.png")

# --- fit and numerical integration to get total force on monopile ---
spline = Spline1D(ω, Fcyl_x)
Fcylx_fit = spline.(ω)

# intergration
area, err = quadgk(spline, minimum(ω), maximum(ω))
println("✅ Total monopile force (spline fit integral): ", round(area, digits=5), " N·rad/(m·s)")
println("Estimated error: ±", round(err, digits=2))

# # plot raw and fitted Prf/Fcyl_x to check if it fits corrctly
# plt = plot(ω, Fcyl_x, label="Prf (Raw)", lw=2, markersize=3)
# plot!(plt, ω, Fcylx_fit, label="Prf (Spline Fit)", lw=2, linestyle=:dash)
# xlabel!(plt, L"\omega\ \mathrm{(rad/s)}")
# ylabel!(plt, L"P_{rf}\ \mathrm{(W/m)}")
# title!(plt, "Reflected Power vs Frequency (Spline Fit)")
# savefig(plt, name * "/mem_plots/Pr_vs_omega_spline.png")
# println("✅ Saved: Pr_vs_omega_spline.png")


# --- Multi-subplot: Power Balance (Kr, Kt, Ka, Error) ---
plt1 = plot(ω, Kr, lw=2, label="", legend= false)
xlabel!(plt1, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt1, L"K_R")
title!(plt1, "Reflection coefficient")

plt2 = plot(ω, Kt, lw=2, label="", legend= false)
xlabel!(plt2, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt2, L"K_T")
title!(plt2, "Transmission coefficient")

plt3 = plot(ω, Ka, lw=2, label="", legend= false)
xlabel!(plt3, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt3, L"K_A")
title!(plt3, "Absorption coefficient")

plt4 = plot(ω, Err, lw=2, label="", legend= false)
xlabel!(plt4, L"\omega\ \mathrm{(rad/s)}")
ylabel!(plt4, L"\mathrm{Error}\ \%")
title!(plt4, "Power Relative Error")
ylims!(plt4, (-0.05, 1.05))

fullplot = plot(plt1, plt2, plt3, plt4, layout=(2,2), size=(1200,700),
    plot_title="Power Balance",
    leftmargin = 10mm,
    topmargin = 4mm )

savefig(fullplot, name * "/mem_plots/mem_powerBalance_grid.png")
println("✅ Saved: mem_powerBalance_grid.png")


end # module plot energy coefficients


#### 新增：绘制一系列 damping情况下Fx对比图 （运行一次之后就没用了）
module plot_design_results

using JLD2, Plots, LaTeXStrings, Measures, Printf
using QuadGK, Dierckx

# === 设定仿真路径 & 参数 ===
ζ_list = [0.0, 0.02, 0.05, 0.1, 0.5]
base_dir = "data/sims_202509/mem_lrmm_ζ"
Dc, Lc, ρ = 1.0, 10.0, 1025.0  # non-dimensional force parameters

# === 储存数据结构 ===
results = Dict{Float64, Dict}()

# === 读取数据 ===
for ζ in ζ_list
    path = joinpath(base_dir, "ζ_$(round(ζ; digits=3))", "mem_data.jld2")
    if isfile(path)
        data = load(path)
        ω = data["ω"]
        k = data["k"]
        prbPow = Matrix(data["prbPow"])
        prbForce = Matrix(data["prbForce"])

        Pin = prbPow[:, 1]
        Prf = prbPow[:, 2]
        Ptr = prbPow[:, 3]
        Pd  = prbPow[:, 4]
        Pdm = prbPow[:, 7]

        Fx = abs.(prbForce[:, 2])
        Fx_nd = abs.(Fx) ./ (ρ .* Dc^2 .* Lc .* k .* ω.^2)

        results[ζ] = Dict(
            :ω => ω,
            :Fx_nd => Fx_nd,
            :Kr => Prf ./ Pin,
            :Kt => Ptr ./ Pin,
            :Ka => (Pd .+ Pdm) ./ Pin,
        )
    else
        @warn "Result file not found for ζ = $ζ"
    end
end

# === 绘图模块1：Fx_nd vs ω 所有 ζ 一张图 ===
function plot_Fx_vs_omega(results::Dict)
    plt = plot()
    for (ζ, r) in sort(collect(results); by=first)
        plot!(plt, r[:ω], r[:Fx_nd],
            label = "ζ = $(ζ)",
            lw = 2)
    end
    xlabel!(plt, L"\omega\ \mathrm{(rad/s)}")
    ylabel!(plt, L"\frac{F_x}{\rho D_c^2 L_c k \omega^2}")
    title!(plt, "Non-dimensional Monopile Force")
    savefig(plt, base_dir * "/Fx_nd_comparison.png")
    println("✅ Saved: Fx_nd_comparison.png")
end

# === 绘图模块2：每个 ζ 单独画 Power Balance 图 ===
function plot_power_balance_per_case(results::Dict)
    for (ζ, r) in results
        ω, Kr, Kt, Ka = r[:ω], r[:Kr], r[:Kt], r[:Ka]
        Err = 1 .- (Kr .+ Kt .+ Ka)

        plt1 = plot(ω, Kr, lw=2, legend=false, title="Reflection", ylabel=L"K_R")
        plt2 = plot(ω, Kt, lw=2, legend=false, title="Transmission", ylabel=L"K_T")
        plt3 = plot(ω, Ka, lw=2, legend=false, title="Absorption", ylabel=L"K_A")
        plt4 = plot(ω, Err, lw=2, legend=false, title="Error", ylabel=L"\mathrm{Err}", ylims=(-0.05, 0.05))

        full = plot(plt1, plt2, plt3, plt4, layout=(2,2),
            size=(1200,700),
            plot_title = @sprintf("Power Balance (ζ = %.2f)", ζ),
            leftmargin = 10mm, topmargin = 4mm)

        savefig(full, base_dir * "/ζ_$(round(ζ; digits=3))/power_balance_grid.png")
        println("✅ Saved: power_balance_grid (ζ = $ζ)")
    end
end

# === 主调用 ===
plot_Fx_vs_omega(results)
plot_power_balance_per_case(results)

end # module

