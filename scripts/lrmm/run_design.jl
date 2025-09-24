module run_lrmm

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
include(srcdir("lrmm","mem_freq_lrmm_dsn_fnc.jl"))

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

ω, η₀, α = generate_wave_parameters(0.4, 2.5, 0.95; plotloc=resDir*"/mem")

# # === Run simulation with Membrane + LRMM ===
# params = Memb2D.Memb_params(;
#     name = "data/sims_202509/newcorefunc",
#     ω = ω,
#     η₀ = η₀,
#     α = α,
#     T = 2π ./ ω,
# )

# tick()
# Memb2D.main(params)
# tock()

# === Run simulation with Membrane only (no LRMM) ===
params_onlymem = Memb2D.Memb_params_nolrmm(;
    name = "data/sims_202509/onlymem",
    ω = ω,
    η₀ = η₀,
    α = α,
    T = 2π ./ ω,
)

Memb2D.main(params_onlymem)
end



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

name::String = "data/sims_202509/newcorefunc"
file = name * "/mem_data.jld2"
# file = name * "/mem_data_merged.jld2"
data = load(file)
prbPow = Matrix(data["prbPow"])
ω = data["ω"]

Pin = prbPow[:, 1]
Prf = prbPow[:, 2]
Ptr = prbPow[:, 3]
Pd_total  = prbPow[:, 4] + prbPow[:, 7]

Kr = Prf ./ Pin
Kt = Ptr ./ Pin
Ka = Pd_total  ./ Pin
Err = 1 .- (Kr .+ Kt .+ Ka)

# # 反射系数 Kr
# plot(ω, Kr, label="Kr (Reflected)", lw=2)
# xlabel!("ω (rad/s)")
# ylabel!("Kr")
# title!("Reflection Coefficient vs Frequency")
# savefig( name * "/mem_plots/Kr_vs_omega.png")
# println("✅ Saved: Kr_vs_omega.png")

# # 透射系数 Kt
# plot(ω, Kt, label="Kt (Transmitted)", lw=2)
# xlabel!("ω (rad/s)")
# ylabel!("Kt")
# title!("Transmission Coefficient vs Frequency")
# savefig( name * "/mem_plots/Kt_vs_omega.png")
# println("✅ Saved: Kt_vs_omega.png")

# # 吸收系数 Ka
# plot(ω, Ka, label="Ka (Absorbed)", lw=2)
# xlabel!("ω (rad/s)")
# ylabel!("Ka")
# title!("Absorption Coefficient vs Frequency")
# savefig( name * "/mem_plots/Ka_vs_omega.png")
# println("✅ Saved: Ka_vs_omega.png")

# 反射能量 Pr
plot(ω, Prf, label="Pr (Reflected)", lw=2)
xlabel!(L"\omega\ \mathrm{(rad/s)}")
ylabel!(L"P_{rf}\ \mathrm{(W/m)}")
title!("Reflection Power vs Frequency")
savefig( name * "/mem_plots/Pr_vs_omega.png")
println("✅ Saved: Pr_vs_omega.png")

# --- fit and numerical integration to get total reflected power ---
spline = Spline1D(ω, Prf)
Prf_fit = spline.(ω)

# # plot raw and fitted Prf
# plt = plot(ω, Prf, label="Prf (Raw)", lw=2, markersize=3)
# plot!(plt, ω, Prf_fit, label="Prf (Spline Fit)", lw=2, linestyle=:dash)
# xlabel!(plt, L"\omega\ \mathrm{(rad/s)}")
# ylabel!(plt, L"P_{rf}\ \mathrm{(W/m)}")
# title!(plt, "Reflected Power vs Frequency (Spline Fit)")
# savefig(plt, name * "/mem_plots/Pr_vs_omega_spline.png")
# println("✅ Saved: Pr_vs_omega_spline.png")


# intergration
area, err = quadgk(spline, minimum(ω), maximum(ω))
println("✅ Total reflected power (spline fit integral): ", round(area, digits=5), " W·rad/(m⋅s)")
println("Estimated error: ±", round(err, digits=2))
# --- end fit and numerical integration ---

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
    margin = 10mm )

savefig(fullplot, name * "/mem_plots/mem_powerBalance_grid.png")
println("✅ Saved: mem_powerBalance_grid.png")


end # module plot energy coefficients
