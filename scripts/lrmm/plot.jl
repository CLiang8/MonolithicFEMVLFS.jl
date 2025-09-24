module plot_response_contour

using Gridap
using JLD2
using Plots
using Measures
using DataFrames
using WaveSpec.Constants
using LaTeXStrings
using Colors

export plot_contour

function plot_contour(name::String)
    filename = joinpath(name, "mem_data.jld2")
    data = load(filename)
    name_tag = splitpath(name)[end]   

    ω = data["ω"]
    η₀ = data["η₀"]
    H0 = 10
    Lm = 20.0
    xκ = [p[1] for p in data["prxΓκ"]] ./ Lm
    xη = [p[1] for p in data["prxΓη"]] ./ Lm
    Zκ = abs.(Matrix(data["prbDaΓκ"])) ./ η₀  # κr and κt
    Zη = abs.(Matrix(data["prbDaΓη"])) ./ η₀  # η over membrane
    Zκ_r = abs.(Matrix(data["prbDaΓκr"]))./ η₀  # κr normalized by η₀
    # @show Zη[26, :]

    # stitching complete coordinates and data
    x_total = vcat(xκ[(xκ .>= 3.0) .& (xκ .< 4.0)], xη, xκ[(xκ .> 5.0) .& (xκ .<= 6.0)])
    Z_total = hcat(Zκ_r[:, (xκ .>= 3.0) .& (xκ .< 4.0)], Zη, Zκ[:, (xκ .> 5.0) .& (xκ .<= 6.0)])

    # Create a grid for the contour plot
    # gray10 = cgrad([
    #     RGB(0.95, 0.95, 0.95),  # very light gray
    #     RGB(0.88, 0.88, 0.88),
    #     RGB(0.81, 0.81, 0.81),
    #     RGB(0.74, 0.74, 0.74),
    #     RGB(0.67, 0.67, 0.67),
    #     RGB(0.60, 0.60, 0.60),
    #     RGB(0.48, 0.48, 0.48),
    #     RGB(0.36, 0.36, 0.36),
    #     RGB(0.26, 0.26, 0.26),
    #     RGB(0.18, 0.18, 0.18)   # deep gray, but not black
    # ])

    # contour plot
    plt = heatmap(
        x_total, ω, Z_total,
        xlabel = L"x / L_m",                    
        ylabel = L"\omega\ (rad/s)",
        c = cgrad([RGB(0.95,0.95,0.95), RGB(0.2,0.2,0.2)]),
        # c = gray10,
        clims = (0, 1.0),
        xlims = (3, 6),
        ylims = (1.0, maximum(ω)),
        # ylims = (minimum(ω), maximum(ω)),
        aspect_ratio = 2/3,
        title = "Normalized Amplitude $(name_tag)",
        titlefont = font(11),
        colorbar = true,
        colorbar_title = "values",
        framestyle = :box,
        dpi = 1000,
        top_margin = 6mm
    )

    # add dash line
    vline!([4.0, 5.0], lw=1, linestyle=:dot, color=:white, label="")

    # 添加顶部注释文字
    annotate!(3.5, maximum(ω) + 0.1, text(L"|\kappa_r| / \kappa_0", :black, 8, :center))
    annotate!(4.5, maximum(ω) + 0.1, text(L"|\eta| / \kappa_0", :black, 8, :center))
    annotate!(5.5, maximum(ω) + 0.1, text(L"|\kappa_t| / \kappa_0", :black, 8, :center))

    # （可选）添加横向频率分段线
    # ω_ticks = range(minimum(ω), stop=maximum(ω), length=5)
    ω_ticks = range(1, stop=maximum(ω), length=5)
    hline!(ω_ticks, lw=1, linestyle=:dot, color=:white, label="")

    # 添加红色箭头指示振子位置，默认90.0
    xr_match = match(r"xr_(\d+\.?\d*)", name)
    xr = xr_match !== nothing ? parse(Float64, xr_match.captures[1]) : 90.0
    xr_norm = xr / Lm
    arrow_y = minimum(ω)-0.02  # 贴近 x 轴

    scatter!([xr_norm], [arrow_y],
    markershape = :utriangle,
    markercolor = :red,
    markersize = 4,
    label = "")

    savefig(plt, joinpath(name, "response_contour_zoned.png"))
    println("✅ Saved: response_contour_zoned.png")
end
end # module plot_response_contour


# ---------------------- Plot Energy Coefficients ----------------------
module plot_energy_coefficients

using JLD2
using Plots
using DataFrames

export plot_coefficients

function plot_coefficients(name::String)
    file = joinpath(name, "mem_data.jld2")
    name_tag = splitpath(name)[end] 
    data = load(file)
    println("Loading and plotting energy coefficients...")

    prbPow = Matrix(data["prbPow"])
    ω = data["ω"]

    Pin = prbPow[:, 1]
    Prf = prbPow[:, 2]
    Ptr = prbPow[:, 3]
    Pd_total  = prbPow[:, 4]

    Kr = Prf ./ Pin
    Kt = Ptr ./ Pin
    Ka = Pd_total  ./ Pin
    Err = 1 .- (Kr .+ Kt .+ Ka)

    # 反射系数 Kr
    plot(ω, Kr, label="Kr_$name_tag", lw=2)
    xlabel!("ω (rad/s)")
    ylabel!("Kr")
    title!("Reflection Coefficient vs Frequency")
    savefig(joinpath(name, "Kr_vs_omega.png"))
    println("✅ Saved: Kr_vs_omega.png")

    # 透射系数 Kt
    plot(ω, Kt, label="Kt (Transmitted)", lw=2)
    xlabel!("ω (rad/s)")
    ylabel!("Kt")
    title!("Transmission Coefficient vs Frequency")
    savefig(joinpath(name, "Kt_vs_omega.png"))
    println("✅ Saved: Kt_vs_omega.png")

    # 吸收系数 Ka
    plot(ω, Ka, label="Ka (Absorbed)", lw=2)
    xlabel!("ω (rad/s)")
    ylabel!("Ka")
    title!("Absorption Coefficient vs Frequency")
    savefig(joinpath(name, "Ka_vs_omega.png"))
    println("✅ Saved: Ka_vs_omega.png")

    # 4 in 1 plot
    plot(ω, Kr, label="Kr (Reflected)", lw=2)
    plot!(ω, Kt, label="Kt (Transmitted)", lw=2)
    plot!(ω, Ka, label="Ka (Absorbed)", lw=2)
    plot!(ω, Err, label="Energy Error", lw=1, linestyle=:dash)

    xlabel!("ω (rad/s)")
    ylabel!("Energy Coefficient")
    title!("Energy Coefficients vs Frequency")
    savefig( name * "/energy_coeffs.png")
    println("Plot saved to: $(name)/energy_coeffs.png")


    # # 探针点响应图（随频率）
    # # --------------------------
    # ω = data["ω"]
    # k = data["k"]
    # H0 = 10.0  # 你模拟中用的 H0，用于计算 kh
    # prbxy = data["prbxy"]
    # prbDa = data["prbDa"]
    # prbDa_x = data["prbDa_x"]
    # prbPow = data["prbPow"]

    #   for lprb in 1:length(prbxy)
    #     xloc = prbxy[lprb][1]

    #     plt1 = plot(k .* H0, abs.(prbDa[:, lprb]), lw=2,
    #         xlabel = "kh", ylabel = "|A| (m)", title = "Amplitude")

    #     plt2 = plot(k .* H0, abs.(prbDa_x[:, lprb]), lw=2,
    #         xlabel = "kh", ylabel = "|dA/dx|", title = "Slope Magnitude")

    #     plt3 = plot(k .* H0, angle.(prbDa[:, lprb]), lw=2,
    #         xlabel = "kh", ylabel = "∠A (rad)", title = "Phase")

    #     plt4 = plot(k .* H0, angle.(prbDa_x[:, lprb]), lw=2,
    #         xlabel = "kh", ylabel = "∠dA/dx (rad)", title = "Slope Phase")

    #     plt_all = plot(plt1, plt2, plt3, plt4, layout=(4,1), size=(600, 800),
    #         plot_title = "Probe at x = $xloc")

    #     savefig(plt_all, name * "/_dxPrb_$lprb.png")
    #   end
end
end # module plot energy coefficients