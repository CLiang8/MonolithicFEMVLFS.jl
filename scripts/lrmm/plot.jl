module plot_response_contour

using Gridap
using JLD2
using Plots
using Measures
using DataFrames
using WaveSpec.Constants
using LaTeXStrings
using Colors
using Measures

# load data
name = "data/sims_mem_freq_lrmm/Coefficients"
filename = name * "/mem_data.jld2"
data = load(filename)

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
    ylims = (1, 5),
    aspect_ratio = 2/3,
    title = "Normalized Amplitude",
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
ω_ticks = range(minimum(ω), stop=maximum(ω), length=5)
hline!(ω_ticks, lw=1, linestyle=:dot, color=:white, label="")


savefig(plt, joinpath(name, "response_contour_zoned.png"))
println("✅ Saved: response_contour_zoned.png")
end