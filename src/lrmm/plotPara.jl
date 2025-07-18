module Plotting

using DrWatson
using Gridap
using DelimitedFiles
using DataFrames:DataFrame
using DataFrames:Matrix
using Revise
using Plots
using Plots.PlotMeasures
using LinearAlgebra
using Test
using Printf
using LaTeXStrings

# export plotRes

# function plotRes()
# global κₕ


# function read_case(plt1; mfac = 0.9, tfac = 0.1, modeN=1, lsvar=1)

#   caseName = "ten" * @sprintf("%0.2f", tfac) *"_mass" * @sprintf("%0.2f", mfac)
#   name = "data/sims_202403/mem_modes_res_free_eigen/mem_modes_"*caseName
#   filename = name*"/mem_modesdata.jld2"

#   da = wload(filename)
  
#   #xp = da["xp"]
#   ω = da["ωₙ"]
#   V = da["V"]

#   # Plots Wet
#   i = modeN
#   ωi = ω[i]
#   Vi = V[i]
#   Vi = Vi * sign(Vi[4])
#   Vi = real.(Vi)/norm(real.(Vi))
#   ωi_rnd = round(ωi; digits=2)

#   # SingleFieldFunction Wet
#   nx = length(Vi)-1
#   fva(x) = 0.0
#   reffe = ReferenceFE(lagrangian,Float64,2)
#   memModel = CartesianDiscreteModel((0.0, Lm), nx/2 )
#   fv = interpolate_everywhere( fva, FESpace( memModel, reffe ) )

#   # [ println( (i,ivi) ) for (i,ivi) in zip( range(1,nx+1), Vi ) ] 
  
#   @test length(fv.free_values) == length(Vi)
#   for i in 1:nx+1
#     fv.free_values[i] = Vi[i]
#   end
  
#   xp = range(0, Lm, nx+1)
#   xpPoi = [Point(ixp) for ixp in xp]

#   tlbl = tfac/4
#   mlbl = mfac / Lm

#   plot!(plt1, (xp.+80)/Lm, abs.(fv.(xpPoi))/Lm*50 .+ 0.25, 
#     lw=5, label=nothing, lc=:green, ls = :dash)    
  
  
# end

filename = "data/sims_202507/mono_freq_free/mem_data.jld2"
da = wload(filename)
ηₕ = da["ηₕ"]
κₕ = da["κₕ"]
κr = da["κr"]
κin = da["κin"]

dx = 0.2
amp = 0.1

setlw = 5

# # Contour
# px = 50:0.5:130
# py = -10:0.5:0
# pp = @.Point(px',py) #wow
# contourf!(plt1, px, py, abs.(ϕₕ.(pp)),
#   color = :turbo, levels=20)

plt1 = plot()

# read_case(plt1, modeN=3)

# Membrane
xka = 80:dx:100
prbxka = Point.(xka, 0.0)
plot!(plt1, xka / Lm, abs.(ηₕ.(prbxka)) / amp, 
    lc=:red, lw = setlw, label = false)
  

# Reflected
xka = 50:dx:80
prbxka = Point.(xka, 0.0)
plot!(plt1, xka / Lm, abs.(κr.(prbxka)) / amp, 
  lc=:blue, lw = setlw, label = false)


# Incident
xka = 50:dx:80
prbxka = Point.(xka, 0.0)
plot!(plt1, xka / Lm, abs.(κin.(prbxka)) / amp, 
    lc=:blue, lw = setlw, ls = :dot, label = false)

xka = 100:dx:130
prbxka = Point.(xka, 0.0)
plot!(plt1, xka / Lm, abs.(κin.(prbxka)) / amp, 
    lc=:blue, lw = setlw, ls = :dot, label = false)
  

# Free-surface
xka = 50:dx:80
prbxka = Point.(xka, 0.0)
plot!(plt1, xka / Lm, abs.(κₕ.(prbxka)) / amp, 
    lc=:black, lw = setlw, label = false)

xka = 100:dx:130
prbxka = Point.(xka, 0.0)
plot!(plt1, xka / Lm, abs.(κₕ.(prbxka)) / amp, 
    lc=:black, lw = setlw, label = false)


lgndlw = 2
xlbl = 0:1:1
ylbl = xlbl * 0
plot!(plt1, xlbl, ylbl,   
  lc = :black, lw = lgndlw, label = L"|\, \kappa \, |")
plot!(plt1, xlbl, ylbl,   
  lc = :red, lw = lgndlw, label = L"| \, \eta \, |")
plot!(plt1, xlbl, ylbl,   
  lc = :blue, lw = lgndlw, label = L"| \, \kappa_r \, |")
plot!(plt1, xlbl, ylbl,   
  lc = :blue, lw = lgndlw, ls = :dot, label = L"| \, \kappa_{in} \, |")
plot!(plt1, xlbl, ylbl,   
  lc = :green, lw = lgndlw, ls = :dash,
  label = L"|" )
  # label = latexstring(L"|"," Mode 2 + Rigid Heave ",L"|") )


plot!(plt1, 
  xlim = (2.5, 6.5),
  ylim = (0,2),
  xticks = 2.5:0.5:6.5,
  yticks = 0.0:0.5:2.0)

plot!(plt1, 
  xlabel = L"x/L_m",
  ylabel = L"| \eta | /\kappa_0",
  title = L"\textit{MEMB1} \, | \, ω = 2.4 \, rad \, s^{-1} \, | \, \tau = 0.0")

vline!(plt1, (2.5, 2.5), lw = 3, color=:black, label=false )
hline!(plt1, (0, 0), lw = 3, color=:black, label=false )

plot!( plt1, 
  titlefontsize=30,
  tickfontsize=25, 
  labelfontsize=30,
  legendfontsize = 25, 
  left_margin=16mm, right_margin=10mm,
  bottom_margin=18mm, top_margin=10mm)    

plot!(plt1, 
  legend = :topright, legendcolumns = 4,
  dpi = 150, size = (2400, 700),
  grid=:true, gridcolor=:black, 
  gridalpha=0.5, gridlinestyle=:dot,
  gridlinewidth = 3)


savefig(plt1,filename*"_ParaPlot.png")
# end

end