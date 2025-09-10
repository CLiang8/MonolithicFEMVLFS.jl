module Membrane_modes

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using LinearAlgebra
using TickTock
using DataFrames


function run_freq(ω)

  k = dispersionRelAng(H0, ω; msg=false)
  @show ω, k

  # Weak form: ω dependent
  c14(q, v) = (im*ω*rC - rkᵨ) * δ_p(v * (q⋅î1))
  c41(η, ξ) = -(-im*ω*rC + rkᵨ) * δ_p((ξ⋅î1) * η)

  # Global matrices: ω dependent
  C14 = get_matrix(AffineFEOperator(c14, l1, U_Γq, V_Γη))
  C41 = get_matrix(AffineFEOperator(c41, l4, U_Γη, V_Γq))

  # # Solution
  # tick()
  # Rosc = -C14 * ((-ω^2 * M44 + K44) \ C41)
  # Ktot = K11 + K11r + Rosc
  # tock()

  # # Eigen values
  # λ = LinearAlgebra.eigvals(M11\Matrix(Ktot))
  # V = LinearAlgebra.eigvecs(M11\Matrix(Ktot))  
  
  # solution as a GEP
  A = [K11 + K11r     C14;
       C41            K44]
  B = [M11            0*C14;
       0*C41          M44]

  tick()
  λ, V = eigen(Matrix(B) \ Matrix(A))  # or eigen(A, B)
  tock()

  #@show real.(λ[1:nωₙ])
  @show sum(imag.(λ))
  ωₙ = real.(sqrt.(λ))
  return(ωₙ[1:nωₙ], V[:,1:nωₙ])
    
end


name::String = "data/sims_memmodes/mem_modes_dry_lrmm"
order::Int = 1
vtk_output::Bool = true
filename = name*"/mem"
mkpath(name)

ρw = 1025 #kg/m3 water
H0 = 10 #m #still-water depth

# Membrane parameters
Lm = 2*H0 #m
@show g #defined in .Constants
mᵨ = 0.9 #mass per unit area of membrane / ρw
Tᵨ = 0.1*g*H0*H0 #T/ρw

# Resonator parameters
rM = 1000   # kg
# rK = 5.9e3  # N/m
rMᵨ = rM/ρw
rkᵨ = 2.4*2.4*rMᵨ  # N/m
@show sqrt(rkᵨ/rMᵨ)  # rad/s
ζ = 0.0     # damping ratio
rC = 2 * ζ * sqrt(rkᵨ * rMᵨ)  # N·s/m

# Excitation wave parameters
ω = 1.0


# Domain 
nx = 300
ny = 20
mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
LΩ = 6*H0 
x₀ = 0.0
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xm₀ = x₀ + 2*H0
xm₁ = xm₀ + Lm
@show Lm
@show LΩ
@show domain
@show partition
@show (xm₀, xm₁)
@show isinteger(Lm/LΩ*nx)
@show LΩ/nx
@show H0/ny
#@show Ld*k/2/π
#@show cosh.(k*H0*0.5)./cosh.(k*H0)
println()


# Mesh
function f_y(y, r, n, H0; dbgmsg = false)
  # Mesh along depth as a GP
  # Depth is 0 to -H0    
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
    if(dbgmsg)
      ln = 0:n
      ly = -a0 / (r-1) * (r.^ln .- 1)         
      @show hcat( ly, [ 0; ly[1:end-1] - ly[2:end] ] )
    end
    
    if y ≈ 0
      return 0.0
    end
    j = abs(y) / H0 * n  
    return -a0 / (r-1) * (r^j - 1)
  end
end
map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0) )
model = CartesianDiscreteModel(domain,partition,map=map)


# Labelling
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)


# Triangulations
Ω = Interior(model) #same as Triangulation()
Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
Γin = Boundary(model,tags="inlet")
Γot = Boundary(model,tags="outlet")


# Auxiliar functions
function is_mem(xs) # Check if an element is inside the beam1
  n = length(xs)
  x = (1/n)*sum(xs)
  (xm₀ <= x[1] <= xm₁ ) * ( x[2] ≈ 0.0)
end


# Masking and Beam Triangulation
xΓ = get_cell_coordinates(Γ)
Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask))
Γη = Triangulation(Γ, findall(Γm_to_Γ_mask))
Γκ = Triangulation(Γ, findall(!,Γm_to_Γ_mask))


# Construct the tag for membrane boundary
Λmb = Boundary(Γm)
xΛmb = get_cell_coordinates(Λmb)
xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
new_entity = num_entities(labels_Ω) + 1
labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
add_tag!(labels_Ω, "mem_bnd", [new_entity])


writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γ,filename*"_G")
  writevtk(Γm,filename*"_Gm")  
  writevtk(Γfs,filename*"_Gfs")
  writevtk(Λmb,filename*"_Lmb")  
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓm = Measure(Γm,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)
dΓot = Measure(Γot,degree)
dΛmb = Measure(Λmb,degree)

# Dirac delta
xr::Float64 = 30.0
δ_p = DiracDelta(Γ, [Point(xr, 0.0)])


# Normals
@show nΛmb = get_normal_vector(Λmb)


# Dirichlet Fnc
gη(x) = ComplexF64(0.0)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
# V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
#   vector_type=Vector{ComplexF64},
#   dirichlet_tags=["mem_bnd"]) #diri
V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γκ = TrialFESpace(V_Γκ)
# U_Γη = TrialFESpace(V_Γη, gη) #diri
U_Γη = TrialFESpace(V_Γη)

V_Γq = ConstantFESpace(Ω, vector_type=Vector{ComplexF64}, 
  field_type=VectorValue{1,ComplexF64})
U_Γq = TrialFESpace(V_Γq)
î1 = VectorValue(1.0)

# Weak form
∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
m11(η,v) = ∫( mᵨ*v*η )dΓm
k11(η,v) = ∫( Tᵨ*∇(v)⋅∇(η) )dΓm #+  
            #∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

# Spring-mass-damper oscillator coupling terms
k11r(η, v) = - (im*ω*rC - rkᵨ) * δ_p(v * η)
m44(q, ξ) = rMᵨ * δ_p(q⋅ξ)
k44(q, ξ) = (rkᵨ - im*ω*rC) * δ_p(q⋅ξ)

l1(v) = ∫( 0*v )dΓm
zero_vec = VectorValue(0.0+0im)
l4(ξ) = ∫( zero_vec ⋅ ξ )dΩ
println("[MSG] Done Weak form")


# Global matrices
M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
K11r = get_matrix(AffineFEOperator(k11r, l1, U_Γη, V_Γη))

M44 = get_matrix(AffineFEOperator(m44, l4, U_Γq, V_Γq))
K44 = get_matrix(AffineFEOperator(k44, l4, U_Γq, V_Γq))
println("[MSG] Done Global matrices")

#xp = range(xm₀, xm₁, size(V,2)+2)

nωₙ = 10
da_ωₙ = zeros(Float64, 1, nωₙ)
@show ωₙ=zeros(Float64, 1, nωₙ) .+ ω
da_V = []

# # For index=1 not looping coz ωₙ[1] = 0.0
# i = 1
# ωₙ, V = run_freq(ω)
# da_ωₙ[i] = ωₙ[1]
# push!(da_V, V[:,i])

for i in 1:nωₙ
  global da_ωₙ, da_V  
  global ωₙ, ω
  local V
  Δω = 1
  ω = ωₙ[i]
  while Δω > 1e-4
    global ω, ωₙ
    ωₙ, V = run_freq(ω)
    Δω = abs(ωₙ[i] - ω)
    if(i==1)
      ω = 0.2 * ωₙ[i] + 0.8*ω
      # ω = 0.0
      # Δω = 0.0
      # V = V*0.0
    else
      ω = 0.8 * ωₙ[i] + 0.2*ω
    end
    @show ωₙ
    @show i, ω, Δω
    @show "--------------next iter--------------"
  end
  da_ωₙ[i] = ω
  push!(da_V, V[:,i])
end

println(da_ωₙ)

# xp = range(xm₀, xm₁, length(da_V[1]))

# data = Dict(
#   "xp" => xp,
#   "ωₙ" => da_ωₙ,
#   "V" => da_V  
# )

# solution 2 data saving
nη = size(M11, 1)
nq = size(M44, 1)

xp = range(xm₀, xm₁, nη)
η_all = [v[1:nη] for v in da_V]
q_all = [v[nη+1:end] for v in da_V]

data = Dict(
  "xp" => xp,
  "ωₙ" => da_ωₙ,
  "V" => η_all,
  "q_modes" => q_all,
)

# wsave("$(filename)_modesdata_m=$(rM).jld2", data)
wsave(filename * "_modesdata_lrmm.jld2", data)

end


#------------------------------plotting--------------------------------------------
using Plots, JLD2

rM = 1000.0  # 根据你的绘图脚本中的 rM 值
base_dir = "data/sims_memmodes/mem_modes_dry_lrmm"
filename = "mem"
path = "$(base_dir)/$(filename)_modesdata_m=$(rM).jld2"

# === Load from Dict ===
data = load(path)
xp = data["xp"]
ωₙ = data["ωₙ"]
η_all = data["V"]
q_all = data["q_modes"]

N = 6  # 模态总数（通常为你保存的 nωₙ）

# === 1. 多模态 η(x) subplot 绘图 ===
plt1 = plot(layout = (N, 1), size = (800, 250 * N), legend = false)

for n in 1:N
    ηn = real.(η_all[n])./ maximum(abs.(η_all[n]))  # 可换成 abs.() 查看振幅模式
    ω_str = "ω = $(round(ωₙ[n], digits=3)) rad/s"
    plot!(
        plt1[n], xp, ηn,
        lw = 2,
        xlabel = "x (m)",
        ylabel = "η(x)",
        ylims = (-1.1, 1.1),
        title = "Mode $n — $ω_str",
        grid = true,
    )
end

savefig(plt1, filename * "_modes_subplot.png")
display(plt1)

# === 2. q 的 modal participation 柱状图 ===
plt2 = bar(abs.([q[1] for q in q_all]),
    xlabel = "Mode number",
    ylabel = "|q|",
    title = "Modal participation of resonator q",
    legend = false,
    xticks = 1:N,
    bar_width = 0.6,
    lw = 0.5,
    grid = true,
)

savefig(plt2, filename * "_q_modal_participation.png")
display(plt2)
