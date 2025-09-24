module Membrane_modes_lrmm

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using LinearAlgebra
using TickTock
using DataFrames
using Printf


function run_case(rMfac = 1, rωfac = 2.4)

  function run_freq(ω)

    k = dispersionRelAng(H0, ω; msg=false)
    @show ω, k

    # Weak form: ω dependent
    k22(ϕ,w) = ∫( ∇(w)⋅∇(ϕ) )dΩ +
        ∫( -w * im * k * ϕ )dΓin + ∫( -w * im * k * ϕ )dΓot
        
    c23(κ,w) = ∫( im*ω*w*κ )dΓfs 
    c32(ϕ,u) = ∫( -im*ω*u*ϕ )dΓfs 
      
    # Weak form ：ω dependent resonator
    k11r(η, v) = - (im*ω*rC - rKᵨ) * δ_p(v * η)
    k44(q, ξ) = (rKᵨ - im*ω*rC) * δ_p(q⋅ξ)
    c14(q, v) = (im*ω*rC - rKᵨ) * δ_p(v * (q⋅î1))
    c41(η, ξ) = -(-im*ω*rC + rKᵨ) * δ_p((ξ⋅î1) * η)

    # Global matrices: ω dependent
    K11r = get_matrix(AffineFEOperator(k11r, l1, U_Γη, V_Γη))
    K22 = get_matrix(AffineFEOperator( k22, l2, U_Ω, V_Ω ))
    C23 = get_matrix(AffineFEOperator( c23, l2, U_Γκ, V_Ω ))
    C32 = get_matrix(AffineFEOperator( c32, l3, U_Ω, V_Γκ ))
    K44 = get_matrix(AffineFEOperator(k44, l4, U_Γq, V_Γq))
    C14 = get_matrix(AffineFEOperator(c14, l1, U_Γq, V_Γη))
    C41 = get_matrix(AffineFEOperator(c41, l4, U_Γη, V_Γq))

    K11r = real.(K11r)
    K22 = real.(K22)
    C23 = real.(C23)
    C32 = real.(C32)
    K44 = real.(K44)
    C14 = real.(C14)
    C41 = real.(C41)

    # tick()
    Mϕ = K22 - ( C23 * (Matrix(K33) \ C32) )
    Mhat = C12 * (Mϕ \ C21)

    Mtot = M11 + Mhat

    # # Spring-mass coupling: this doesn't work
    # Rosc = -C14 * ((-ω^2 * M44 + K44) \ C41)
    # Ktot = K11 + K11r + Rosc
    # Meff = Mtot - (1/ω^2) * Rosc
    # tock()

    # Eigen values
    # λ = LinearAlgebra.eigvals(Mtot\Matrix(Ktot))
    # V = LinearAlgebra.eigvecs(Mtot\Matrix(Ktot))

    # couple with: K44 this way does not make sense
    # meff = diag((V[1:end-1, 1:nωₙ])' * Meff * V[1:end-1, 1:nωₙ])
    
    A = [K11 + K11r     C14;
         C41            K44] 
    B = [Mtot          0*C14;
         0*C41          M44]

    λ, V = eigen(Matrix(B) \ Matrix(A))
    # @show sum(imag.(λ))

    #this is not right cuz now we solve EVP system A x = λ B x ,B is the mass
    # meff = diag((V[1:end-1, 1:nωₙ])' * Mtot * V[1:end-1, 1:nωₙ])

    meff = diag(transpose(V[:, 1:nωₙ])* Matrix(B) * V[:, 1:nωₙ])

    # # H transpose
    # Mherm = real(Mtot) #  should work the same if syymmetric
    # meff = diag((V[1:end-1, 1:nωₙ])' * Mherm * V[1:end-1, 1:nωₙ])

    # ωₙ = real.(sqrt.(Complex.(λ))) #damped frequency
    ωₙ = real.(sqrt.(λ)) # undamped frequency
    # α = -imag.( sqrt.(Complex.(λ)) ) # decay rate

    # return(ωₙ[1:nωₙ], V[:,1:nωₙ])
    # return(ωₙ[1:nωₙ], α[1:nωₙ], V[:,1:nωₙ])
    return(ωₙ[1:nωₙ], V[:,1:nωₙ], meff)    
  end
  


  caseName = "rω" * @sprintf("%0.2f", rωfac) *"_rmass" * @sprintf("%0.2f", rMfac)
  name::String = "data/sims_memmodes/lrmm_modes_wet_realpart/lrmm_modes_"*caseName
  order::Int = 2
  vtk_output::Bool = true
  filename = name*"/mem"

  # 跳过已存在的文件夹的计算
  # if( isdir(name) )
  #   return
  # end

  mkpath(name)
  @info "▶ Running case rMᵨ = $rMfac, rω = $rωfac"

  ρw = 1025 #kg/m3 water
  H0 = 10 #m #still-water depth

  # Membrane parameters
  Lm = 2*H0 #m
  @show g #defined in .Constants
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1*g*H0*H0 #T/ρw
  
  # Resonator parameters
  # rM = 1.0e3   # kg
  # rK = 5.9e3  # N/m
  rMᵨ = rMfac
  rω = rωfac
  rKᵨ = rω*rω*rMᵨ  # N/m
  @show rKᵨ 
  ζ = 0.0     # damping ratio
  rC = 2 * ζ * sqrt(rKᵨ * rMᵨ)  # N·s/m
  mass_ratio = rMᵨ / (ρw*Lm*mᵨ)

  # Excitation wave parameters
  ω = 1.0

  # Domain 
  nx = 120
  ny = 10
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


  # Weak form:  ω independent constant matrices
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  m11(η,v) = ∫( mᵨ*v*η )dΓm
  k11(η,v) = ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm #+  
              #∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

  c12(ϕ,v) = ∫( v*ϕ )dΓm
  c21(η,w) = ∫( w*η )dΓm  

  k33(κ,u) = ∫( u*g*κ )dΓfs

  # Spring-mass-damper oscillator coupling terms
  m44(q, ξ) = rMᵨ * δ_p(q⋅ξ)

  l1(v) = ∫( 0*v )dΓm
  l2(w) = ∫( 0*w )dΩ
  l3(u) = ∫( 0*u )dΓfs # + ∫( 0*u )dΓd1 + ∫( 0*u )dΓd2
  zero_vec = VectorValue(0.0+0im)
  l4(ξ) = ∫( zero_vec ⋅ ξ )dΩ
  println("[MSG] Done Weak form")

  # Global matrices
  M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
  K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
  C12 = get_matrix(AffineFEOperator( c12, l1, U_Ω, V_Γη ))
  C21 = get_matrix(AffineFEOperator( c21, l2, U_Γη, V_Ω ))
  K33 = get_matrix(AffineFEOperator( k33, l3, U_Γκ, V_Γκ ))
  M44 = get_matrix(AffineFEOperator(m44, l4, U_Γq, V_Γq))
  println("[MSG] Done Global matrices")
  println(K11 == transpose(K11))

  # Iterative algorithm
  maxIter = 20
  nωₙ = 7
  da_ωₙ = zeros(Float64, 1, nωₙ)
  # @show ωₙ=zeros(Float64, 1, nωₙ) .+ ω
  @show ωₙ = [0.5, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
  println("[MSG] Starting iterative solution for wet natural frequencies")
  da_V = []
  da_meff=[]


  for i in 1:nωₙ
    # global da_ωₙ, da_V  
    # global ωₙ, ω
    local V, lIter, meff
    lIter = 0    
    Δω = 1
    ω = ωₙ[i]
    # ω = 1.0
    tick()
    while ((Δω > 1e-3) && (lIter < maxIter))
      # global ω, ωₙ
      
      ωₙ, V, meff = run_freq(ω)
      ωₒ = ω
      ωᵣ = ωₙ[i]     
      # ωᵣ = real(sqrt(λ[i]))

      # if(i==1)
      #   #ω = 0.2 * ωₙ[i] + 0.8*ω
      #   ω = 0.0
      #   Δω = 0.0
      #   V = V*0.0
      # elseif(i==4)
      #   ω = 0.4 * ωᵣ + 0.6*ωₒ
      #   Δω = abs(ω - ωₒ)/ωₒ
      # else
      #   ω = 0.5 * ωᵣ + 0.5*ωₒ
      #   Δω = abs(ω - ωₒ)/ωₒ
      # end      
    
      ω = 0.8 * ωᵣ + 0.2*ωₒ
      Δω = abs(ωᵣ - ωₒ)/ωₒ
      
      lIter += 1
      @show ωₙ
      @show i, ω, Δω, lIter
    end
    tock()
    println("---- Found ωₙ[$i] = $ω in $lIter iterations ----")

    da_ωₙ[i] = ω
    push!(da_V, V[:,i])
    push!(da_meff, meff[i])
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
  "meff" => da_meff,
  )

  wsave(filename*"_modesdata.jld2", data)
end


# rMfac = [100, 500, 3000]
# rωfac = [1.01, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50 ]

# rMfac = 1000
# rωfac = [1.01, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00 ]

# rMfac = [0.1, 0.3, 0.5]
# rωfac = [1.5, 3.5, 4.5]

rMfac = [0.1, 0.3, 0.5, 1.0]
rωfac = [2.4]

# rMfac = 0.1:0.1:1.0
# rωfac = 2.40

for irMfac in rMfac
  for irωfac in rωfac
    run_case(irMfac, irωfac)
  end
end

end