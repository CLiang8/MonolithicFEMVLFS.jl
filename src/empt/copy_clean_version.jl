module Empty2D_design

using Gridap
using Plots
using DrWatson
using WaveSpec
using WaveSpec.Constants
using WaveSpec.Jonswap
using DataFrames:DataFrame
using DataFrames:Matrix
using TickTock
using Parameters
using Printf
using MonolithicFEMVLFS.Resonator



function main(params)

  function run_freq(ω, η₀, α)

    tick()
    # Wave parameters
    k = dispersionRelAng(H0, ω; msg=false)
    λ = 2π/k   
    T = 2π/ω
    ηᵢₙ(x) = η₀*exp(im*k*x[1] + im*α)
    ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(H0 + x[2])) / 
      sinh(k*H0))*exp(im*k*x[1] + im*α)
    vxᵢₙ(x) = (η₀*ω)*(cosh(k*(H0 + x[2])) / 
      sinh(k*H0))*exp(im*k*x[1] + im*α)
    vzfsᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1] + im*α) #???
    @show ω, T
    @show λ
    @show η₀
    @show H0, H0/λ

    # Numeric constants
    αₕ = -im*ω/g * (1-βₕ)/βₕ
    @show αₕ
    println()  

    # Damping
    # Ldw = min( 5.0*λ, Ld )
    μ₀ = 2.5#maximum([2.5, 5.24/(ω^0.922)])#2.5
    μ₁ᵢₙ(x) =  μ₀*(1.0 - sin(π/2 * (x[1]-x₀)/Ld ))
    # μ₁ᵢₙ(x) =  μ₀*(1.0 - sin(π/2 * min( (x[1]-x₀)/Ldw, 1.0 ) ))
    μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
    # μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
    ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
    ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???

    # Weak form
    ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)

    function a(trialVars, testVars)
      (ϕ,κ) = trialVars
      (w,u) = testVars

      val = ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
        ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
        ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
          - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
        ∫( -w * im * k * ϕ )dΓot 
        # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
        #   - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
        # ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
          # - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  #+ 
        #∫(- Tᵨ*(1-im*ω*τ)*v*∇(η)⋅nΛmb )dΛmb

      # for (qi, ξi, δi, rSi) in zip(q, ξ, δ_p_Arr, rS)
      #   val += 
      #     (+im*ω*rSi.C -rSi.K)/ρw*δi( v*( (qi⋅î1) - η ) ) +    
      #     ∫( (ξi⋅qi)* 0.0 )dΩ + 
      #     # ∫( -rM/cnstFEArea*ω^2*(q⋅ξ) + rK/cnstFEArea*(ξ⋅q) )dΩ +
      #     -rSi.M*ω^2*δi(qi⋅ξi) +
      #     (-im*ω*rSi.C + rSi.K)*δi(qi⋅ξi - (ξi⋅î1)*η)    
      # end

      return val
    end

    l((w,u)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1


    # Solution
    op = AffineFEOperator(a,l,X,Y)
    (ϕₕ, κₕ) = solve(op)

    # Function for inlet phase
    κin = interpolate_everywhere(ηᵢₙ, 
      FESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64}))
    
    # Wave energy flux
    # ηrf = abs(κr(Point(prbPowx[1],0.0)))
    # ηtr = abs(κₕ(Point(prbPowx[2],0.0)))
    kh = k*H0
    wave_n = 0.5*(1 + 2*kh/sinh(2*kh))
    Pin = (0.5*ρw*g*η₀*η₀)*(ω/k)*wave_n
    # Prf = (0.5*ρw*g*ηrf*ηrf)*(ω/k)*wave_n
    # Ptr = (0.5*ρw*g*ηtr*ηtr)*(ω/k)*wave_n
    # PErr = Pin - Prf - Ptr - Pd - sum( Pd_r )
    println("Power In \t ",Pin,"  W/m")
    # println("Power Ref \t ",Prf," W/m")
    # println("Power Trans \t ",Ptr," W/m")
    # println("Power Abs \t ",Pd," W/m")
    # println("Power Abs Resonator \t ",Pd_r," W")
    # println("Error \t ",PErr," W/m")    

    # Interpolation on prboes
    # prb_κ = zeros(ComplexF64, 1, length(prbxy))
    # prb_κ_x = zeros(ComplexF64, 1, length(prbxy))  
    

    # prb_κ = κₕ(prbxy)
    # prb_κ_x= (∇(κₕ)⋅VectorValue(1.0,0.0))(prbxy)
  
    # push!(prbDa, prb_κ)  
    # push!(prbDa_x, prb_κ_x)  

    # push!(prbDaΓη, ηₕ(prxΓη))
    # push!(prbDaΓκ, κₕ(prxΓκ))
    # push!(prbDaΓκr, κr(prxΓκ))

    # push!(prbPow, [Pin, Prf, Ptr, Pd, PErr, 0.0, Pd_r...])

    # 新增 cylinder pressure/force 计算
    P_cyl = ρw * im * ω * ϕₕ
    Fcyl = sum(∫( P_cyl * nΓot )dΓot)
    # @show P_cyl
    @show Fcyl
    push!(prbForce, [ω, Fcyl[1], Fcyl[2]]) # x, y 分量

    # p_cyl = -im*ω*ρw*ϕₕ( Point(x_cyl, 0.0) )
    # 如果你想只记录 Fx 的模长，也可以用：
    # push!(prbForce, [ω, abs(Fcyl[1])])
    
    # VTK Output
    # ---------------------Start---------------------
    if vtk_output == true
      
      freqName = filename*"_ω_" * @sprintf("%.3f", ω)

      mkpath(freqName)

      writevtk(Ω, freqName * "/mem_O_sol.vtu",
        cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),
        "phi_abs" => abs(ϕₕ), "phi_ang" => angle∘(ϕₕ)])
      
      writevtk(Γκ, freqName * "/mem_Gk_sol.vtu",
        cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ),
        "eta_abs" => abs(κₕ), "eta_ang" => angle∘(κₕ),
        # "etaR_re" => real(κr),"etaR_im" => imag(κr),
        # "etaR_abs" => abs(κr), "etaR_ang" => angle∘(κr),
        "ηin_abs" => abs(κin), "ηin_ang" => angle∘(κin)])

      # writevtk(Γη, freqName * "/mem_Ge_sol.vtu",
      #   cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ),
      #   "eta_abs" => abs(ηₕ), "eta_ang" => angle∘(ηₕ)])
    end
    # ----------------------End----------------------
    
    tock()
    return 0
  end
  


  @unpack name, order, vtk_output = params
  @show name
  @show order
  @show vtk_output
  filename = name*"/mem"
  mkpath(name)

  @unpack H0, ω, T, η₀, α = params 
  k = dispersionRelAng.(H0, ω; msg=false)
  
  ρw = 1025 #kg/m3 water
  @show H0  #m #still-water depth
  @show ω

  # Peak Wave
  ωₚ, indp = findmax(η₀)
  @show ωₚ = ω[indp]
  kₚ = dispersionRelAng(H0, ωₚ; msg=false)
  println("Peak Wave T, L ", 2*pi/ωₚ, " ", 2*pi/kₚ)


  # Membrane parameters
  @unpack Lm, mᵨ, Tᵨ, τ = params
  @show Lm  #m
  @show g #defined in .Constants
  @show mᵨ #mass per unit area of membrane / ρw
  @show Tᵨ #T/ρw
  @show τ #damping coeff


  # Domain 
  @unpack nx, ny, mesh_ry, Ld, LΩ, x₀ = params
  @unpack domain, partition, xdᵢₙ, xm₀, xm₁ = params
  @show Lm
  @show LΩ, Ld
  @show domain
  @show partition
  @show mesh_ry
  @show (xm₀, xm₁)
  @show isinteger(Lm/LΩ*nx)
  @show LΩ/nx
  @show H0/ny
  println()


  # Numeric constants
  h = LΩ / nx
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  # αₕ = -im*ω/g * (1-βₕ)/βₕ
  @show h
  @show βₕ
  # @show αₕ
  println()



  # Mesh
  # ---------------------Start---------------------
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
  map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0; dbgmsg=false) )
  model = CartesianDiscreteModel(domain,partition,map=map)
  # ----------------------End----------------------


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
  function is_damping1(xs) # Check if an element is inside the damping zone 1
    n = length(xs)
    x = (1/n)*sum(xs)
    (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
  end
  # function is_damping2(xs) # Check if an element is inside the damping zone 2
  #   n = length(xs)
  #   x = (1/n)*sum(xs)
  #   (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
  # end
  
  # # 添加 monopile 边界定义代码
  # @show x_cyl = 130 + 33H0/(2*nx) # xm₁ + 3H0  # 可以设为 xm₁ + 0.5*H0 等
  # function is_monopile(xs)
  #     n = length(xs)
  #     x = (1/n) * sum(xs)
  #     return abs(x[1] - x_cyl) < 0.05
  # end
  # xΓ = get_cell_coordinates(Γ)
  # Γcyl_to_Γ_mask = lazy_map(is_monopile, xΓ)
  # Γcyl = Triangulation(Γ, findall(Γcyl_to_Γ_mask))
  # if vtk_output == true
  #   writevtk(Γcyl, filename * "_Gcyl")
  # end


  # Masking and Beam Triangulation
  xΓ = get_cell_coordinates(Γ)
  Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
  Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
  #Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
  Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
  #Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask .| 
    Γd1_to_Γ_mask ))# .| Γd2_to_Γ_mask))
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
    writevtk(Γd1,filename*"_Gd1")
    # writevtk(Γd2,filename*"_Gd2")
    writevtk(Γfs,filename*"_Gfs")
    writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree)
  dΓd1 = Measure(Γd1,degree)
  # dΓd2 = Measure(Γd2,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  dΛmb = Measure(Λmb,degree)
  # dΓcyl = Measure(Γcyl,degree)


  # Normals
  # @show nΛmb = get_normal_vector(Λmb)
  nΓot = get_normal_vector(Γot)
  @show nΓot


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
  # U_Γη = TrialFESpace(V_Γη, gη)
  U_Γη = TrialFESpace(V_Γη)


  X = MultiFieldFESpace([U_Ω, U_Γκ])
  Y = MultiFieldFESpace([V_Ω, V_Γκ])

  # Probes
  # @unpack prbx = params
  # prbxy = Point.(prbx, 0.0)

  # 存储表（首行占位，最后丢弃）
  # lDa = zeros(ComplexF64, 1, length(prbxy))
  # prbDa   = DataFrame(lDa, :auto)   # η
  # prbDa_x = DataFrame(lDa, :auto)   # ∂η/∂x（用势替代）
  prbForce = DataFrame(zeros(ComplexF64, 1, 3), :auto)


  # Remove old vtk files  
  for entry in readdir(name)
    if startswith(entry, "mem_ω")
      rm(joinpath(name, entry); force=true, recursive=true)
    end
  end    

  # Run weak-form for each freq
  # ---------------------Start---------------------
  run_freq.(ω, η₀, α)
  # ----------------------End----------------------

  # prbDa = prbDa[2:end, :]
  # prbDa_x = prbDa_x[2:end, :]
  # prbDaΓη = prbDaΓη[2:end,:]
  # prbDaΓκ = prbDaΓκ[2:end,:]
  # prbDaΓκr = prbDaΓκr[2:end,:]
  # prbPow = prbPow[2:end,:]
  prbForce = prbForce[2:end,:]

  k = dispersionRelAng.(H0, ω; msg=false)

  data = Dict("ω" => ω,
              "η₀" => η₀,
              "k" => k,
              # "prbxy" => prbxy,
              # "prbDa" => prbDa,            
              # "prbDa_x" => prbDa_x,
              # "prxΓκ" => prxΓκ,
              # "prxΓη" => prxΓη,
              # "prbDaΓκ" => prbDaΓκ,
              # "prbDaΓη" => prbDaΓη,
              # "prbDaΓκr" => prbDaΓκr,
              # "prbPow" => prbPow,
              "prbForce" => prbForce,
              # "rS" => rS
              )

  wsave(filename*"_data.jld2", data)

end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct Memb_params_warmup
  name::String = "data/sims_202509/empty"
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [2*π/2.53079486745378, 2*π/2.0]
  # η₀ = [0.25, 0.25]
  ω = 1:0.05:5
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)
  # k = dispersionRelAng.(H0, ω; msg=false)

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff


  # Domain 
  nx = 1650
  ny = 20
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  Ld = 15*H0 #damping zone length
  LΩ = 18*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 8*H0
  xm₁ = xm₀ + Lm


  # Resonator parameters
  rS = Resonator.Array1D(
    1, 
    1e3, 
    5.9e3, 
    0.0,
    [Point(xm₀ + Lm/2.0,0.0)]
  )


  # Probes
  prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
          52.7, 53.7, 55, 60.0, 80.0, 
          85.0, 90.0, 95.0, 100.0, 120.0, 
          125.0, 140.0, 160.0, 180.0 ]
  prbPowx=[ 55.0, 125.0 ]

end




# ----------------------------------------------------------execution-----------------------------------------------------------

# params = Memb_params_warmup()
# params = Memb_params()
# main(params)


end