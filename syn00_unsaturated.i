# mpiexec -np 4 archaea-opt -i syn00_unsaturated.i # => Num Local DOFs: 10300]
# 9950 sec. [710 steps] with KT, VanLeer
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -12000
    xmax =  12000
    ymin = -15000
    ymax =  -3000
    elem_type = TRI3
    nx = 200
    ny = 100
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 -9.81 0'
[]

[Variables]
  [porepressure]
    family = LAGRANGE
    order = FIRST
  []
   [temperature]            # [K] parsed initial condition in the ICs block below
    family = LAGRANGE
    order = FIRST
    scaling = 1E-08         # this variable scaling brings the residual R_T to the same order of magnitude that the Pressure residual R_P 
  []
[]

[PorousFlowUnsaturated]
  porepressure = porepressure
  temperature = temperature
  coupling_type = ThermoHydro
  fp = tabulated_water95
  #fp = the_simple_fluid
  relative_permeability_exponent = 2 # default
  relative_permeability_type = Corey # 
  residual_saturation = 0.1
  van_genuchten_alpha = 1E-06
  van_genuchten_m = 0.6
  #stabilization = Full # default: full upwinding
  stabilization = KT
  flux_limiter_type = VanLeer
[]

[AuxVariables]
  [mass_flux_top_sink]
    family = LAGRANGE
    order = FIRST
  []
  [mass_flux_top_source]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_adv_top_sink]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_adv_top_source]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_dif_top_sink]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_dif_top_source]
    family = LAGRANGE
    order = FIRST
  []
  [ptop_auxvar]
    family = LAGRANGE
    order = FIRST
  []
  [perm_auxvar]
    family = MONOMIAL # beware: it is possible than making this monomial makes the code faster
    order = CONSTANT
  []
  [water_kg_per_s_left]
    family = LAGRANGE
    order = FIRST
  []
  [water_kg_per_s_right]
    family = LAGRANGE
    order = FIRST
  []
   [heat_J_per_s_left]
    family = LAGRANGE
    order = FIRST
  []
  [heat_J_per_s_right]
    family = LAGRANGE
    order = FIRST
  []
[]

[Functions]
  [temperature00]
    type = ParsedFunction
    value = '273.15 + ttop + (tbot-ttop)*exp(-3*(y-ybot)/rate)'
    vars = 'ybot   rate ttop tbot'
    vals = '-15000 2000 5    800'
  []
  [porepressure00]
    type = ParsedFunction
    value = '101325 + 9.81*(-y)*1000'          # note: top of domain at -3000 m depth
  []
  [permeability_fun]
    type = ParsedFunction
    value = 'vmin + (vmax-vmin) * exp(-((x-xc)/xsd)^2)'
    vars = 'xc  vmin    vmax    xsd'
    vals = '0.0 1.0E-17 1.0E-14 5000'
  []
[]

[ICs]                                          
  [porepressure_IC]
    type = FunctionIC
    variable = porepressure                    # [Pa]
    function = porepressure00                  # (-y+1000) => atmospheric pressure at the top corners <- (PHY.yseq = [-11., 1]*km)
  []
  [temperature_IC]
    type = FunctionIC
    variable = temperature
    function = temperature00
  []
  [ptop_auxvar_IC]
    type = FunctionIC
    variable = ptop_auxvar
    function = porepressure00
  []
  [perm_auxvar_IC]
    type = FunctionIC
    variable = perm_auxvar
    function = permeability_fun
  []
[]

[BCs]
  [Ptop_sink]
    type = PorousFlowPiecewiseLinearSink  # IntegratedBC
    variable = porepressure
    boundary = top
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    use_mobility = true
    use_relperm = true
    flux_function = 100                    # ~minimum where plumes still look natural
    fluid_phase = 0
    save_in = mass_flux_top_sink
  []
  [Ptop_source]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = top
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure
    pt_vals = '-1e9 0'
    multipliers = '-1e9 0'
    use_mobility = true
    use_relperm = false
    flux_function = 100                    #
    fluid_phase = 0
    save_in = mass_flux_top_source
  []
  [Tbot]
    type = FunctionDirichletBC 
    variable = temperature
    boundary = bottom
    function = temperature00
  []
  [Ttop_adv_sink] # advective heat flux = enthalpy * darcy Flux
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = ptop_auxvar
    boundary = top
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    use_mobility = true
    use_relperm = true
    use_enthalpy = true # shoud this be use_internal_energy instead?
    flux_function = 100 
    fluid_phase = 0
    save_in = heat_flux_adv_top_sink
  []
  [Ttop_adv_source] # advective heat flux = enthalpy * darcy Flux
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = ptop_auxvar
    boundary = top
    pt_vals = '-1e9 0'
    multipliers = '-1e9 0'
    use_mobility = true
    use_relperm = false
    use_enthalpy = true # shoud this be use_internal_energy instead?
    flux_function = 100 
    fluid_phase = 0
    save_in = heat_flux_adv_top_source
  []

  [Ttop_diff_sink] # diffusive heat flux = - lambda * grad_T
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = 278.15                            # [K]    5ºC
    boundary = top
    pt_vals     = '0 1000'                 # x coordinates defining g
    multipliers = '0 1000'                 # y coordinates defining g
    use_thermal_conductivity = true
    flux_function = 0.01                  # allow for hydrothermal plumes to go up to ~100 meters in the ocean
    save_in = heat_flux_dif_top_sink 
  []
  [Ttop_diff_source] # diffusive heat flux = - lambda * grad_T
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = 278.15                            # [K]    5ºC
    boundary = top
    pt_vals     = '-1000 0'
    multipliers = '-1000 0'
    use_thermal_conductivity = true
    flux_function = 0.1
    save_in = heat_flux_dif_top_source 
  []
   [right_water_component0]
    type = PorousFlowOutflowBC
    boundary = right
    variable = porepressure
    mass_fraction_component = 0
    include_relperm = true
    save_in = water_kg_per_s_right
  []
  [left_water_component0]
    type = PorousFlowOutflowBC
    boundary = left
    variable = porepressure
    mass_fraction_component = 0
    include_relperm = true
    save_in = water_kg_per_s_left
  []
  [right_heat_component0]
    type = PorousFlowOutflowBC
    boundary = right
    variable = temperature
    flux_type = heat
    mass_fraction_component = 0
    include_relperm = true
    save_in = heat_J_per_s_right
  []
  [left_heat_component0]
    type = PorousFlowOutflowBC
    boundary = left
    variable = temperature
    flux_type = heat
    mass_fraction_component = 0
    include_relperm = true
    save_in = heat_J_per_s_left
  []
[]

[Modules]
  [FluidProperties]
    #[the_simple_fluid]
    #  type = SimpleFluidProperties
    #  bulk_modulus = 2E9
    #  viscosity = 1.0E-3
    #  density0 = 1000.0
    #[]
    [true_water97]
      type = Water97FluidProperties    # IAPWS-IF97
    []
    [tabulated_water95]                # tabulation is the only way to make this effective
      type = TabulatedFluidProperties  # the range 273.15 K <= T <= 1073.15 K for p <= 100 MPa should be OK [e.g. 800ºC at 10 km depth]
      # interpolated_properties = 'density enthalpy internal_energy viscosity k cp cv entropy'
      fp = true_water97
      fluid_property_file = water_IAPWS95.csv
    []
  []
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability] # has to range from 14-14 in the center to 1e-17 towards the sides
    type = PorousFlowPermeabilityTensorFromVar
    perm = perm_auxvar
  []
  [internal_energy]                # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    density = 3000.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 880.0        # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [lambda_mantle]                         # 'at_nodes=False' by default
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '1.5 0 0 0 1.5 0 0 0 3.3'
  []
[]

[Preconditioning]
  active = smp_lu_mumps
  [smp_lu_mumps]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  # start_time = 0             # [s] 0 year | env['time_vals'].max() from archaea03_exodus.e
  #end_time = 50               # [s] 1000 year [i.e. run for another 4000 kyr]
  end_time = 315576000000      # [s] 10000 yr
  #dtmax = 3.2E10              # [s] ~1000 year. advanced parameter. Maximum timestep size in an adaptive run. Default to 1E30   
  dtmax = 1.0E9                # 158 year
  nl_max_its = 25              # solver parameter. Max Nonlinear Iterations. Default to 50
  l_max_its = 500              # solver parameter. Max Linear Iterations. Default to 10000
  nl_abs_tol = 1E-06           # solver parameter. Nonlinear absolute tolerance. Default to 1E-50
  nl_rel_tol = 1E-08           # solver parameter. Nonlinear Relative Tolerance. Default to 1E-08 
  scheme = 'implicit-euler'    # the default TimeIntegrator [also known as backwards Euler method]
  # fixed_point_algorithm = 'picard' # the fixed point algorithm to converge the sequence of problems
  # line_search = 'default'
  [TimeStepper]                # TimeStepper subsystem [block always nested within the Executioner block]
    type = IterationAdaptiveDT # adjust the timestep based on the number of iterations
    optimal_iterations = 6     # target number of nonlinear iterations [from 6 to 10 this does not seem to affect convergence]
    #dt = 1                     # ~[27 h]
    dt = 1024
    growth_factor = 2
    cutback_factor = 0.5
  []
[]

[VectorPostprocessors]         # sample along a vector of coordinates [a transect]
  [temp_vtransect_x0]
    type = LineValueSampler
    start_point = '0 -15000 0'
    end_point =   '0 -3000  0'
    num_points = 500
    sort_by = y
    variable = temperature
  []
[]

[Outputs]
  [checkpoint]
    type = Checkpoint
    num_files = 4
    interval = 5
  []
  [exodus]
    type = Exodus
    output_material_properties = false
    interval = 5                          # a bit more than once per year for the initial dt
    append_date = false                   # datetime of creation of output file
    #execute_on = 'linear nonlinear failed timestep_end'
    execute_on = 'timestep_end'
    #start_step = 33
    #end_step = 34
    #output_linear = true
    #output_nonlinear = true
  []
  [csv]
    type = CSV                            # output for postprocessors, vector postprocessors, and scalar variables
    time_data = true
  []
[]
