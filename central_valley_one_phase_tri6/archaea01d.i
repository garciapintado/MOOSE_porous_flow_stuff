# run as:
# porousflow-opt -i archaea01d.i --mesh-only
#
# mpiexec -np 6 porousflow-opt -i archaea01d.i 
#                                          
# simple synthetic setup of hydrothermal circulation model with central valley topography
# unsaturated TH simulation
# step 01 simply attempts to brings the model to equilibrium with a standard geothermal gradient 
#
# domain: xlim = [-40000,40000] m
#         ylim = [-15000,1000] m with central valley reaching y=-2000m at the center
#                                and maximum topography reaching y=1000 at the top corners
# ocean at y=0

[Mesh]
  construct_side_list_from_node_list=true
  [gen]
    type = FileMeshGenerator
    file = mesh_000000_exomerge.e
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 -9.81 0'
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure temperature'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]                                      # // Add the capillary pressure UserObject
    type = PorousFlowCapillaryPressureVG
    alpha = 1E-6
    m = 0.6
  []
  #[pp_KT_advectiveflux_onecomp_userobj]  # // add Advective Flux calculator UserObjects
  #  type = PorousFlowAdvectiveFluxCalculatorUnsaturated
  #  phase = 0
  #  flux_limiter_type = VanLeer
  #  multiply_by_density = true # default. Implies the advective flux is multiplied by density, so it is a mass flux
  #[]
  #[heat_KT_advectiveflux_userobj]
  #  type = PorousFlowAdvectiveFluxCalculatorUnsaturatedHeat
  #  flux_limiter_type = VanLeer
  #  multiply_by_density = true # default. Implies the advective flux is multiplied by density, so it is a heat flux
  #[]
[]

[Variables]
  [porepressure]
    family = LAGRANGE
    order = SECOND
  []
   [temperature]            # [K] parsed initial condition in the ICs below
    family = LAGRANGE
    order = SECOND
    scaling = 1E-08
  []
[]

[AuxVariables]
  [darcy_vel_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [mass_flux_top_submarine]
    family = LAGRANGE
    order = SECOND
  []
  [mass_flux_top_subaerial_source]
    family = LAGRANGE
    order = SECOND
  []
  [heat_flux_adv_top_submarine]
    family = LAGRANGE
    order = SECOND
  []
  [heat_flux_adv_top_subaerial_source]
    family = LAGRANGE
    order = SECOND
  []
  [heat_flux_dif_top]
    family = LAGRANGE
    order = SECOND
  []
  [ptop_auxvar]
    family = LAGRANGE
    order = SECOND
  []
  [perm_auxvar]
    family = MONOMIAL      # I saw that MONOMIAL makes the model more stable. Don't know why.
    order = CONSTANT
  []
  [sat_auxvar]
    family = MONOMIAL
    order = CONSTANT
  [] 
[]

[AuxKernels]
  [darcy_vel_x_kernel]
    type = PorousFlowDarcyVelocityComponent
    component = x
    variable = darcy_vel_x
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
  []
  [darcy_vel_y_kernel]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_y
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
  []
  [saturation_kernel]
    type = PorousFlowPropertyAux
    variable = sat_auxvar
    property = saturation
  []
[]

[Kernels] # as added by the action [PorousFlowFullySaturated]
  [pp_time_derivative]                         # one kernel for each fluid component
    type = PorousFlowMassTimeDerivative        # these kernels lump the fluid-component mass to the nodes to ensure superior numerical stabilization
    variable = porepressure
    fluid_component = 0 
  []
  [pp_advectiveflux_kernel]  # if (_stabilization == StabilizationEnum::Full)
    type = PorousFlowAdvectiveFlux # full upwinded advective flux of the fluid component
    variable = porepressure
  []
  #[pp_KT_advectiveflux_kernel] # if (_stabilization == StabilizationEnum::KT)
  #  type = PorousFlowFluxLimitedTVDAdvection    # one kernel for each fluid component
  #  advective_flux_calculator = pp_KT_advectiveflux_onecomp_userobj  # PorousFlowAdvectiveFluxCalculator UserObjectName
  #  variable = porepressure
  #[]
  # Heat equation
  [heat_time_derivative]
    type = PorousFlowEnergyTimeDerivative      # this kernel lumps the heat energy-density to the nodes to ensure superior numerical stabilization
    variable = temperature
  []
  [heat_advectiveflux_kernel] # if (_stabilization == StabilizationEnum::Full)
    type = PorousFlowHeatAdvection
    variable = temperature
  []
  #[heat_KT_advectiveflux_kernel] # if (_stabilization == StabilizationEnum::KT)
  #  type = PorousFlowFluxLimitedTVDAdvection
  #  variable = temperature
  #  advective_flux_calculator = heat_KT_advectiveflux_userobj         # PorousFlowAdvectiveFluxCalculator UserObjectName
  #[]
  [heat_conduction]
    type = PorousFlowHeatConduction            # weak form of $-\div (\lambda \grad T)$
    variable = temperature
  []
[]

[Functions]
  [temperature00] # temperature ICs
    type = ParsedFunction
    value = '273.15 + 5 + (150-5)*exp(-3*(y-ybot)/rate)'             # 150ºC at the bottom and exponential decay toward the surface
    vars = 'ybot   rate'
    vals = '-15000 6000'
  []
  [porepressure00]
    type = ParsedFunction
    value = 'if(y > sea_level, atm, atm + 9.81*(-y + sea_level)*1000)'    # top corners [y=1000m] at atmospheric pressure
    vars = 'sea_level atm'
    vals = '0.0 0.0'
  []
  [permeability00]
    type = ParsedFunction 
    value = 'if(y > sea_level, 1e-15, vmin + (vmax-vmin) * exp(-((x-xc)/xsd)^2))'
    vars = 'sea_level xc  vmin    vmax    xsd'
    vals = '0.0      0.0  1.0E-17 1.0E-14 5000'
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
    function = permeability00
  []
[]

[BCs]
  [Ptop_submarine] # sink & source
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = top_submarine
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 1e9'
    multipliers = '-1e9 1e9'
    use_mobility = true
    use_relperm = true
    flux_function = 100
    fluid_phase = 0
    save_in = mass_flux_top_submarine
  []
  [Ttop_submarine_adv] # sink & source: advective heat flux = enthalpy * darcy Flux
    type = PorousFlowPiecewiseLinearSink         # this is a NodalBC
    variable = temperature
    PT_shift = ptop_auxvar
    boundary = top_submarine
    pt_vals = '-1e9 1e9'
    multipliers = '-1e9 1e9'
    use_mobility = true
    use_enthalpy = true
    use_relperm = true
    flux_function = 100 
    fluid_phase = 0
    save_in = heat_flux_adv_top_submarine
  []
  [Ptop_subaerial_source]                    # small continuous rainfall rate
    # 100 [mm.yr-1] => 3.2e-06 [kg.m-2.s-1] 
    type = PorousFlowSink
    variable = porepressure
    boundary = top_subaerial
    #PT_shift = ptop_auxvar
    #pt_vals = '-1e9 0'
    #multipliers = '-1e9 0'
    use_mobility = true
    use_relperm = false
    flux_function = -3.2e-04                   # to be calibrated: should allow for emptying of porespace
    fluid_phase = 0
    save_in = mass_flux_top_subaerial_source
  []
  [Ttop_adv_subaerial_source]                    # advective heat flux = enthalpy * darcy Flux
    type = PorousFlowSink                        # this is a NodalBC
    variable = temperature
    boundary = top_subaerial
    # PT_shift = ptop_auxvar
    # pt_vals = '-1e9 0'
    # multipliers = '-1e9 0'
    use_mobility = true
    use_enthalpy = true
    use_relperm = false
    flux_function = -3.2e-04 
    fluid_phase = 0
    save_in = heat_flux_adv_top_subaerial_source
  []
   [Ttop_diff] # diffusive heat flux = - lambda * grad_T
    type = PorousFlowPiecewiseLinearSink          # this is a NodalBC
    variable = temperature
    PT_shift = 278.15                             # [K]    5ºC
    boundary = top
    pt_vals     = '-1000 0 1000'                  # x coordinates defining g
    multipliers = '-10000 0 1000'                 # y coordinates defining g
    use_thermal_conductivity = true 
    flux_function = 0.01
    save_in = heat_flux_dif_top 
   []
  [Tbot]
    type = FunctionDirichletBC 
    variable = temperature
    boundary = bottom
    function = temperature00
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2E9
      viscosity = 1.0E-3
      density0 = 1000.0
    []
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
  [saturation_calculator]
    type = PorousFlow1PhaseP
    porepressure = porepressure
    capillary_pressure = pc
  []
  [temperature_material]
    type = PorousFlowTemperature            # at_nodes=false by default
    temperature = temperature
  []
  [massfrac] # although it is trivially 1.0 for this single-phase water-only problem                               
    type = PorousFlowMassFraction           # at_nodes=false by default
  []
   [fluid_eos]
    type = PorousFlowSingleComponentFluid # see documentation for this material regarding the choice of units
    fp = the_simple_fluid                          # this Material is at_nodes=false by default
    #fp = tabulated_water95
    phase = 0
  []
  [effective_fluid_pressure]                   # create effective fluid pressure [requested by PorousFlowPorosity even it has not mechanical coupling]
    type = PorousFlowEffectiveFluidPressure
  []
  # materials not explicitly added by the Action PorousFlowUnsaturated:
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability] # has to range from 14-14 in the center to 1e-17 towards the sides
    type = PorousFlowPermeabilityTensorFromVar
    perm = perm_auxvar
  []
   [internal_energy]                      # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    density = 3000.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1000.0       # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [lambda]                                           # 'at_nodes=False' by default
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '1.5 0 0 0 1.5 0 0 0 3.3'
  []
  [relperm] # required by PorousFlowDarcyVelocityComponent AuxKernels
    type = PorousFlowRelativePermeabilityCorey       # atNodes=false by default
     n = 3           # Corey exponent of the phase
    s_res = 0.1     # residual saturation of the phase
    sum_s_res = 0.1 # [>= s_res]
    phase = 0
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
  end_time = 3600              # [s] 1000 year [i.e. run for another 4000 kyr]
  #end_time = 315576000000     # [s] 10000 yr - one timestep of our mechanical model
  dtmax = 5e9                  # [s] ~31.7 year. advanced parameter. Maximum timestep size in an adaptive run. Default to 1E30   
  nl_max_its = 25              # solver parameter. Max Nonlinear Iterations. Default to 50
  l_max_its = 500              # solver parameter. Max Linear Iterations. Default to 10000
  nl_abs_tol = 1.1E-06         # solver parameter. Nonlinear absolute tolerance. Default to 1E-50
  nl_rel_tol = 1E-08           # solver parameter. Nonlinear Relative Tolerance. Default to 1E-08 
  scheme = 'implicit-euler'    # the default TimeIntegrator [also known as backwards Euler method]
  # fixed_point_algorithm = 'picard' # the fixed point algorithm to converge the sequence of problems
  # line_search = 'default'
  [TimeStepper]                # TimeStepper subsystem [block always nested within the Executioner block]
    type = IterationAdaptiveDT # adjust the timestep based on the number of iterations
    optimal_iterations = 6     # target number of nonlinear iterations [from 6 to 10 this does not seem to affect convergence]
    dt = 1                     # ~[27 h]
    growth_factor = 2
    cutback_factor = 0.5
  []
[]

[VectorPostprocessors]         # sample along a vector of coordinates [a transect]
  [temp_vtransect_x0]
    type = LineValueSampler
    start_point = '0 -15000 0'
    end_point =   '0 -2000  0'
    num_points = 500
    sort_by = y
    variable = temperature
  []
  [temp_vtransect_left]
    type = LineValueSampler
    start_point = '-39000 -15000 0'
    end_point =   '-39000   980  0'
    num_points = 500
    sort_by = y
    variable = porepressure
  []
[]

[Outputs]
  #[checkpoint]
  #  type = Checkpoint
  #  num_files = 4
  #  interval = 5
  #[]
  [exodus]
    type = Exodus
    output_material_properties = true
    interval = 1                          # a bit more than once per year for the initial dt
    append_date = false                   # datetime of creation of output file
  []
  [csv]
    type = CSV # output for postprocessors, vector postprocessors, and scalar variables
    time_data = true
  []
[]
