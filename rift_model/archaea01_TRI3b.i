# run as:
# archaea-opt -i archaea01_TRI3b.i --mesh-only
#
# mpiexec -np 4 archaea-opt -i archaea01_TRI3b.i # Num Local DOFs: 5380
# [output: archaea01_TRI3b_exodus.e]
#
# +++ purpose +++
# Test: domain inherited from rift2ridge2D
# - 'nopore' subdomains: only thermal conduction
# in this "b" alternative, the "unsat" domain takes the same parameters as the upper crust, 
# and an unsaturated TH simulation is attempted
# 

[Mesh]
  construct_side_list_from_node_list=true
  [gen]
    type = FileMeshGenerator
    file = archaea_000563_i_exomerge_TRI3.e # I could alway rename input to mesh_000000_exomerge.e if all the rest does not depend on timestep
    use_for_exodus_restart = true      # needed to be able to read the ICs from the file [see below]
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
  [pc] # // Add the capillary pressure UserObject
    type = PorousFlowCapillaryPressureVG
    block = 'mantle lcrust ucrust sediment unsat'
    alpha = 1E-6
    m = 0.6
  []
  [pp_KT_advectiveflux_onecomp_userobj] # // add Advective Flux calculator UserObjects
    type = PorousFlowAdvectiveFluxCalculatorUnsaturated
    #type = PorousFlowAdvectiveFluxCalculatorUnsaturatedMultiComponent   # multi-phase, multi-comp one kernel for each fluid component
    block = 'mantle lcrust ucrust sediment unsat'
    phase = 0
    #fluid_component = 0
    flux_limiter_type = VanLeer
    multiply_by_density = true # default. Implies the advective flux is multiplied by density, so it is a mass flux
  []
  [heat_KT_advectiveflux_userobj]
    type = PorousFlowAdvectiveFluxCalculatorUnsaturatedHeat
    block = 'mantle lcrust ucrust sediment unsat'
    flux_limiter_type = VanLeer
    multiply_by_density = true # default. Implies the advective flux is multiplied by density, so it is a mass flux
  []
  [pp_solution]
    type = SolutionUserObject
    execute_on = INITIAL
    mesh = archaea_000563_i_exomerge_TRI3.e
    timestep = LATEST
    system_variables = P
  []
  [temp_solution]
    type = SolutionUserObject
    execute_on = INITIAL
    mesh = archaea_000563_i_exomerge_TRI3.e
    timestep = LATEST
    system_variables = T
  []
[]

[Variables]
  [porepressure]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
    initial_from_file_var = P # this is handled by CopyNodalVarsAction
    initial_from_file_timestep = LATEST # LATEST [note the ICs block must not appear, or would overrride this]
  []
   [temperature]            # [K] parsed initial condition in the ICs block below
    # exists over the entire mesh: governed by PorousFlow on lower block ['sat'], and diffusion on top block ['unsat'] 
    family = LAGRANGE
    order = FIRST
    initial_from_file_var = T # this is handled by CopyNodalVarsAction
    initial_from_file_timestep = LATEST # LATEST [note the ICs block must not appear, or would overrride this]
    scaling = 1E-08         # this variable scaling hopefully brings the residual R_T to the same order of magnitude that the Pressure residual R_P 
  []
[]

[AuxVariables]
  [darcy_vel_x]
    block = 'mantle lcrust ucrust sediment unsat'
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_y]
    block = 'mantle lcrust ucrust sediment unsat'
    family = MONOMIAL
    order = CONSTANT
  []
  [temperature_ini]
    family = LAGRANGE
    order = FIRST
    initial_from_file_var = T # this is handled by CopyNodalVarsAction
    initial_from_file_timestep = LATEST #
  []
  [mass_flux_top_submarine_sink]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [mass_flux_top_submarine_source]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [mass_flux_top_subaerial_sink]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [mass_flux_top_subaerial_source]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_adv_top_submarine_sink]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_adv_top_submarine_source]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_adv_top_subaerial_sink]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_adv_top_subaerial_source]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_dif_top_submarine]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_dif_top_subaerial]
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [ptop_auxvar] # used for top porepressure boundary conditions [from input file] - assumes the unsatured zone adds a pressure of 1atm/100m
    block = 'mantle lcrust ucrust sediment unsat'
    family = LAGRANGE
    order = FIRST
  []
  [perm_auxvar_el] # It seems this can be also a LAGRANGE, which is worth to look at in in the case we need to resort to time-constant parameters
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = permeability_el
  []
  # auxiliary variables for constant input materials
  # with 'initial_from_file_var' variables are initialized irrespective on whether thery are eventually used
  [porosity_auxvar_el]
    family = MONOMIAL
    order = CONSTANT
    initial_from_file_var = porosity_el # warning: this is defined globally, but for block 0. Let'see it works as input
  []
  #[vol_strain_auxv_el]
  #  family = MONOMIAL
  #  order = CONSTANT
  #  initial_from_file_var = vol_strain_el
  #[]
  #[foo_temperature_auxvar_el] # great this works, while the input file is TRI6!
  #  family = LAGRANGE
  #  order = FIRST
  #  initial_from_file_var = T
  #[]
  [sat_auxvar]
    block = 'mantle lcrust ucrust sediment unsat'
    family = MONOMIAL
    order = CONSTANT
  []  
[]

[AuxKernels]
  [darcy_vel_x_kernel]
    type = PorousFlowDarcyVelocityComponent
    block = 'mantle lcrust ucrust sediment unsat'
    component = x
    variable = darcy_vel_x
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
  []
  [darcy_vel_y_kernel]
    type = PorousFlowDarcyVelocityComponent
    block = 'mantle lcrust ucrust sediment unsat'
    component = y
    variable = darcy_vel_y
    fluid_phase = 0                             # OPTIONAL for single-phase
    execute_on = TIMESTEP_END
  []
  [saturation_kernel]
    type = PorousFlowPropertyAux
    block = 'mantle lcrust ucrust sediment unsat'
    variable = sat_auxvar
    property = saturation
  []
[]

[Kernels] # as added by the action [PorousFlowFullySaturated]
  [pp_time_derivative]                         # one kernel for each fluid component
    type = PorousFlowMassTimeDerivative        # these kernels lump the fluid-component mass to the nodes to ensure superior numerical stabilization
    block = 'mantle lcrust ucrust sediment unsat'
    variable = porepressure
    fluid_component = 0 
  []
  #[pp_advectiveflux_kernel]  # if (_stabilization == StabilizationEnum::Full)
  #  type = PorousFlowAdvectiveFlux # full upwinded advective flux of the fluid component
  #  block = 'mantle lcrust ucrust sediment unsat'
  #  variable = porepressure
  #[]
  [pp_KT_advectiveflux_kernel] # if (_stabilization == StabilizationEnum::KT)
    type = PorousFlowFluxLimitedTVDAdvection    # one kernel for each fluid component
    block = 'mantle lcrust ucrust sediment unsat'
    advective_flux_calculator = pp_KT_advectiveflux_onecomp_userobj  # PorousFlowAdvectiveFluxCalculator UserObjectName
    variable = porepressure
  []
  # Heat equation
  [heat_time_derivative]
    type = PorousFlowEnergyTimeDerivative      # this kernel lumps the heat energy-density to the nodes to ensure superior numerical stabilization
    block = 'mantle lcrust ucrust sediment unsat'
    variable = temperature
  []
  #[heat_advectiveflux_kernel] # if (_stabilization == StabilizationEnum::Full)
  #  type = PorousFlowHeatAdvection
  #  block = 'mantle lcrust ucrust sediment unsat'
  #  variable = temperature
  #[]
  [heat_KT_advectiveflux_kernel] # if (_stabilization == StabilizationEnum::KT)
    type = PorousFlowFluxLimitedTVDAdvection
    block = 'mantle lcrust ucrust sediment unsat'
    variable = temperature
    advective_flux_calculator = heat_KT_advectiveflux_userobj         # PorousFlowAdvectiveFluxCalculator UserObjectName
  []
  [heat_conduction]
    type = PorousFlowHeatConduction            # weak form of $-\div (\lambda \grad T)$
    block = 'mantle lcrust ucrust sediment unsat'
    variable = temperature
  []
  # kernels for subDomains with only thermal conduction
  [heat_time_derivative_nopore]
    type = SpecificHeatConductionTimeDerivative
    block = 'nopore'
    variable = temperature
    #lumping = true | lumping here prevents convergence! Don't know why
    density = density_nopore
    specific_heat = specific_heat_nopore
  []
  [heat_conduction_nopore]
    type = HeatConduction
    block = 'nopore'
    variable = temperature
    diffusion_coefficient = thermal_conductivity_nopore
  []
[]

[Functions]
  [porepressure_f]
    type = SolutionFunction
    from_variable = P
    solution = pp_solution
  []
  [temperature_f]
    type = SolutionFunction
    from_variable = T
    solution = temp_solution
  []
[]

[ICs]
  [ptop_auxvar_IC] # 
    block = 'mantle lcrust ucrust sediment unsat'
    type = FunctionIC
    variable = ptop_auxvar
    function = porepressure_f
  []
[]

[BCs]
  [Ptop_submarine_sink]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = top_submarine               # nodes touching the ocean bottom
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure: 1 atm + water column
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    use_mobility = true
    use_relperm = true
    flux_function = 100
    fluid_phase = 0
    save_in = mass_flux_top_submarine_sink
  []
  [Ptop_submarine_source]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = top_submarine
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 0'
    multipliers = '-1e9 0'
    use_mobility = true
    use_relperm = false
    flux_function = 100
    fluid_phase = 0
    save_in = mass_flux_top_submarine_source
  []
  [Ptop_subaerial_sink]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = top_subaerial               # nodes top-bounding the subDomain for porous flow simulation 
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure: 1 atm in subaerial nodeset
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    use_mobility = true
    use_relperm = true
    flux_function = 100                     # similar to Ptop_submarine_sink to prevent sharp fronts
    fluid_phase = 0
    save_in = mass_flux_top_subaerial_sink
  []
  [Ptop_subaerial_source]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = top_subaerial
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 0'
    multipliers = '-1e9 0'
    use_mobility = true
    use_relperm = false
    flux_function = 0.001                   # to be calibrated: should allow for emptying of porespace
    fluid_phase = 0
    save_in = mass_flux_top_subaerial_source
  []
  [Ttop_submarine_adv_sink]
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    boundary = top_submarine               # nodes touching the ocean bottom
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure: 1 atm + water column
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    use_mobility = true
    use_enthalpy = true
    use_relperm = true
    flux_function = 100
    fluid_phase = 0
    save_in = heat_flux_adv_top_submarine_sink
  []
  [Ttop_submarine_adv_source]
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    boundary = top_submarine               # nodes touching the ocean bottom
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure: 1 atm + water column
    pt_vals = '-1e9 0'
    multipliers = '-1e9 0'
    use_mobility = true
    use_enthalpy = true
    use_relperm = false
    flux_function = 100
    fluid_phase = 0
    save_in = heat_flux_adv_top_submarine_source
  []
  [Ttop_subaerial_sink]
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    boundary = top_subaerial               # nodes top-bounding the subDomain for porous flow simulation 
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure: 1 atm + water column
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    use_mobility = true
    use_enthalpy = true
    use_relperm = true
    flux_function = 100 # similar to Ptop_submarine_sink to avoid shocks
    fluid_phase = 0
    save_in = heat_flux_adv_top_subaerial_sink
  []
  [Ttop_subaerial_source]
    type = PorousFlowPiecewiseLinearSink  # 
    variable = temperature
    boundary = top_subaerial               # nodes top-bounding the subDomain for porous flow simulation 
    PT_shift = ptop_auxvar                 # [Pa] reference environmental pressure: 1 atm + water column
    pt_vals = '-1e9 0'
    multipliers = '-1e9 0'
    use_mobility = true
    use_enthalpy = true
    use_relperm = false
    flux_function = 0.001 # to be calibrated: should allow for emptying of porespace
    fluid_phase = 0
    save_in = heat_flux_adv_top_subaerial_source
  []
  [Ttop_submarine_dif] # diffusive heat flux = - lambda * grad_T
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = 278.15                            # [K] 5ºC
    boundary = top_submarine
    pt_vals     = '-1000 0 1000'
    multipliers = '-10000 0 1000'
    use_thermal_conductivity = true 
    flux_function = 0.01
    save_in = heat_flux_dif_top_submarine
  []
  [Ttop_subaerial_dif]
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = 278.15
    boundary = top_subaerial
    pt_vals     = '-1000 0 1000'
    multipliers = '-1000 0 1000'
    use_thermal_conductivity = true 
    flux_function = 1
    save_in = heat_flux_dif_top_subaerial
  []
  [Tbot]
    type = FunctionDirichletBC 
    variable = temperature
    boundary = bottom
    function = temperature_f
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
      # fluid_property_file = water97_tabulated_extrap_allsides.csv
      fluid_property_file = water_IAPWS95.csv
    []
  []
[]

[Materials]
  [materials_nopore]
    type = GenericConstantMaterial
    block = nopore
    prop_names = 'thermal_conductivity_nopore specific_heat_nopore density_nopore' 
    prop_values = '3.3 1200.0 3360.0' # assume nopore is mantle material for the physical properties
  []
  #[materials_unsat]
  #  type = GenericConstantMaterial
  #  block = unsat
  #  prop_names = 'thermal_conductivity_unsat specific_heat_unsat density_unsat' 
  #  prop_values = '2.3 1200.0 2700.0' # assume unsat is upper crust material for the physical properties
  #[]
  # materials added by the Action PorousFlowUnsaturated:
  [saturation_calculator]
    type = PorousFlow1PhaseP
    block = 'mantle lcrust ucrust sediment unsat'
    porepressure = porepressure
    capillary_pressure = pc
  []
  [temperature_material]
    type = PorousFlowTemperature            # at_nodes=false by default
    block = 'mantle lcrust ucrust sediment unsat' # despite the variable temperature exists in all blocks, this material only applies to the PorousFlow subdomain
    temperature = temperature
  []
  [massfrac] # although it is trivially 1.0 for this single-phase water-only problem                               
    type = PorousFlowMassFraction           # at_nodes=false by default
    block = 'mantle lcrust ucrust sediment unsat'
  []
  [fluid_eos]
    type = PorousFlowSingleComponentFluid      # see documentation for this material regarding the choice of units
    block = 'mantle lcrust ucrust sediment unsat'
    fp = the_simple_fluid                     # this Material is at_nodes=false by default
    #fp = tabulated_water95
    phase = 0
  []
  [effective_fluid_pressure]                   # create effective fluid pressure [requested by PorousFlowPorosity even it has not mechanical coupling]
    type = PorousFlowEffectiveFluidPressure
    block = 'mantle lcrust ucrust sediment unsat'
  []
  [nearest_qp]                                 # according to the example multiblock.i, this seems not needed
    type = PorousFlowNearestQp                 # atNodes=false by default
    block = 'mantle lcrust ucrust sediment unsat'
  []
  # materials not explicitly added by the Action PorousFlowUnsaturated:
  [porosity]
    type = PorousFlowPorosityConst
    block = 'mantle lcrust ucrust sediment unsat'
    porosity = porosity_auxvar_el
  []
  [permeability]
    type = PorousFlowPermeabilityTensorFromVar
    block = 'mantle lcrust ucrust sediment unsat'
    perm = perm_auxvar_el # warning: this makes the model halt if input variable is nodal [LAGRANGE] 
  []
  [internal_energy_mantle]                # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = 'mantle'
    density = 3360.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0       # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [internal_energy_lcrust]                       # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = 'lcrust'
    density = 2850.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0        # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [internal_energy_ucrust]                       # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = 'ucrust'
    density = 2700.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0        # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [internal_energy_sediment]                # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = 'sediment'
    density = 2450.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0        # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [internal_energy_unsat]                # 'at_nodes=True' by default | internal_energy [J.K-1.m-3] = density*specific_heat_capacity
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = 'unsat'
    density = 2700.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0        # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  # ocrust subDomain does not exist in this simulation. This test indicates the lack of the indicated subdomain leads to a crash by this input lock 
  #[internal_energy_ocrust]                
  #  type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
  #  block = 'ocrust'
  #  density = 2990.0                      # [kg.m-3] density of rock grains
  #  specific_heat_capacity = 1200.0        # [J.kg-1.K-1] specific heat capacity of rock grains
  #[]
  [lambda_mantle]                         # 'at_nodes=False' by default
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = 'mantle'
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '3.3 0 0 0 3.3 0 0 0 3.3'
  []
  [lambda_lcrust]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = 'lcrust'
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '2.50 0 0 0 2.50 0 0 0 2.50'
  []
  [lambda_ucrust]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = 'ucrust'
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '2.30 0 0 0 2.30 0 0 0 2.30'
  []
  [lambda_sediment]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = 'sediment'
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '3.00 0 0 0 3.00 0 0 0 3.00'
  []
  [lambda_unsat]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = 'unsat'
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '2.30 0 0 0 2.30 0 0 0 2.30'
  []
  [relperm] # required by PorousFlowDarcyVelocityComponent AuxKernels
    type = PorousFlowRelativePermeabilityCorey       # atNodes=false by default
    block = 'mantle lcrust ucrust sediment unsat'    #
    n = 3           # Corey exponent of the phase
    s_res = 0.1     # residual saturation of the phase
    sum_s_res = 0.1 # [>= s_res]
    phase = 0
  []
  #[vol_strain_mat] # defined by Javier [this is actually not used for constant input porosity here]
  #  type = PorousFlowVolumetricStrainConst
  #  block = 'mantle lcrust ucrust sediment unsat'    #
  #  vol_strain_auxvar = vol_strain_auxv_el
  #[]
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
  #end_time = 2                # [s] 1000 year [i.e. run for another 4000 kyr]
  end_time = 315576000000     # [s] 10 kyr
  #end_time =  315576000000     # [s] 0.1 Myr
  dtmax = 1e9                  # [s] 31.7 year. Maximum timestep size in an adaptive run. Default to 1E30   
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
    end_point =   '0 -5800  0' # in general this range should be automated
    num_points = 500
    sort_by = y
    variable = temperature
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
    output_material_properties = false
    interval = 5                          # a bit more than once per year for the initial dt
    append_date = false                   # datetime of creation of output file
    execute_on = 'timestep_end'
    #execute_on = 'nonlinear failed timestep_end'
    #output_nonlinear = true
  []
  [csv]
    type = CSV # output for postprocessors, vector postprocessors, and scalar variables
    time_data = true
  []
[]
