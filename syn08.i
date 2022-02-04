# mpiexec -np 4 ../../archaea-opt -i syn08.i
# construction of the problem in syn07.i but without the use of actions
# as a comparison test with 07.i, both {MPI ppn:4, grad_T_y = 600/150000, simple_fluid}
# 07.i: [376.1 s]
# 08.i: [380.6 s]


[Mesh]
  [2D]
    type = FileMeshGenerator
    file = syn00_in.e
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 -9.81 0'              # make global, so no need to add it to the userobject requiring it
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator 
    porous_flow_vars = 'porepressure temperature'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pp_KT_advectiveflux_onecomp_userobj]
    type = PorousFlowAdvectiveFluxCalculatorSaturated       # [1-phase, 1-comp, fully saturated] one kernel for each fluid component
    # block                                                 # OPTIONAL, list of blocks that this object will be applied to
    # fe_family   # needed only if the variables in the dictator have different families
    # fe_order    # " " " " " "
    flux_limiter_type = superbee
    multiply_by_density = true # default. Implies the advective flux is multiplied by density, so it is a mass flux
    # block
  []
  [heat_KT_advectiveflux_userobj]
    type = PorousFlowAdvectiveFluxCalculatorSaturatedHeat
    # block                                                 # OPTIONAL, list of blocks that this object will be applied to
    # fe_family   # needed only if the variables in the dictator have different families
    # fe_order    # " " " " " "
    flux_limiter_type = superbee
    multiply_by_density = true # default. Implies the advective flux is multiplied by density, so it is a mass flux.
    # block
  []
[]


[Variables]
  [porepressure]
    # active = 'ucrust lcrust'              # blocks where this variable exists
    family = LAGRANGE
    order = SECOND   
  []
   [temperature]                     # [K] parsed initial condition in the ICs block below
    family = LAGRANGE
    order = SECOND
    scaling = 1E-08                  # this variable scaling brings the residual R_T to the same order of magnitude that the Pressure residual R_P 
  []
[]

[AuxVariables]
  [darcy_vel_x]
    # active = 'ucrust lcrust'
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_y]
    # active = 'ucrust lcrust'
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
    # block = 
  []
  [darcy_vel_y_kernel]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_y
    fluid_phase = 0                             # OPTIONAL for single-phase
    # block = 
  []
[]

[Kernels] # as added by the action [PorousFlowFullySaturated]
  # note: kernels for chemistry not added in this simulation
  # fluid-flow equations
  [pp_time_derivative]                         # one kernel for each fluid component
    type = PorousFlowMassTimeDerivative        # these kernels lump the fluid-component mass to the nodes to ensure superior numerical stabilization
    variable = porepressure 
    # block
  []
  #[pp_mass_vol_expansion]                     # one kernel for each fluid component
  #  type = PorousFlowMassVolumetricExpansion  # only for transient simulations with mechanical coupling_type
  #  variable = porepressure
  #[]
  [pp_KT_advectiveflux_kernel]
    type = PorousFlowFluxLimitedTVDAdvection                              # one kernel for each fluid component
    advective_flux_calculator = pp_KT_advectiveflux_onecomp_userobj  # PorousFlowAdvectiveFluxCalculator UserObjectName
    variable = porepressure
    # block                                                           # OPTIONAL, list of blocks that this object will be applied to
  []
  # Heat equation
  [heat_time_derivative]
    type = PorousFlowEnergyTimeDerivative      # this kernel lumps the heat energy-density to the nodes to ensure superior numerical stabilization
    variable = temperature
    # save_in = # name of auxiliary variable to save this Kernel residual contributions to.
  []
  [heat_conduction]
    type = PorousFlowHeatConduction
    variable = temperature
  []
  [heat_KT_advectiveflux_kernel]
    type = PorousFlowFluxLimitedTVDAdvection                              # one kernel for each fluid component
    advective_flux_calculator = heat_KT_advectiveflux_userobj         # PorousFlowAdvectiveFluxCalculator UserObjectName
    variable = temperature
  []
[]


[ICs]
  [porepressure_IC]
    type = FunctionIC
    variable = porepressure               # [Pa]
    function = '9.81*(-y)*1000'
  []
  [temperature_IC]
    type = FunctionIC
    variable = temperature
    #function = '273.15+5+(-y)*42/150000'     # => 320.15 = 273.15+5+42 at the bottom
    #function = '273.15+5+(-y)*300/150000'    # => 578.15 = 273.15+5+300 at the bottom
    function = '273.15+5+(-y)*600/150000'     # => 878.15 = 273.15+5+600 at the bottom
    #function = '273.15+5+(-y)*1300/150000'   # => 1578.15 = 273.15+5+1300 at the bottom
  []
[]

[BCs]
  [Ptop]
    type = DirichletBC                     # this is a NodalBC
    variable = porepressure
    value = 101325.0                       # [Pa] 1 atm
    boundary = top
  []
  [Ttop]
    type = DirichletBC
    variable = temperature
    value = 278.15                         #
    boundary = top 
  []
  [Tbot]
    type = DirichletBC
    variable = temperature
    #value = 320.15
    #value = 578.15                         # [K]    100ºC > ICs
    value = 878.15
    #value = 1578.15
    boundary = bottom 
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
    #[true_water97]
    #  type = Water97FluidProperties    # IAPWS-IF97
    #[]
    #[tabulated_water97]                # tabulation is the only way to make this effective
    #  type = TabulatedFluidProperties  # the range 273.15 K <= T <= 1073.15 K for p <= 100 MPa should be OK [e.g. 800ºC at 10 km depth]
      # interpolated_properties = 'density enthalpy internal_energy viscosity k cp cv entropy'
    #  fp = true_water97
    #  fluid_property_file = water97_tabulated_extrap.csv
    #[]
  []
[]

[Materials]
  [porosity_mantle]
    type = PorousFlowPorosity                  # 'at_nodes=False' by default
    block = mantle
    porosity_zero = 0.03
    thermal = true
    thermal_expansion_coeff = 3.96E-05
    reference_temperature = 273.15
    fluid = true
    solid_bulk = 110.E09
    biot_coefficient = 0.8
    biot_coefficient_prime = 0.8
    mechanical = false
  []
  [porosity_lcrust]
    type = PorousFlowPorosity
    block = lcrust
    porosity_zero = 0.05
    thermal = true
    thermal_expansion_coeff = 3.66E-05
    reference_temperature = 273.15
    fluid = true
    solid_bulk = 10.E09
    biot_coefficient = 0.8
    biot_coefficient_prime = 0.8
    mechanical = false
  []
  [porosity_ucrust]
    type = PorousFlowPorosity
    block = ucrust
    porosity_zero = 0.05
    thermal = true
    thermal_expansion_coeff = 3.38E-05
    reference_temperature = 273.15
    fluid = true
    solid_bulk = 10.E09
    biot_coefficient = 0.8
    biot_coefficient_prime = 0.8
    mechanical = false
  []
  [porosity_unsat]
    type = PorousFlowPorosity
    block = unsat
    porosity_zero = 0.05
    thermal = true
    thermal_expansion_coeff = 3.38E-05
    reference_temperature = 273.15
    fluid = true
    solid_bulk = 10.E09
    biot_coefficient = 0.8
    biot_coefficient_prime = 0.8
    mechanical = false
  []
  [permeability_mantle]
    type = PorousFlowPermeabilityKozenyCarman    # 'at_nodes=False' by default
    block = mantle
    n = 3
    m = 2
    phi0 = 0.03
    k0 = 1.0E-22
    poroperm_function = 'kozeny_carman_phi0'
  []
  [permeability_lcrust]
    type = PorousFlowPermeabilityKozenyCarman
    block = lcrust
    n = 3
    m = 2
    phi0 = 0.05
    k0 = 9.0E-16
    poroperm_function = 'kozeny_carman_phi0'
  []
  [permeability_ucrust]
    type = PorousFlowPermeabilityKozenyCarman
    block = ucrust
    n = 3
    m = 2
    phi0 = 0.05
    k0 = 1.0E-15
    poroperm_function = 'kozeny_carman_phi0'
  []
  [permeability_unsat]
    type = PorousFlowPermeabilityKozenyCarman
    block = unsat
    n = 3
    m = 2
    phi0 = 0.05
    k0 = 1.0E-22
    poroperm_function = 'kozeny_carman_phi0'
  []                
  [internal_energy_mantle]                # 'at_nodes=True' by default
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = mantle                        # internal_energy [J.K-1.m-2¡3] = density*specific_heat_capacity
    density = 3360.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0       # [J.kg-1.K-1] specific heat capacity of rock grains
  [] 
  [internal_energy_lcrust]
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = lcrust
    density = 2850.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0       # [J.kg-1.K-1] specific heat capacity of rock grains
  [] 
  [internal_energy_ucrust]
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = ucrust
    density = 2700.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0       # [J.kg-1.K-1] specific heat capacity of rock grains
  [] 
  [internal_energy_unsat]
    type = PorousFlowMatrixInternalEnergy # this Material calculated the internal energy of solid rock grains
    block = unsat
    density = 2700.0                      # [kg.m-3] density of rock grains
    specific_heat_capacity = 1200.0       # [J.kg-1.K-1] specific heat capacity of rock grains
  []
  [lambda_mantle]                         # 'at_nodes=False' by default
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = mantle
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '3.3 0 0 0 3.3 0 0 0 3.3'
  []
  [lambda_lcrust]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = lcrust
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '2.5 0 0 0 2.5 0 0 0 2.5'
  []
  [lambda_ucrust]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = ucrust
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '2.3 0 0 0 2.3 0 0 0 2.3'
  []
  [lambda_unsat]
    type = PorousFlowThermalConductivityFromPorosity # rock-fluid combined thermal conductivity by weighted sum of rock and fluid conductivities
    block = unsat
    lambda_f = '0.56 0 0 0 0.56 0 0 0 0.56'
    lambda_s = '2.3 0 0 0 2.3 0 0 0 2.3'
  []
  # material added by PorousFlowFullySaturated Action
  [porepressure_material]
    type = PorousFlow1PhaseFullySaturated # at_nodes=false by default
    porepressure = porepressure
  []
  [temperature_material]
    type = PorousFlowTemperature          # at_nodes=false by default
    temperature = temperature  
  []
  [massfrac]
    type = PorousFlowMassFraction         # at_nodes=false by default
    # mass_fraction_vars # nod needed when num_phases=num_component
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid      # see documentation for this material regarding the choice of units
    fp = the_simple_fluid                      # this Material is at_nodes=false by default
    phase = 0
  []
  [effective_fluid_pressure] # create effective fluid pressure [is requested by PorousFlowPorosity even it has not mechanical coupling]
    type = PorousFlowEffectiveFluidPressure
  []
  # type = PorousFlowVolumetricStrain          # only for mechanically-coupled model. But I should replace this by an input Auxiliary variable!
  [nearest_qp]
    # block # I am not sure about this!
    type = PorousFlowNearestQp                 # atNodes=false by default
  []
  [relperm]
    type = PorousFlowRelativePermeabilityConst # atNodes=false by default
    phase = 0
    kr = 1                                     # default, anyway
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
  end_time = 315576000000      # [s] 10000 year
  #end_time = 1E6              # [s] ~[12 days]
  dtmax = 3.2E10               # [s] ~1000 year. advanced parameter. Maximum timestep size in an adaptive run. Default to 1E30   
  nl_max_its = 25              # solver parameter. Max Nonlinear Iterations. Default to 50
  l_max_its = 100              # solver parameter. Max Linear Iterations. Default to 10000
  nl_abs_tol = 1E-06           # solver parameter. Nonlinear absolute tolerance. Default to 1E-50
  nl_rel_tol = 1E-08           # solver parameter. Nonlinear Relative Tolerance. Default to 1E-08 
  scheme = 'implicit-euler'    # the default TimeIntegrator
  [TimeStepper]                # TimeStepper subsystem [block always nested within the Executioner block]
    type = IterationAdaptiveDT # adjust the timestep based on the number of iterations
    optimal_iterations = 10    # target number of nonlinear iterations
    dt = 1E5                   # ~[27 h]
    growth_factor = 2
    cutback_factor = 0.5
  []
[]

[Outputs]
  [ou]
    type = Exodus
    #output_material_properties = true
  []
[]
