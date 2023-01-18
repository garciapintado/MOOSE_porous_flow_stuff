# gas injection
# non-isothermal miscible 2-phase 2-component model with UserObject PorousFlowWaterNCG
# - 2D domain: rectangular column [width:25m, height:200m]
# - no mechanics
# - standard gravity
# - gas injection increases linearly with time to allow for initial stabilization 

ybot = -200
xmin = 0
xmax = 25
density0_phase0 = 1000

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${xmin}
    xmax = ${xmax}
    ymin = ${ybot}
    ymax =  0
    elem_type = TRI3
    nx = 50
    ny = 100
  []
  [caps]
    type = SubdomainBoundingBoxGenerator
    input = 'gen'
    block_id = 1
    bottom_left = '-50 -50 0'
    top_right = '50 0 0'
  []
  [injection_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > -150 & y < -125 & x < ${xmin}+1e-05'
    included_subdomains = 0
    new_sideset_name = 'injection_area'
    input = 'caps'
  []
  [rename]
    type = RenameBlockGenerator
    old_block = '0 1'
    new_block = 'aquifer caps'
    input = 'injection_area'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  biot_coefficient = 1.0
  gravity = '0 -9.81 0'
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi temperature'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 1E-6
    m = 0.6
  []
  [fs]
    type = PorousFlowWaterNCG
    water_fp = true_water
    gas_fp = true_methane
    capillary_pressure = pc
  []
[]

[Variables]
  [pgas]
  []
  [zi]
    initial_condition = 0
  []
  [temperature]
    scaling = 1E-5
  []
[]

[Kernels]
  [mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pgas
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pgas
  []
  [mass_methane_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = zi
  []
  [flux_methane]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = zi
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = temperature
  []
  [heat_advection]
    type = PorousFlowHeatAdvection
    variable = temperature
  []
  [heat_conduction]
    type = PorousFlowHeatConduction
    variable = temperature
  []
[]

[AuxVariables]
  [darcy_vel_phase0_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_phase0_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_phase1_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [darcy_vel_phase1_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [pwater_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [pgas_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
  []
  [sgas]
    family = MONOMIAL
    order = CONSTANT
  []
  [mfrac_0_0]   # _\beta_\k
    order = CONSTANT
    family = MONOMIAL
  []
  [mfrac_0_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [mfrac_1_0]
    order = CONSTANT
    family = MONOMIAL
  []
  [mfrac_1_1]
    order = CONSTANT
    family = MONOMIAL
  []
  [density_0_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [density_1_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [enthalpy_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [enthalpy_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [relperm_0_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [relperm_1_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [porosity_aux]
    family = MONOMIAL
    order = CONSTANT
  []
  [mass_water_flux_top]
    family = LAGRANGE
    order = FIRST
  []
  [mass_methane_flux_top]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_water_adv_top]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_gas_adv_top]
    family = LAGRANGE
    order = FIRST
  []
  [heat_flux_dif_top]
    family = LAGRANGE
    order = FIRST
  []
  [ptop_auxvar]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxKernels]
  [darcy_vel_phase0_x]
    type = PorousFlowDarcyVelocityComponent
    component = x
    variable = darcy_vel_phase0_x
    fluid_phase = 0
  []
  [darcy_vel_phase0_y]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_phase0_y
    fluid_phase = 0
  []
  [darcy_vel_phase1_x]
    type = PorousFlowDarcyVelocityComponent
    component = x
    variable = darcy_vel_phase1_x
    fluid_phase = 1
  []
  [darcy_vel_phase1_y]
    type = PorousFlowDarcyVelocityComponent
    component = y
    variable = darcy_vel_phase1_y
    fluid_phase = 1
  []
  [pwater_aux]
    type = PorousFlowPropertyAux
    variable = pwater_aux
    property = pressure
    phase = 0
  []
  [pgas_aux]
    type = PorousFlowPropertyAux
    variable = pgas_aux
    property = pressure
    phase = 1
  []
  [swater]
    type = PorousFlowPropertyAux
    variable = swater
    property = saturation
    phase = 0
  []
  [sgas]
    type = PorousFlowPropertyAux
    variable = sgas
    property = saturation
    phase = 1
  []
  [mfrac_0_0]
    type = PorousFlowPropertyAux
    variable = mfrac_0_0
    property = mass_fraction
    phase = 0
    fluid_component = 0
  []
  [mfrac_0_1]
    type = PorousFlowPropertyAux
    variable = mfrac_0_1
    property = mass_fraction
    phase = 0
    fluid_component = 1
  []
  [mfrac_1_0]
    type = PorousFlowPropertyAux
    variable = mfrac_1_0
    property = mass_fraction
    phase = 1
    fluid_component = 0
  []
  [mfrac_1_1]
    type = PorousFlowPropertyAux
    variable = mfrac_1_1
    property = mass_fraction
    phase = 1
    fluid_component = 1
  []
  [density_0]
    type = PorousFlowPropertyAux
    variable = density_0_aux
    property = density
    phase = 0
  []
  [density_1]
    type = PorousFlowPropertyAux
    variable = density_1_aux
    property = density
    phase = 1
  []
  [viscosity_water]
    type = PorousFlowPropertyAux
    variable = viscosity_water
    property = viscosity
    phase = 0
  []
  [viscosity_gas]
    type = PorousFlowPropertyAux
    variable = viscosity_gas
    property = viscosity
    phase = 1
  []
  [enthalpy_water]
    type = PorousFlowPropertyAux
    variable = enthalpy_water
    property = enthalpy
    phase = 0
  []
  [enthalpy_gas]
    type = PorousFlowPropertyAux
    variable = enthalpy_gas
    property = enthalpy
    phase = 1
  []
  [relperm_0]
    type = PorousFlowPropertyAux
    variable = relperm_0_aux
    property = relperm
    phase = 0
  []
  [relperm_1]
    type = PorousFlowPropertyAux
    variable = relperm_1_aux
    property = relperm
    phase = 1
  []
  [porosity]
    type = PorousFlowPropertyAux
    variable = porosity_aux
    property = porosity
  []
[]

[Functions]
  [pressure00]
    type = ParsedFunction
    value = '101325*10 + 9.81*(-y)*density0'
    vars = 'density0'
    vals = '${density0_phase0}'
  []
  [temperature00]
    type = ParsedFunction
    value = '273.15 + ttop + dTdz*(-y)'
    vars = 'ttop dTdz'
    vals = '5   0.1' # geothermal gradient dT/dz=100K/1km
  []
  [gas_injection]
    type = ParsedFunction
    value = 'vmax * min(t,1e6)/1e6'
    vars = 'vmax'
    vals = '-1e-10'
  []
[]

[ICs]
  [pgas_IC]
    type = FunctionIC
    variable = pgas
    function = pressure00    
  []
  [temperature_IC]
    type = FunctionIC
    variable = temperature
    function = temperature00
  []
  [ptop_auxvar_IC]
    type = FunctionIC
    variable = ptop_auxvar
    function = pressure00
  []
[]

[BCs]
  [constant_gas_injection]
    type = PorousFlowSink
    boundary = injection_area
    variable = zi
    flux_function = gas_injection
  []
  [mass_water_top]
    type = PorousFlowPiecewiseLinearSink
    boundary = top
    variable = pgas
    fluid_phase = 0
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 1e9'
    multipliers = '-1e9 1e9'
    use_mobility = true
    use_relperm = true
    flux_function = 100
    save_in = mass_water_flux_top
  []
  [heat_phase0_adv_top]
    type = PorousFlowPiecewiseLinearSink
    boundary = top
    variable = temperature
    fluid_phase = 0
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 1e9'
    multipliers = '-1e9 1e9'
    use_mobility = true
    use_relperm = true
    use_enthalpy = true
    flux_function = 100
    save_in = heat_flux_water_adv_top
  []
  [mass_methane_top]
    type = PorousFlowPiecewiseLinearSink
    boundary = top
    variable = pgas
    fluid_phase = 1
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 1e9'
    multipliers = '-1e9 1e9'
    use_mobility = true
    use_relperm = true
    flux_function = 100
    save_in = mass_methane_flux_top
  []
  [heat_phase1_adv_top]
    type = PorousFlowPiecewiseLinearSink
    boundary = top
    variable = temperature
    fluid_phase = 1
    PT_shift = ptop_auxvar
    pt_vals = '-1e9 1e9'
    multipliers = '-1e9 1e9'
    use_mobility = true
    use_relperm = true
    use_enthalpy = true
    flux_function = 100
    save_in = heat_flux_gas_adv_top
  []
  [heat_diff_top]
    type = PorousFlowPiecewiseLinearSink
    variable = temperature
    PT_shift = 278.15
    boundary = top
    pt_vals     = '-1000 0 1000'
    multipliers = '-1000 0 100'
    use_thermal_conductivity = true 
    flux_function = 0.1
    save_in = heat_flux_dif_top 
  []
  [Tbot]
    type = FunctionDirichletBC 
    variable = temperature
    boundary = bottom
    function = temperature00
  []
[]

[FluidProperties]
  [true_water]
    type = Water97FluidProperties
  []
  [true_methane]
    type = MethaneFluidProperties
  []
[]
[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [waterncg]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature = temperature
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosity
    fluid = true
    mechanical = false
    thermal = true
    porosity_zero = 0.1
    reference_temperature = 330
    reference_porepressure = 20E5
    thermal_expansion_coeff = 15E-6
    solid_bulk = 8E9
  []
  [effective_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []

  [permeability_aquifer]
    type = PorousFlowPermeabilityKozenyCarman
    block = aquifer
    poroperm_function = kozeny_carman_phi0
    phi0 = 0.1
    n = 2
    m = 2
    k0 = 1E-12
  []
  [permeability_caps]
    type = PorousFlowPermeabilityKozenyCarman
    block = caps
    poroperm_function = kozeny_carman_phi0
    phi0 = 0.1
    n = 2
    m = 2
    k0 = 1E-15
    k_anisotropy = '1 0 0  0 1 0  0 0 1'
  []

  [relperm_phase0] # liquid
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
    s_res = 0.1
    sum_s_res = 0.1
  []
  [relperm_phase1] # gas
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    s_res = 0.05
    sum_s_res = 0.15
  []
  [rock_internal_energy]
    type = PorousFlowMatrixInternalEnergy
    density = 2500
    specific_heat_capacity = 1100
  []
  [rock_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.5 0 0  0 2.5 0  0 0 2.5'
  []
[]

[Postprocessors]
[]

[Preconditioning]
  active = lu_mumps
  [basic]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres     asm      lu           NONZERO                   2'
  []
  [lu_mumps]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 315576000 # 10 years
  dtmax = 1e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1E2
    growth_factor = 1.2
    optimal_iterations = 10
  []
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-5
[]

[Outputs]
  exodus = true
[]
