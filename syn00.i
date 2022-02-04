# ../../archaea-opt -i syn00.i --mesh-only
# peacock -r syn00_in.e                                          # requires 'conda activate moose'
# In 2D GeneratedMeshGenerator automatically creates the boundaries bottom = 0, right = 1, top = 2, left = 3
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -200000
    xmax =   200000
    ymin = -150000
    ymax =   0
    elem_type = TRI6
    nx = 100
    ny = 75
  []
  #[mantle]
  #  type = SubdomainBoundingBoxGenerator
  #  block_id = 1
  #  bottom_left = '-200000 -150000 0'
  #  top_right = '200000 -40000 0'
  #  input = gen
  #[]
  [lcrust]
    type = SubdomainBoundingBoxGenerator
    input = 'gen'
    block_id = 1
    #block_name = lcrust
    bottom_left = '-200000 -40000 0'
    top_right = '200000 -20000 0'
  []
  [ucrust]
    type = SubdomainBoundingBoxGenerator
    input = 'lcrust'
    block_id = 2
    #block_name = ucrust
    bottom_left = '-200000 -20000 0'
    top_right = '200000 -10000 0'
  []
  [unsat]
    type = SubdomainBoundingBoxGenerator
    input = 'ucrust'
    block_id = 3
    #block_name = unsat
    bottom_left = '-200000 -10000 0'
    top_right = '200000 0 0'
  []
  [rename]
    type = RenameBlockGenerator
    input = 'unsat'
    old_block = '0 1 2 3'
    new_block = 'mantle lcrust ucrust unsat'
  []
  #[top]
  #  type = ParsedGeneratedSideset
  #  input = 'rename'
  #  combinatorial_geometry = 'y > -0.1'
  #  new_sideset_name = 'top'
  #[]
  #
  # SidesetsBetweenSubdomains
[]

[Variables]
  [dummy_var]
  []
[]
[Kernels]
  [dummy_diffusion]
    type = Diffusion
    variable = dummy_var
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
