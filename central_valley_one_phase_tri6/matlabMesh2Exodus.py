# main for exportMeshPy()
# run as:
# python3 matlabMesh2Exodus.py mesh_000000.mat 1

# the following entities are present in an ExodusII file:
# nodes           :: a location in 3D space
# element blocks  :: an element block is a collection of elements of the same type.
# elements        :: each element block in the model can contain any number of elements
# side sets       :: a side set is a collection of element faces: internally these are stored by listing a pair of the form (element_index, face_index) 
# node sets       :: a node set is a collection of nodes
# timesteps       :: a file may have any number of timesteps, including zero. Field values are defined for each timestep.
# global variables:: a global variable has a single alue for each timestep. 
# node fields     :: a node field lists the values of the field at each node for each timestep.
# element fields  :: an element field lists the value of the field at each element within a given block for each timestep.
#                    Element fields may not be defined on every element block, but if the are defined, a value is defined for each element within that block.
# side set fields :: similar to element fields, but for side sets.
# node set fields :: similar to element fields, but for node sets.
# information records :: list of strings wihtin the file.
# QA information  :: a list of quality assurance may be defined in the file. These hold the names and version numbers of eah program used to create or modify the file.

# An EXODUS file is a netCDF file, an application program can access data via the EXODUS API or via netCDF API function calls directly
# The EXODUS-II models in described in in exodusII-new.pdf and the exomerge python package used here to fill the input 
# EXODUS file into PorousFlow is described in exomerge.pdf
#
# The EXODUS model is described by data which are static (do not change through time). This data includes nodal coordinates,
# element connectivity (node lists for each element), element attributes, and node sets and side sets (used to apply
# loading conditions and boundary conditions).
# The number of nodes, element blocks, side sets and other entities are fixed and once the file is created these cannot change.
# The results are optional and include five types of variables --nodal, element, nodeset, sideset, and global--, each of which is
# stored through time.
#
# the input assumes a [0-based] connectivity matrix, either "TRI3" or "TRI6" as in libmesh (i.e. Hughes' notation)

# import importlib
# import matplotlib.tri as tri
import numpy as np
import os
import scipy.io as sio
import sys
import meshio

# HOME = os.environ['HOME']
# WORK = os.environ['WORK']
SEACAS_DIR = os.environ['SEACAS_DIR']
sys.path.append(os.path.join(SEACAS_DIR,"lib"))         # {exodus.py, exomerge.py}
#import exodus
import exomerge

# raise ValueError("enough is enough")
def main():
    fnamei = sys.argv[1]                                       #fnamei = "exportMeshPy_mesh_1380_i.mat"; verify=True
    verify = bool(int(sys.argv[2]))
    # import os  
    #os.chdir('/Users/javier/scratch/rift2ridge2/SIU3P/data/siu3p.nchem_hydr_melt_serp_spro.tests.mesh.central_valley_one_phase_tri6.001.001')  
    #fnamei = "mesh_000000.mat"; verify=True
    env = sio.loadmat(fnamei, verify_compressed_data_integrity=verify, squeeze_me=True)

    GCOORD = env['MESH']['GCOORD'].item()     # [ndim=2,nnod]
    EL2NOD = env['MESH']['EL2NOD'].item() - 1 # [nnodel,nel] {Tri3,Tri6}

    nnodel, nel = EL2NOD.shape
    eltype = "triangle" if nnodel == 3 else "triangle6"

    mesh = meshio.Mesh(                           # meshio mesh object
            GCOORD.T,                             # nodes                           -> mesh.points     :: numpy.ndarray [nnod,ndim]
            [(eltype,EL2NOD.T)]                   # cells [list ot 'tuple' objects] -> mesh.cell       :: list of [nblk] CellBlock objects
        )
    meshio_meshfile = fnamei.replace(".mat","_meshio.exo")
    mesh.write(meshio_meshfile)                            # meshio: export raw mesh. This raw mesh has one timestep set to 0.0

    #model = exomerge.import_model('meshio_out.exo', timesteps='none') # only if we want to disregard the meshio default created timestep=0.0
    #model.create_timestep(0.0) # required to start adding variables   # and create it manually
    model = exomerge.import_model(meshio_meshfile)                    # == model = ExodusModel(); model.import_model(meshio_meshfile) 
    os.remove(meshio_meshfile)

    # retrieve some information
    #model.timestep_exists(0.0)                                       # True
    #model.summarize()                                                # np.array(model.get_connectivity())[0:6] == EL2NOD[:,0]
    #model._detailed_summary()
    #model.get_closest_node_distance()
    #model._get_element_edge_indices('tri6')                          #
    PhaseID = env['MESH']['PhaseID'].item()                           # (nel,) dtype('uint32') or dtype('uint8) - depending on mesh size
    
    eltype = "TRI3" if nnodel == 3 else "TRI6"

    isr2r = True
    phids = np.unique(PhaseID)
    for pid in phids:
        model.create_element_block(pid, [eltype, np.sum(PhaseID == pid),nnodel,0],
            EL2NOD[:,PhaseID == pid].T.flatten().tolist()) # see exomerge.pdf::p.31
        if isr2r:
            # phid:        1        2        3        4          5        6        7
            elblknames = ["mantle","lcrust","ucrust","sediments","ocrust","nopore","unsat"]
            model.rename_element_block(pid, elblknames[pid-1])

    model.delete_element_block(0)

    # commonly there will be the following four nodesets:
    # ["bottom","right","top","left"]
    # plus two nodesets distinguising the top ubmarine and subaerial parts
    # ["top_subaerial","top_submarine"]
    BCnodeset_names = env['MESH']['BCnodesets'].item().dtype.names # tuple
    BCnodesets = {}
    for name in BCnodeset_names: # commonly ['bottom','right','top','letf'] + top subsets
        BCnodesets[name] = env['MESH']['BCnodesets'].item()[name].item() - 1 # dtype('uint32') map to 0-based python indexing

    for id, key in enumerate(BCnodesets.keys()):
        model.create_node_set(id, BCnodesets[key].tolist())          # unnamed 0-based indexing
        model.rename_node_set(id,key)                                # in the seacas program explore. They can now be listed as EXPLORE>LIST NSETS 

    model.create_node_field('temperature', 298.15)                   # example: initialize

    exomerge_meshfile = fnamei.replace(".mat","_exomerge.e")         # archaea_??????_i.exomerge.exo
    model.export_model(exomerge_meshfile)                            # write out the current model to an ExodusII file

if __name__ == "__main__":
    main()
