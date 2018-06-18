waLBerla_link_files_to_builddir( "*.prm" )

waLBerla_add_executable ( NAME FluidizedBedMEM
                          FILES FluidizedBedMEM.cpp
                          DEPENDS blockforest boundary core domain_decomposition field lbm pe pe_coupling postprocessing timeloop vtk )