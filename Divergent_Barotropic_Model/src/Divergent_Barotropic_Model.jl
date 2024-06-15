################################################################################
# Custom module
module Divergent_Barotropic_Model
"""
TODO : Seperate SW and generalize method
"""

################################################################################
# Default package
import FFTW
import JLD2

################################################################################
# Shared
include("kernels/Atmo_Data.jl")
include("kernels/Gauss_And_Legendre.jl")
include("kernels/Spectral_Spherical_Mesh.jl")
include("kernels/Vert_Coordinate.jl")
include("kernels/Time_Integrator.jl")
include("kernels/Semi_Implicit.jl")
include("kernels/Press_And_Geopot.jl")

################################################################################
# Shallow Water System
include("kernels/Dyn_Data.jl")
include("kernels/Output_Manager.jl")
include("kernels/Shallow_Water_Physics.jl")
include("kernels/Shallow_Water_Dynamics.jl")
include("kernels/Initial_Condition_Controller.jl")

################################################################################
# ?

end
