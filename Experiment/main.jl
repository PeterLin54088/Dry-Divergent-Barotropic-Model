# Version
versioninfo()

# Environment initialization
using Pkg: Pkg
cd()
Pkg.activate("Dynamical_Core/Divergent_Barotropic_Model")
cd("Dynamical_Core")

# Model initialization
include(joinpath(pwd(), "Experiment/__main.jl"))

# Run
Shallow_Water_Main()
