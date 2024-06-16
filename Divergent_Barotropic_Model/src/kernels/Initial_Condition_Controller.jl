################################################################################
export Background_Vorticity_Strip
export Isolated_Vorticity_Blob
export Apply_InverseLaplacian!
################################################################################

function Background_Vorticity_Strip(; mesh::Spectral_Spherical_Mesh,
                                    atmo_data::Atmo_Data,
                                    dyn_data::Dyn_Data,
                                    vor_amp::Float64=1.0,
                                    vor_lat::Float64=0.0,
                                    vor_width::Float64=1.0,
                                    DEG_TO_RAD::Float64=pi / 180)
    """
    Generate background vorticity and corresponding balance height.
    The actual vorticity amplitude is reduced due to implementation,
    so need some guess.
    """
    ########################################################
    ### Formatter
    nλ, nθ, nd = dyn_data.nλ, dyn_data.nθ, dyn_data.nd
    num_fourier, num_spherical = dyn_data.num_fourier, dyn_data.num_spherical

    ########################################################
    ### Main variable
    grid_u = zeros(Float64, nλ, nθ, nd)
    grid_v = zeros(Float64, nλ, nθ, nd)
    grid_vor = zeros(Float64, nλ, nθ, nd)
    grid_div = zeros(Float64, nλ, nθ, nd)
    grid_h = zeros(Float64, nλ, nθ, nd)

    spe_vor_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_div_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_h_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    ########################################################
    ### Temporaray variables
    θc, λc = mesh.θc, mesh.λc
    lats = reshape(θc' .* ones(nλ), (nλ, nθ, 1))

    grid_absvor = zeros(Float64, nλ, nθ, nd)
    grid_kin = zeros(Float64, nλ, nθ, nd)
    spe_kin = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δu = zeros(Float64, nλ, nθ, nd)
    grid_δv = zeros(Float64, nλ, nθ, nd)
    spe_δvor = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_δdiv = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    ########################################################
    ### Perturbation (first guess)
    grid_vor .= vor_amp
    grid_vor .*= exp.(-((lats .- vor_lat * DEG_TO_RAD) / (vor_width * DEG_TO_RAD)) .^ 2)
    grid_div .= 0.0

    ########################################################
    ### Convergence condition
    # Due to convergence condition required in 
    # Compute_Ucos_Vcos_From_Vor_Div!
    # the wind field in meridional boundary must vanish, 
    # meaning u = 0, v = 0 at θ = +- pi/2
    # Thus, given perturbations in vorticity form leads to
    # complex transformation and filtering.

    # ***using asymptotic heaviside step function
    Trans_Grid_To_Spherical!(mesh, grid_vor, spe_vor_c)
    Trans_Grid_To_Spherical!(mesh, grid_div, spe_div_c)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
    grid_u .*= (0.5 .+ 0.5 * tanh.(4.0 * (lats .+ (-10) * DEG_TO_RAD)))
    grid_u .*= (0.5 .+ 0.5 * tanh.(-4.0 * (lats .+ (-30) * DEG_TO_RAD)))

    ########################################################
    ### Output stage (for momentum field)
    Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)
    spe_div_c .= 0.0 # Non-divergent field
    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
    Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)

    ########################################################
    ### Balance height (thermal wind relation)
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    grid_δu .+= grid_absvor .* grid_v
    grid_δv .+= -grid_absvor .* grid_u
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    grid_kin .= 0.5 * (grid_u .^ 2 + grid_v .^ 2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
    spe_h_c .= spe_δdiv - spe_kin

    ########################################################
    ### Output stage (for mass field)
    spe_h_c[1, 1] = 0.0 # NO mean thickness
    Trans_Spherical_To_Grid!(mesh, spe_h_c, grid_h)

    return (grid_u, grid_v, grid_vor, grid_div, grid_h, spe_vor_c, spe_div_c, spe_h_c)
end

function Isolated_Vorticity_Blob(; mesh::Spectral_Spherical_Mesh,
                                 atmo_data::Atmo_Data,
                                 dyn_data::Dyn_Data,
                                 vor_amp::Float64=1.0,
                                 vor_lon::Float64=0.0,
                                 vor_lat::Float64=0.0,
                                 vor_width::Float64=1.0,
                                 DEG_TO_RAD::Float64=pi / 180)
    """
    Generate perturbation vorticity and corresponding balance height.
    The actual vorticity amplitude is reduced due to implementation,
    so need some guess.
    """
    ########################################################
    ### Formatter
    nλ, nθ, nd = dyn_data.nλ, dyn_data.nθ, dyn_data.nd
    num_fourier, num_spherical = dyn_data.num_fourier, dyn_data.num_spherical

    ########################################################
    ### Main variable
    grid_u = zeros(Float64, nλ, nθ, nd)
    grid_v = zeros(Float64, nλ, nθ, nd)
    grid_vor = zeros(Float64, nλ, nθ, nd)
    grid_div = zeros(Float64, nλ, nθ, nd)
    grid_h = zeros(Float64, nλ, nθ, nd)

    spe_vor_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_div_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_h_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    ########################################################
    ### Temporaray variables
    θc, λc = mesh.θc, mesh.λc
    lons = reshape(ones(nθ)' .* λc, (nλ, nθ, 1))
    lats = reshape(θc' .* ones(nλ), (nλ, nθ, 1))

    grid_absvor = zeros(Float64, nλ, nθ, nd)
    grid_kin = zeros(Float64, nλ, nθ, nd)
    spe_kin = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δu = zeros(Float64, nλ, nθ, nd)
    grid_δv = zeros(Float64, nλ, nθ, nd)
    spe_δvor = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_δdiv = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    grid_δh = zeros(Float64, nλ, nθ, nd)

    ########################################################
    ### Perturbation (first guess)
    grid_vor .= vor_amp
    grid_vor .*= exp.(-((lons .- vor_lon * DEG_TO_RAD) / (vor_width * DEG_TO_RAD)) .^ 2)
    grid_vor .*= exp.(-((lats .- vor_lat * DEG_TO_RAD) / (vor_width * DEG_TO_RAD)) .^ 2)
    grid_div .= 0.0

    ########################################################
    ### Convergence condition
    # Due to convergence condition required in 
    # Compute_Ucos_Vcos_From_Vor_Div!
    # the wind field in meridional boundary must vanish, 
    # meaning u = 0, v = 0 at θ = +- pi/2
    # Thus, given perturbations in vorticity form leads to
    # complex transformation and filtering.

    # ***using asymptotic heaviside step function
    Trans_Grid_To_Spherical!(mesh, grid_vor, spe_vor_c)
    Trans_Grid_To_Spherical!(mesh, grid_div, spe_div_c)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
    grid_u .*= (0.5 .+ 0.5 * tanh.(4.0 * (lats .+ (-10) * DEG_TO_RAD)))
    grid_u .*= (0.5 .+ 0.5 * tanh.(-4.0 * (lats .+ (-30) * DEG_TO_RAD)))
    grid_v .*= (0.5 .+ 0.5 * tanh.(4.0 * (lats .+ (-10) * DEG_TO_RAD)))
    grid_v .*= (0.5 .+ 0.5 * tanh.(-4.0 * (lats .+ (-30) * DEG_TO_RAD)))

    ########################################
    ### Output stage (for momentum field)

    # spe_vor, spe_div
    Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)
    spe_div_c .= 0.0 # Non-divergent field
    # grid u, grid v
    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)

    # grid_vor, grid_div
    Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)

    ########################################################
    ### Balance height (thermal wind relation)
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    grid_δu .+= grid_absvor .* grid_v
    grid_δv .+= -grid_absvor .* grid_u
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
    Apply_InverseLaplacian!(mesh, spe_δdiv)
    grid_kin .= 0.5 * (grid_u .^ 2 + grid_v .^ 2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
    spe_h_c .= spe_δdiv - spe_kin

    ########################################################
    ### Output stage (for mass field)
    spe_h_c[1, 1] = 0.0 # NO mean thickness
    Trans_Spherical_To_Grid!(mesh, spe_h_c, grid_h)

    return (grid_u, grid_v, grid_vor, grid_div, grid_h, spe_vor_c, spe_div_c, spe_h_c)
end

function Apply_InverseLaplacian!(mesh::Spectral_Spherical_Mesh,
                                 spherical_u::Array{ComplexF64,3})
    """
    Inverse operations of Apply_Laplacian! in Spectral_Spherical_Mesh.jl
    """
    eig = mesh.laplacian_eig
    spherical_u .= spherical_u ./ eig
    spherical_u[isnan.(spherical_u)] .= 0.0
    return nothing
end
