################################################################################
export Shallow_Water_Dynamics!
export Implicit_Correction!
################################################################################

function Shallow_Water_Dynamics!(; mesh::Spectral_Spherical_Mesh,
                                 atmo_data::Atmo_Data,
                                 dyn_data::Dyn_Data,
                                 integrator::Filtered_Leapfrog,
                                 h_0::Float64)
    """
    Execute dynamical process in simulation in each time step.

    Reference from GFDL's idealized dry GCM guide, you can find it online.
    https://www.gfdl.noaa.gov/idealized-spectral-models-quickstart/
    """
    """
    TODO:: 
    For horizontal advection, spe_h_c != spe_h_n
    """
    ########################################
    ### Grid space

    # zonal wind
    grid_u_n = dyn_data.grid_u_n
    grid_u = dyn_data.grid_u_c
    grid_δu = dyn_data.grid_δu

    # meridional wind
    grid_v_n = dyn_data.grid_v_n
    grid_v = dyn_data.grid_v_c
    grid_δv = dyn_data.grid_δv

    # vorticity
    grid_vor = dyn_data.grid_vor

    # divergence
    grid_div = dyn_data.grid_div

    # surface height
    grid_h_n = dyn_data.grid_ps_n
    grid_h = dyn_data.grid_ps_c
    grid_δh = dyn_data.grid_δlnps

    ########################################
    ### Spectral space

    # vorticity
    spe_vor_n = dyn_data.spe_vor_n
    spe_vor_c = dyn_data.spe_vor_c
    spe_vor_p = dyn_data.spe_vor_p
    spe_δvor = dyn_data.spe_δvor

    # divergence
    spe_div_n = dyn_data.spe_div_n
    spe_div_c = dyn_data.spe_div_c
    spe_div_p = dyn_data.spe_div_p
    spe_δdiv = dyn_data.spe_δdiv

    # surface height
    spe_h_n = dyn_data.spe_lnps_n
    spe_h_c = dyn_data.spe_lnps_c
    spe_h_p = dyn_data.spe_lnps_p
    spe_δh = dyn_data.spe_δlnps

    ########################################
    ### Calculation

    # compute (f + ζ)v and −(f + ζ)u on the grid at time t, 
    # and add to δu and δv respectively
    grid_absvor = dyn_data.grid_absvor
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    grid_δu .+= grid_absvor .* grid_v
    grid_δv .+= -grid_absvor .* grid_u

    # compute the divergence and curl of (δu, δv) to obtain δζ and δD
    # in the spectral domain
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)

    # compute the kinetic energy, transform to spectral domain,
    # add to Φ at t on the spectral, take the Laplacian, and 
    # add −∇^2E to the spectral divergence tendency δD;
    grid_kin, spe_kin = dyn_data.grid_d_full1, dyn_data.spe_d1
    spe_energy = dyn_data.spe_energy
    grid_kin .= 0.5 * (grid_u .^ 2 + grid_v .^ 2)
    Trans_Grid_To_Spherical!(mesh, grid_kin, spe_kin)
    spe_energy .= spe_kin + spe_h_c
    Apply_Laplacian!(mesh, spe_energy)
    spe_δdiv .-= spe_energy

    # add to δh on the grid the horizontal advection term −∇ · (vT), 
    # then convert δT to spectral domain
    Add_Horizontal_Advection!(mesh, spe_h_c, grid_u, grid_v, grid_δh)
    grid_δh .-= grid_h .* grid_div
    Trans_Grid_To_Spherical!(mesh, grid_δh, spe_δh)

    # correct (δT, δD, δ ln ps) tendencies to slow down gravity wave
    Implicit_Correction!(integrator, spe_div_c, spe_div_p,
                         spe_h_c, spe_h_p, h_0,
                         spe_δdiv, spe_δh)

    # add the harmonic damping to δζ, δD, δh in the spectral domain,
    # treating the damping implicitly
    Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
    Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
    Compute_Spectral_Damping!(integrator, spe_h_c, spe_h_p, spe_δh)

    # use leapfrog to generate spectral t + ∆t values of ζ, D, h, ln ps and
    # apply Robert filter to the t values
    Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
    Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
    Filtered_Leapfrog!(integrator, spe_δh, spe_h_p, spe_h_c, spe_h_n)

    # generate grid (u, v) from spectral (ζ, D) at t + ∆t
    UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)

    # compute grid (ζ, D, h), at t + ∆t
    Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
    return Trans_Spherical_To_Grid!(mesh, spe_h_n, grid_h_n)
end

function Implicit_Correction!(integrator::Filtered_Leapfrog,
                              spe_div_c::Array{ComplexF64,3},
                              spe_div_p::Array{ComplexF64,3},
                              spe_h_c::Array{ComplexF64,3},
                              spe_h_p::Array{ComplexF64,3},
                              h_0::Float64,
                              spe_δdiv::Array{ComplexF64,3},
                              spe_δh::Array{ComplexF64,3})
    """
    Slow down gravity waves.
    """
    implicit_coef = integrator.implicit_coef
    eigen = integrator.laplacian_eigen
    init_step = integrator.init_step

    if init_step
        Δt = integrator.Δt
    else
        Δt = 2.0 * integrator.Δt
    end

    if init_step
        # for the first time step 
        # spe_h_c = spe_h_p and spe_div_c = spe_div_p
    else
        spe_δdiv .+= eigen .* (spe_h_c - spe_h_p)
        spe_δh .+= h_0 * (spe_div_c - spe_div_p)
    end

    μ = implicit_coef * Δt
    μ2 = μ^2
    spe_δdiv .= (spe_δdiv .- μ * eigen .* spe_δh) ./ (1.0 .- μ2 * eigen * h_0)
    return spe_δh .-= μ * h_0 * spe_δdiv
end
