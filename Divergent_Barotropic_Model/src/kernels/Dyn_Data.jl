################################################################################
# 
export Dyn_Data
# Time-forwarding
export Time_Advance!
################################################################################

mutable struct Dyn_Data
    """
    Define prognostic variables, and associated temporary variables.
    """
    #########################################################
    # Meta
    name::String

    #########################################################
    # Resolution
    num_fourier::Int64
    num_spherical::Int64
    nλ::Int64
    nθ::Int64
    nd::Int64
    num_tracers::Int64

    #########################################################
    ### Grid space

    # zonal velocity
    grid_u_n::Array{Float64,3}
    grid_u_c::Array{Float64,3}
    grid_u_p::Array{Float64,3}

    # meridional velocity
    grid_v_n::Array{Float64,3}
    grid_v_c::Array{Float64,3}
    grid_v_p::Array{Float64,3}

    # vorticity 
    grid_vor::Array{Float64,3}

    # divergence
    grid_div::Array{Float64,3}

    # surface height or surface pressure
    grid_ps_n::Array{Float64,3}
    grid_ps_c::Array{Float64,3}
    grid_ps_p::Array{Float64,3}

    # temperature
    grid_t_n::Array{Float64,3}
    grid_t_c::Array{Float64,3}
    grid_t_p::Array{Float64,3}

    # specific humidity
    grid_q_n::Array{Float64,3}
    grid_q_c::Array{Float64,3}
    grid_q_p::Array{Float64,3}

    # passive tracer
    grid_tracers_n::Array{Float64,4}
    grid_tracers_c::Array{Float64,4}
    grid_tracers_p::Array{Float64,4}

    #########################################################
    ### Spectral space

    # vorticity 
    spe_vor_n::Array{ComplexF64,3}
    spe_vor_c::Array{ComplexF64,3}
    spe_vor_p::Array{ComplexF64,3}

    # divergence
    spe_div_n::Array{ComplexF64,3}
    spe_div_c::Array{ComplexF64,3}
    spe_div_p::Array{ComplexF64,3}

    # surface height or surface pressure
    spe_lnps_n::Array{ComplexF64,3}
    spe_lnps_c::Array{ComplexF64,3}
    spe_lnps_p::Array{ComplexF64,3}

    # temperature
    spe_t_n::Array{ComplexF64,3}
    spe_t_c::Array{ComplexF64,3}
    spe_t_p::Array{ComplexF64,3}

    # specific humidity
    spe_q_n::Array{ComplexF64,3}
    spe_q_c::Array{ComplexF64,3}
    spe_q_p::Array{ComplexF64,3}

    # passive tracer
    spe_tracers_n::Array{ComplexF64,4}
    spe_tracers_c::Array{ComplexF64,4}
    spe_tracers_p::Array{ComplexF64,4}

    ############################################################
    ### Temporary variables

    # vor
    spe_δvor::Array{ComplexF64,3}
    grid_δvor::Array{Float64,3}

    # div
    spe_δdiv::Array{ComplexF64,3}
    grid_δdiv::Array{Float64,3}

    # w-e velocity tendency
    spec_δu::Array{ComplexF64,3}
    grid_δu::Array{Float64,3}

    # n-s velocity tendency
    spec_δv::Array{ComplexF64,3}
    grid_δv::Array{Float64,3}

    # pressure     
    spe_δlnps::Array{ComplexF64,3}
    grid_lnps::Array{Float64,3}
    grid_δlnps::Array{Float64,3}

    grid_δps::Array{Float64,3}
    grid_p_full::Array{Float64,3}   # pressure at full level
    grid_p_half::Array{Float64,3}   # pressure at half level
    grid_lnp_full::Array{Float64,3} # ln pressure at full level
    grid_lnp_half::Array{Float64,3} # ln pressure at half level
    grid_Δp::Array{Float64,3}       # pressure difference at each level

    # pressure gradient
    grid_dλ_ps::Array{Float64,3}
    grid_dθ_ps::Array{Float64,3}

    # temperature tendency
    spe_δt::Array{ComplexF64,3}
    grid_δt::Array{Float64,3}

    # absolute vor
    grid_absvor::Array{Float64,3}

    # vertical velocity
    grid_w_full::Array{Float64,3}
    grid_M_half::Array{Float64,3}

    # sum of potential energy and kinematic energy
    spe_energy::Array{ComplexF64,3}
    grid_energy_full::Array{Float64,3}

    # geopotential 
    grid_geopots::Array{Float64,3} # surface geopotential
    grid_geopot_full::Array{Float64,3}
    grid_geopot_half::Array{Float64,3}

    # z coordinate
    grid_z_full::Array{Float64,3}
    grid_z_half::Array{Float64,3}

    # equilibrium temperature in HS_Forcing
    grid_t_eq::Array{Float64,3}

    ## wrapper
    spe_d1::Array{ComplexF64,3}
    spe_d2::Array{ComplexF64,3}
    #
    grid_d_full1::Array{Float64,3}
    grid_d_full2::Array{Float64,3}
    #
    grid_d_half1::Array{Float64,3}
    grid_d_half2::Array{Float64,3}

    spe_zeros::Array{ComplexF64,3}
end

function Dyn_Data(; name::String,
                  num_fourier::Int64,
                  num_spherical::Int64,
                  nλ::Int64,
                  nθ::Int64,
                  nd::Int64,
                  num_tracers::Int64=0)
    """
    Define prognostic variables, and associated temporary variables.
    """
    #########################################################
    ### Grid space

    # zonal velocity
    grid_u_n = zeros(Float64, nλ, nθ, nd)
    grid_u_c = zeros(Float64, nλ, nθ, nd)
    grid_u_p = zeros(Float64, nλ, nθ, nd)

    # meridional velocity
    grid_v_n = zeros(Float64, nλ, nθ, nd)
    grid_v_c = zeros(Float64, nλ, nθ, nd)
    grid_v_p = zeros(Float64, nλ, nθ, nd)

    # vorticity 
    grid_vor = zeros(Float64, nλ, nθ, nd)

    # divergence
    grid_div = zeros(Float64, nλ, nθ, nd)

    # surface height or surface pressure
    grid_ps_n = zeros(Float64, nλ, nθ, 1)
    grid_ps_c = zeros(Float64, nλ, nθ, 1)
    grid_ps_p = zeros(Float64, nλ, nθ, 1)

    # temperature
    grid_t_n = zeros(Float64, nλ, nθ, nd)
    grid_t_c = zeros(Float64, nλ, nθ, nd)
    grid_t_p = zeros(Float64, nλ, nθ, nd)

    # specific humidity
    grid_q_n = zeros(Float64, nλ, nθ, nd)
    grid_q_c = zeros(Float64, nλ, nθ, nd)
    grid_q_p = zeros(Float64, nλ, nθ, nd)

    # grid tracer
    grid_tracers_n = zeros(Float64, nλ, nθ, nd, num_tracers)
    grid_tracers_c = zeros(Float64, nλ, nθ, nd, num_tracers)
    grid_tracers_p = zeros(Float64, nλ, nθ, nd, num_tracers)

    #########################################################
    ### Spectral space

    # vorticity 
    spe_vor_n = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_vor_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_vor_p = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    # divergence
    spe_div_n = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_div_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_div_p = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    # surface height or surface pressure
    spe_lnps_n = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, 1)
    spe_lnps_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, 1)
    spe_lnps_p = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, 1)

    # temperature
    spe_t_n = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_t_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_t_p = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    # specific humidity
    spe_q_n = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_q_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_q_p = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    # passive tracer
    spe_tracers_n = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd, num_tracers)
    spe_tracers_c = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd, num_tracers)
    spe_tracers_p = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd, num_tracers)

    ############################################################
    ### Temporary variables

    # vor
    spe_δvor = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δvor = zeros(Float64, nλ, nθ, nd)

    # div
    spe_δdiv = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δdiv = zeros(Float64, nλ, nθ, nd)

    # w-e velocity tendency
    spec_δu = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δu = zeros(Float64, nλ, nθ, nd)

    # n-s velocity tendency
    spec_δv = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δv = zeros(Float64, nλ, nθ, nd)

    # pressure     
    spe_δlnps = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, 1)
    grid_lnps = zeros(Float64, nλ, nθ, 1)
    grid_δlnps = zeros(Float64, nλ, nθ, 1)

    grid_δps = zeros(Float64, nλ, nθ, 1)

    grid_p_full = zeros(Float64, nλ, nθ, nd)
    grid_p_half = zeros(Float64, nλ, nθ, nd + 1)
    grid_lnp_full = zeros(Float64, nλ, nθ, nd)
    grid_lnp_half = zeros(Float64, nλ, nθ, nd + 1)
    grid_Δp = zeros(Float64, nλ, nθ, nd)

    # pressure gradient
    grid_dλ_ps = zeros(Float64, nλ, nθ, 1)
    grid_dθ_ps = zeros(Float64, nλ, nθ, 1)

    # temperature tendency
    spe_δt = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_δt = zeros(Float64, nλ, nθ, nd)

    # absolute vor
    grid_absvor = zeros(Float64, nλ, nθ, nd)

    # vertical velocity
    grid_w_full = zeros(Float64, nλ, nθ, nd)
    grid_M_half = zeros(Float64, nλ, nθ, nd + 1)

    # sum of potential energy and kinematic energy
    spe_energy = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    grid_energy_full = zeros(Float64, nλ, nθ, nd)

    # geopotential 
    grid_geopots = zeros(Float64, nλ, nθ, 1) # surface geopotential
    grid_geopot_full = zeros(Float64, nλ, nθ, nd)
    grid_geopot_half = zeros(Float64, nλ, nθ, nd + 1)

    # z coordinate
    grid_z_full = zeros(Float64, nλ, nθ, nd)
    grid_z_half = zeros(Float64, nλ, nθ, nd + 1)

    # equilibrium temperature
    grid_t_eq = zeros(Float64, nλ, nθ, nd)
    ## wrapper
    spe_d1 = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    spe_d2 = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    #
    grid_d_full1 = zeros(Float64, nλ, nθ, nd)
    grid_d_full2 = zeros(Float64, nλ, nθ, nd)
    #
    grid_d_half1 = zeros(Float64, nλ, nθ, nd + 1)
    grid_d_half2 = zeros(Float64, nλ, nθ, nd + 1)
    #
    spe_zeros = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)

    ########################################################################
    return Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd, num_tracers,
                    ########################################################################
                    grid_u_n, grid_u_c, grid_u_p,
                    grid_v_n, grid_v_c, grid_v_p,
                    grid_vor, grid_div,
                    grid_ps_n, grid_ps_c, grid_ps_p,
                    grid_t_n, grid_t_c, grid_t_p,
                    grid_q_n, grid_q_c, grid_q_p,
                    grid_tracers_n, grid_tracers_c, grid_tracers_p,
                    #########################################################################
                    spe_vor_n, spe_vor_c, spe_vor_p,
                    spe_div_n, spe_div_c, spe_div_p,
                    spe_lnps_n, spe_lnps_c, spe_lnps_p,
                    spe_t_n, spe_t_c, spe_t_p,
                    spe_q_n, spe_q_c, spe_q_p,
                    spe_tracers_n, spe_tracers_c, spe_tracers_p,
                    #########################################################################
                    spe_δvor, grid_δvor,
                    spe_δdiv, grid_δdiv,
                    spec_δu, grid_δu,
                    spec_δv, grid_δv,
                    spe_δlnps, grid_lnps, grid_δlnps, grid_δps,
                    grid_p_full, grid_p_half, grid_lnp_full, grid_lnp_half, grid_Δp,
                    grid_dλ_ps, grid_dθ_ps,
                    spe_δt, grid_δt,
                    grid_absvor,
                    grid_w_full, grid_M_half,
                    spe_energy, grid_energy_full,
                    grid_geopots, grid_geopot_full, grid_geopot_half,
                    grid_z_full, grid_z_half, grid_t_eq,
                    #########################################################################
                    spe_d1, spe_d2,
                    grid_d_full1, grid_d_full2,
                    grid_d_half1, grid_d_half2,
                    spe_zeros)
    ########################################################################
end

function Time_Advance!(; dyn_data::Dyn_Data)
    """
    Update model states, following leapfrog scheme.
    """

    # forward grid variables
    dyn_data.grid_u_p .= dyn_data.grid_u_c
    dyn_data.grid_u_c .= dyn_data.grid_u_n

    dyn_data.grid_v_p .= dyn_data.grid_v_c
    dyn_data.grid_v_c .= dyn_data.grid_v_n

    dyn_data.grid_ps_p .= dyn_data.grid_ps_c
    dyn_data.grid_ps_c .= dyn_data.grid_ps_n

    dyn_data.grid_t_p .= dyn_data.grid_t_c
    dyn_data.grid_t_c .= dyn_data.grid_t_n

    if dyn_data.num_tracers > 0
        dyn_data.grid_tracers_p .= dyn_data.grid_tracers_c
        dyn_data.grid_tracers_c .= dyn_data.grid_tracers_n
    end

    # forward spectral variables
    dyn_data.spe_vor_p .= dyn_data.spe_vor_c
    dyn_data.spe_vor_c .= dyn_data.spe_vor_n

    dyn_data.spe_div_p .= dyn_data.spe_div_c
    dyn_data.spe_div_c .= dyn_data.spe_div_n

    dyn_data.spe_lnps_p .= dyn_data.spe_lnps_c
    dyn_data.spe_lnps_c .= dyn_data.spe_lnps_n

    dyn_data.spe_t_p .= dyn_data.spe_t_c
    dyn_data.spe_t_c .= dyn_data.spe_t_n

    if dyn_data.num_tracers > 0
        dyn_data.spe_tracers_p .= dyn_data.spe_tracers_c
        dyn_data.spe_tracers_c .= dyn_data.spe_tracers_n
    end
end
