################################################################################
export Shallow_Water_Physics!
################################################################################


function Shallow_Water_Physics!(;dyn_data::Dyn_Data, 
                                kappa_m::Float64 = 1.0,
                                kappa_h::Float64 = 1.0,
                                kappa_chi::Float64 = 1.0,
                                kappa_lambda::Float64 = 1.0)
    """
    Compute the grid tendencies of (u, v, h), and tracers due to 
    non-dynamical processes.
    """
    ########################################
    ### Main variable
    grid_δu = dyn_data.grid_δu
    grid_δv = dyn_data.grid_δv
    grid_δh = dyn_data.grid_δlnps
    # grid_δq = dyn_data.grid_δq
    
    ########################################
    ### Temporary variable
    grid_u = dyn_data.grid_u_c
    grid_v = dyn_data.grid_v_c
    grid_h = dyn_data.grid_ps_c
    # grid_q = dyn_data.grid_q_c
    
    ########################################
    ### Reset
    grid_δu .= 0.0
    grid_δv .= 0.0
    grid_δh .= 0.0
    # grid_δq .= 0.0
    
    ########################################
    ### Rayleigh damping
    # grid_δu .= kappa_m * grid_u
    # grid_δv .= kappa_m * grid_v
    # grid_δh .= kappa_t * grid_h
    # grid_δq .= kappa_c * grid_q
    
    ########################################
    ### Convective process
    # grid_precip .= kappa_chi * (grid_q - q_saturation)
    # grid_precip[grid_precip .< 0.0] .= 0.0
    # grid_δq .+= -grid_precip
    # grid_δh .+= kappa_lambda * grid_δq
    
end