################################################################################
# 
export Output_Manager
# 
export Update_Output_Init!, Update_Output!, Finalize_Output!
# 
export Generate_Output!
# 
export Display_Initial_Background!
export Display_Initial_Perturbation!
export Display_Current_State!
################################################################################


mutable struct Output_Manager
    """
    Define model (tunable) parameters and outputs.
    TODO:: Add relaxation parameter.
    """
    # Meta
    creation_time::String
    datapath::String
    
    # Convert
    hour_to_sec::Int64
    day_to_sec::Int64
    
    # Resolution
    num_fourier::Int64
    num_spherical::Int64
    nλ::Int64
    nθ::Int64
    nd::Int64
    
    # Integrator
        # time
    start_time::Int64
    end_time::Int64
    current_time::Int64
    n_hour::Int64
    n_day::Int64
    Δt::Int64
        # hyper diffusion
    damping_order::Float64
    damping_coef::Float64
    robert_coef::Float64
    
    # Coordinate
    λc::Array{Float64, 1}
    θc::Array{Float64, 1}
    σc::Array{Float64, 1}

    # Dataset (nλ × nθ × nd)
    u_init::Array{Float64, 3}
    v_init::Array{Float64, 3}
    h_init::Array{Float64, 3}
    q_init::Array{Float64, 3}
    vor_init::Array{Float64, 3}
    div_init::Array{Float64, 3}
    
    # Dataset (nλ × nθ × nd × n_hour)
    u_hourly_mean::Array{Float64, 4}
    v_hourly_mean::Array{Float64, 4}
    h_hourly_mean::Array{Float64, 4}
    q_hourly_mean::Array{Float64, 4}
    vor_hourly_mean::Array{Float64, 4}
    div_hourly_mean::Array{Float64, 4}
    hourly_mean_counter::Array{Int64, 1} # data counter
    
    # Dataset (nλ × nθ × nd × n_day)
    u_daily_mean::Array{Float64, 4}
    v_daily_mean::Array{Float64, 4}
    h_daily_mean::Array{Float64, 4}
    q_daily_mean::Array{Float64, 4}
    vor_daily_mean::Array{Float64, 4}
    div_daily_mean::Array{Float64, 4}
    daily_mean_counter::Array{Int64, 1} # data counter
    
end


function Output_Manager(;mesh::Spectral_Spherical_Mesh,
                        integrator::Filtered_Leapfrog,
                        vert_coord::Vert_Coordinate, 
                        start_time::Int64, 
                        end_time::Int64,
                        creation_time::String,
                        datapath::String,
                        hour_to_sec = 3600,
                        day_to_sec = 86400)
    """
    Define model (tunable) parameters and outputs.
    """
    # Resolution
    num_fourier = mesh.num_fourier
    num_spherical = mesh.num_spherical
    nλ = mesh.nλ
    nθ = mesh.nθ
    nd = mesh.nd
    
    # Integrator
        # time
    current_time = start_time
    n_hour = div(end_time - start_time, hour_to_sec)
    n_day = div(end_time - start_time, day_to_sec)
    Δt = integrator.Δt
        # hyper diffusion
    damping_order = integrator.damping_order
    damping_coef = integrator.damping_coef
    robert_coef = integrator.robert_coef
     
    # Coordinate
    λc = mesh.λc
    θc = mesh.θc
    bk = vert_coord.bk
    σc = (bk[2:nd+1] + bk[1:nd])/2.0

    # Dataset (nλ × nθ × nd)
    u_init = zeros(Float64, nλ, nθ, nd)
    v_init = zeros(Float64, nλ, nθ, nd)
    h_init = zeros(Float64, nλ, nθ, 1)
    q_init = zeros(Float64, nλ, nθ, nd)
    vor_init = zeros(Float64, nλ, nθ, nd)
    div_init = zeros(Float64, nλ, nθ, nd)
    
    # Dataset (nλ × nθ × nd × n_hour)
    u_hourly_mean = zeros(Float64, nλ, nθ, nd, n_hour)
    v_hourly_mean = zeros(Float64, nλ, nθ, nd, n_hour)
    h_hourly_mean = zeros(Float64, nλ, nθ, 1, n_hour)
    q_hourly_mean = zeros(Float64, nλ, nθ, nd, n_hour)
    vor_hourly_mean = zeros(Float64, nλ, nθ, nd, n_hour)
    div_hourly_mean = zeros(Float64, nλ, nθ, nd, n_hour)
    hourly_mean_counter = zeros(Int64, n_hour) # data counter
    
    # Dataset (nλ × nθ × nd × n_day)
    u_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    v_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    h_daily_mean = zeros(Float64, nλ, nθ, 1, n_day)
    q_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    vor_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    div_daily_mean = zeros(Float64, nλ, nθ, nd, n_day)
    daily_mean_counter = zeros(Int64, n_day) # data counter
    
    ########################################################################
    output_manager = Output_Manager(creation_time, datapath,
    ########################################################################
                                    hour_to_sec, day_to_sec,
    ########################################################################
                                    num_fourier, num_spherical, nλ, nθ, nd,
    ########################################################################
                                    start_time, end_time, current_time, 
                                    n_hour, n_day, Δt, 
                                    damping_order, damping_coef, robert_coef,
    ########################################################################
                                    λc, θc, σc,
    ########################################################################
                                    u_init, v_init, h_init, q_init,
                                    vor_init, div_init,
    ########################################################################
                                    u_hourly_mean, v_hourly_mean, 
                                    h_hourly_mean, q_hourly_mean, 
                                    vor_hourly_mean, div_hourly_mean, 
                                    hourly_mean_counter,
    ########################################################################
                                    u_daily_mean, v_daily_mean, 
                                    h_daily_mean, q_daily_mean, 
                                    vor_daily_mean, div_daily_mean, 
                                    daily_mean_counter)
    ########################################################################
end


function Update_Output_Init!(;output_manager::Output_Manager, 
                             dyn_data::Dyn_Data, 
                             current_time::Int64)
    """
    Update outputs before any physical or dynamical process. 
    """
    @assert(current_time == output_manager.current_time)

    u_init = output_manager.u_init
    v_init = output_manager.v_init
    h_init = output_manager.h_init
    q_init = output_manager.q_init
    vor_init = output_manager.vor_init
    div_init = output_manager.div_init
    
    u_hourly_mean = output_manager.u_hourly_mean
    v_hourly_mean = output_manager.v_hourly_mean
    h_hourly_mean = output_manager.h_hourly_mean
    q_hourly_mean = output_manager.q_hourly_mean
    vor_hourly_mean = output_manager.vor_hourly_mean
    div_hourly_mean = output_manager.div_hourly_mean
    hourly_mean_counter = output_manager.hourly_mean_counter
    
    u_daily_mean = output_manager.u_daily_mean
    v_daily_mean = output_manager.v_daily_mean
    h_daily_mean = output_manager.h_daily_mean
    q_daily_mean = output_manager.q_daily_mean
    vor_daily_mean = output_manager.vor_daily_mean
    div_daily_mean = output_manager.div_daily_mean
    daily_mean_counter = output_manager.daily_mean_counter

    t_index = 1

    u_init[:,:,:] .+= dyn_data.grid_u_c
    v_init[:,:,:] .+= dyn_data.grid_v_c
    h_init[:,:,1] .+= dyn_data.grid_ps_c[:,:,1]
    q_init[:,:,:] .+= dyn_data.grid_q_c
    vor_init[:,:,:] .+= dyn_data.grid_vor
    div_init[:,:,:] .+= dyn_data.grid_div
    
    u_hourly_mean[:,:,:,t_index] .+= dyn_data.grid_u_c
    v_hourly_mean[:,:,:,t_index] .+= dyn_data.grid_v_c
    h_hourly_mean[:,:,1,t_index] .+= dyn_data.grid_ps_c[:,:,1]
    q_hourly_mean[:,:,:,t_index] .+= dyn_data.grid_q_c
    vor_hourly_mean[:,:,:,t_index] .+= dyn_data.grid_vor
    div_hourly_mean[:,:,:,t_index] .+= dyn_data.grid_div
    hourly_mean_counter[t_index] += 1
    
    u_daily_mean[:,:,:,t_index] .+= dyn_data.grid_u_c
    v_daily_mean[:,:,:,t_index] .+= dyn_data.grid_v_c
    h_daily_mean[:,:,1,t_index] .+= dyn_data.grid_ps_c[:,:,1]
    q_daily_mean[:,:,:,t_index] .+= dyn_data.grid_q_c
    vor_daily_mean[:,:,:,t_index] .+= dyn_data.grid_vor
    div_daily_mean[:,:,:,t_index] .+= dyn_data.grid_div
    daily_mean_counter[t_index] += 1
end


function Update_Output!(;output_manager::Output_Manager, 
                        dyn_data::Dyn_Data, 
                        current_time::Int64)
    """
    Update outputs between each time interval. 
    """
    @assert(current_time > output_manager.current_time)
    hour_to_sec = output_manager.hour_to_sec
    day_to_sec = output_manager.day_to_sec
    start_time = output_manager.start_time
    
    u_hourly_mean = output_manager.u_hourly_mean
    v_hourly_mean = output_manager.v_hourly_mean
    h_hourly_mean = output_manager.h_hourly_mean
    q_hourly_mean = output_manager.q_hourly_mean
    vor_hourly_mean = output_manager.vor_hourly_mean
    div_hourly_mean = output_manager.div_hourly_mean
    hourly_mean_counter = output_manager.hourly_mean_counter
    
    u_daily_mean = output_manager.u_daily_mean
    v_daily_mean = output_manager.v_daily_mean
    h_daily_mean = output_manager.h_daily_mean
    q_daily_mean = output_manager.q_daily_mean
    vor_daily_mean = output_manager.vor_daily_mean
    div_daily_mean = output_manager.div_daily_mean
    daily_mean_counter = output_manager.daily_mean_counter

    t_index_hourly = 1 + div(current_time - start_time, hour_to_sec)
    t_index_daily = 1 + div(current_time - start_time, day_to_sec)
    
    u_hourly_mean[:,:,:,t_index_hourly] .+= dyn_data.grid_u_c
    v_hourly_mean[:,:,:,t_index_hourly] .+= dyn_data.grid_v_c
    h_hourly_mean[:,:,1,t_index_hourly] .+= dyn_data.grid_ps_c[:,:,1]
    q_hourly_mean[:,:,:,t_index_hourly] .+= dyn_data.grid_q_c
    vor_hourly_mean[:,:,:,t_index_hourly] .+= dyn_data.grid_vor
    div_hourly_mean[:,:,:,t_index_hourly] .+= dyn_data.grid_div
    hourly_mean_counter[t_index_hourly] += 1
    
    u_daily_mean[:,:,:,t_index_daily] .+= dyn_data.grid_u_c
    v_daily_mean[:,:,:,t_index_daily] .+= dyn_data.grid_v_c
    h_daily_mean[:,:,1,t_index_daily] .+= dyn_data.grid_ps_c[:,:,1]
    q_daily_mean[:,:,:,t_index_daily] .+= dyn_data.grid_q_c
    vor_daily_mean[:,:,:,t_index_daily] .+= dyn_data.grid_vor
    div_daily_mean[:,:,:,t_index_daily] .+= dyn_data.grid_div
    daily_mean_counter[t_index_daily] += 1
    
    output_manager.current_time = current_time
end


function Finalize_Output!(;output_manager::Output_Manager)
    """
    Calculate "daily" mean. 
    """
    u_hourly_mean = output_manager.u_hourly_mean
    v_hourly_mean = output_manager.v_hourly_mean
    h_hourly_mean = output_manager.h_hourly_mean
    q_hourly_mean = output_manager.q_hourly_mean
    vor_hourly_mean = output_manager.vor_hourly_mean
    div_hourly_mean = output_manager.div_hourly_mean
    hourly_mean_counter = output_manager.hourly_mean_counter
    
    for time = 1:output_manager.n_hour
        u_hourly_mean[:,:,:,time] ./= hourly_mean_counter[time]
        v_hourly_mean[:,:,:,time] ./= hourly_mean_counter[time]
        h_hourly_mean[:,:,1,time] ./= hourly_mean_counter[time]
        q_hourly_mean[:,:,:,time] ./= hourly_mean_counter[time]
        vor_hourly_mean[:,:,:,time] ./= hourly_mean_counter[time]
        div_hourly_mean[:,:,:,time] ./= hourly_mean_counter[time]
        hourly_mean_counter[time] = 1
    end
    
    u_daily_mean = output_manager.u_daily_mean
    v_daily_mean = output_manager.v_daily_mean
    h_daily_mean = output_manager.h_daily_mean
    q_daily_mean = output_manager.q_daily_mean
    vor_daily_mean = output_manager.vor_daily_mean
    div_daily_mean = output_manager.div_daily_mean
    daily_mean_counter = output_manager.daily_mean_counter
    
    for time = 1:output_manager.n_day
        u_daily_mean[:,:,:,time] ./= daily_mean_counter[time]
        v_daily_mean[:,:,:,time] ./= daily_mean_counter[time]
        h_daily_mean[:,:,1,time] ./= daily_mean_counter[time]
        q_daily_mean[:,:,:,time] ./= daily_mean_counter[time]
        vor_daily_mean[:,:,:,time] ./= daily_mean_counter[time]
        div_daily_mean[:,:,:,time] ./= daily_mean_counter[time]
        daily_mean_counter[time] = 1
    end
end


function Generate_Output!(;output_manager::Output_Manager)
    """
    Write outputs into HDF5 file.
    """
    JLD2.jldopen(output_manager.datapath, "a+") do file
        
        # Meta
        file["creation_time"] = output_manager.creation_time
        
        # Resolution
            # Time
        file["num_daily_data"] = output_manager.n_day
            # Horizontal
        file["num_longitude"] = output_manager.nλ
        file["num_latitude"] = output_manager.nθ
        file["num_zonal_spectral_modes"] = output_manager.num_fourier
        file["num_meridional_spectral_modes"] = output_manager.num_spherical
            # Vertical
        file["num_vertical"] = output_manager.nd
        
        # Coordinate
        file["longitude"] = output_manager.λc
        file["latitude"] = output_manager.θc
        file["hybrid_levels"] = output_manager.σc
        
        # Integrator
        file["start_time"] = output_manager.start_time
        file["end_time"] = output_manager.end_time
        file["Δt"] = output_manager.Δt
        file["damping_order"] = output_manager.damping_order
        file["damping_coef"] = output_manager.damping_order
        file["robert_coef"] = output_manager.damping_order
        
        # Dataset
        file["u_hourly_mean"] = output_manager.u_hourly_mean
        file["v_hourly_mean"] = output_manager.v_hourly_mean
        file["h_hourly_mean"] = output_manager.h_hourly_mean
        file["q_hourly_mean"] = output_manager.q_hourly_mean
        file["vor_hourly_mean"] = output_manager.vor_hourly_mean
        file["div_hourly_mean"] = output_manager.div_hourly_mean
        file["u_daily_mean"] = output_manager.u_daily_mean
        file["v_daily_mean"] = output_manager.v_daily_mean
        file["h_daily_mean"] = output_manager.h_daily_mean
        file["q_daily_mean"] = output_manager.q_daily_mean
        file["vor_daily_mean"] = output_manager.vor_daily_mean
        file["div_daily_mean"] = output_manager.div_daily_mean
    end
end


function Display_Initial_Background!(;logpath::String,
                                     grid_u::Array{Float64,3},
                                     grid_v::Array{Float64,3},
                                     grid_h::Array{Float64,3},
                                     grid_vor::Array{Float64,3},
                                     grid_div::Array{Float64,3})
    
    # Display on terminal
    println(repeat("###", 30))
    println("Initial background      zonal wind: ",
            (round(minimum(grid_u); digits = 4), 
            round(maximum(grid_u); digits = 4)),
            " (m/s)")
    println("Initial background meridional wind: ",
            (round(minimum(grid_v); digits = 4), 
            round(maximum(grid_v); digits = 4)),
            " (m/s)")
    println("Initial background   geopot height: ", 
            (round(minimum(grid_h); digits = 4), 
            round(maximum(grid_h); digits = 4)),
            " (m * m*s-2)")
    println("Initial background       vorticity: ", 
            (round(minimum(grid_vor); digits = 9), 
            round(maximum(grid_vor); digits = 9)),
            " (1/s)")
    println("Initial background      divergence: ", 
            (round(minimum(grid_div); digits = 9), 
            round(maximum(grid_div); digits = 9)),
            " (1/s)")
    println(repeat("###", 30))
    
    # Display on file
    open(logpath, "a+") do file
        #
        write(file, repeat("###", 30), "\n")
        #
        write(file, "Initial background      zonal wind: ")
        write(file, "( ", string(round(minimum(grid_u); digits = 4)), ", ")
        write(file, string(round(maximum(grid_u); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial background meridional wind: ")
        write(file, "( ", string(round(minimum(grid_v); digits = 4)), ", ")
        write(file, string(round(maximum(grid_v); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial background   geopot height: ")
        write(file, "( ", string(round(minimum(grid_h); digits = 4)), ", ")
        write(file, string(round(maximum(grid_h); digits = 4)), " )")
        write(file, " (m * m*s-2)\n")
        #
        write(file, "Initial background       vorticity: ")
        write(file, "( ", string(round(minimum(grid_vor); digits = 9)), ", ")
        write(file, string(round(maximum(grid_vor); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, "Initial background      divergence: ")
        write(file, "( ", string(round(minimum(grid_div); digits = 9)), ", ")
        write(file, string(round(maximum(grid_div); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, repeat("###", 30), "\n")
    end
end


function Display_Initial_Perturbation!(;logpath::String,
                                       grid_u::Array{Float64,3},
                                       grid_v::Array{Float64,3},
                                       grid_h::Array{Float64,3},
                                       grid_vor::Array{Float64,3},
                                       grid_div::Array{Float64,3})
    
    # Display on terminal
    println(repeat("###", 30))
    println("Initial perturbation      zonal wind: ",
            (round(minimum(grid_u); digits = 4), 
            round(maximum(grid_u); digits = 4)),
            " (m/s)")
    println("Initial perturbation meridional wind: ",
            (round(minimum(grid_v); digits = 4), 
            round(maximum(grid_v); digits = 4)),
            " (m/s)")
    println("Initial perturbation   geopot height: ", 
            (round(minimum(grid_h); digits = 4), 
            round(maximum(grid_h); digits = 4)),
            " (m * m*s-2)")
    println("Initial perturbation       vorticity: ", 
            (round(minimum(grid_vor); digits = 9), 
            round(maximum(grid_vor); digits = 9)),
            " (1/s)")
    println("Initial perturbation      divergence: ", 
            (round(minimum(grid_div); digits = 9), 
            round(maximum(grid_div); digits = 9)),
            " (1/s)")
    println(repeat("###", 30))
    
    # Display on file
    open(logpath, "a+") do file
        #
        write(file, repeat("###", 30), "\n")
        #
        write(file, "Initial perturbation      zonal wind: ")
        write(file, "( ", string(round(minimum(grid_u); digits = 4)), ", ")
        write(file, string(round(maximum(grid_u); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial perturbation meridional wind: ")
        write(file, "( ", string(round(minimum(grid_v); digits = 4)), ", ")
        write(file, string(round(maximum(grid_v); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial perturbation   geopot height: ")
        write(file, "( ", string(round(minimum(grid_h); digits = 4)), ", ")
        write(file, string(round(maximum(grid_h); digits = 4)), " )")
        write(file, " (m * m*s-2)\n")
        #
        write(file, "Initial perturbation       vorticity: ")
        write(file, "( ", string(round(minimum(grid_vor); digits = 9)), ", ")
        write(file, string(round(maximum(grid_vor); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, "Initial perturbation      divergence: ")
        write(file, "( ", string(round(minimum(grid_div); digits = 9)), ", ")
        write(file, string(round(maximum(grid_div); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, repeat("###", 30), "\n") 
    end
end


function Display_Current_State!(;logpath::String,
                                tick::Int64,
                                grid_u::Array{Float64,3},
                                grid_v::Array{Float64,3},
                                grid_h::Array{Float64,3},
                                grid_vor::Array{Float64,3},
                                grid_div::Array{Float64,3})
    # Display on terminal
    println(repeat("***", 30))
    println("Hour: ", string(tick))
    println("     zonal wind: ",
            (round(minimum(grid_u); digits = 4), 
            round(maximum(grid_u); digits = 4)),
            " (m/s)")
    println("meridional wind: ",
            (round(minimum(grid_v); digits = 4), 
            round(maximum(grid_v); digits = 4)),
            " (m/s)")
    println("  geopot height: ", 
            (round(minimum(grid_h); digits = 4), 
            round(maximum(grid_h); digits = 4)),
            " (m * m*s-2)")
    println("      vorticity: ", 
            (round(minimum(grid_vor); digits = 9), 
            round(maximum(grid_vor); digits = 9)),
            " (1/s)")
    println("     divergence: ", 
            (round(minimum(grid_div); digits = 9), 
            round(maximum(grid_div); digits = 9)),
            " (1/s)")
    println(repeat("***", 30))
    
    # Display on file
    open(logpath, "a+") do file
        #
        write(file, repeat("***", 30), "\n")
        #
        write(file, "Hour: ", string(tick), "\n")
        #
        write(file, "     zonal wind: ")
        write(file, "( ", string(round(minimum(grid_u); digits = 4)), ", ")
        write(file, string(round(maximum(grid_u); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "meridional wind: ")
        write(file, "( ", string(round(minimum(grid_v); digits = 4)), ", ")
        write(file, string(round(maximum(grid_v); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "  geopot height: ")
        write(file, "( ", string(round(minimum(grid_h); digits = 4)), ", ")
        write(file, string(round(maximum(grid_h); digits = 4)), " )")
        write(file, " (m * m*s-2)\n")
        #
        write(file, "      vorticity: ")
        write(file, "( ", string(round(minimum(grid_vor); digits = 9)), ", ")
        write(file, string(round(maximum(grid_vor); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, "     divergence: ")
        write(file, "( ", string(round(minimum(grid_div); digits = 9)), ", ")
        write(file, string(round(maximum(grid_div); digits = 9)), " )")
        write(file, " (1/s)\n")
        write(file, repeat("***", 30), "\n")
    end
end