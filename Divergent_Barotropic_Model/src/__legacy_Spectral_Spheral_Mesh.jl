"""
export
"""


function test_fourier()
    num_fourier, nθ, nd = 31, 48, 2
    num_spherical = num_fourier + 1
    nλ = 2nθ
    
    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, 1.0)
    
    # test 1 use num_fourier = num_spherical
    pfield0 = rand(Float64, nλ, nθ, nd)
    pfield1 = zeros(Float64, nλ, nθ, nd)
    pfield2 = zeros(Float64, nλ, nθ, nd)
    snm0 = zeros(Complex{Float64}, num_fourier + 1, num_fourier + 1, nd)
    snm1 = zeros(Complex{Float64}, num_fourier + 1, num_fourier + 1, nd)
    
    Trans_Grid_To_Spherical!(mesh, pfield0,  snm0)
    Trans_Spherical_To_Grid!(mesh, snm0,  pfield1)
    Trans_Grid_To_Spherical!(mesh, pfield1,  snm1)
    Trans_Spherical_To_Grid!(mesh, snm1,  pfield2)
    
    @info norm(snm0 - snm1), norm(snm0[:, 1:num_spherical,:] - snm1[:, 1:num_spherical,:])
    @info norm(pfield1 - pfield2)
end


function velocity_profile(IC::String = "rigid_rotation")
    # derivative will lose accuracy, which depends on the last two spherical modes)
    num_fourier, nθ, nd = 85, 128, 3 # 21, 32
    #num_fourier, nθ = 127, 128 # 21, 32
    num_spherical = num_fourier + 1
    nλ = 2nθ
    radius = 5.0
    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
    grid_u0, grid_v0 = zeros(Float64, nλ, nθ, nd), zeros(Float64, nλ, nθ, nd)
    grid_vor0, grid_div0 = rand(Float64, nλ, nθ, nd), rand(Float64, nλ, nθ, nd)
    cosθ = mesh.cosθ
    sinθ = mesh.sinθ
    if IC == "rigid_rotation"
        @info "Initialize the rigid_rotation initial condition"
        U, β = 2 * pi / 256.0, pi / 2.0
        λc = mesh.λc
        sinλc = sin.(λc)
        for k = 1:nd
            for i = 1:nλ
                grid_u0[i, :, k] .= U .* (cos(β) .* cosθ .+ sin(β) * cos(λc[i]) .* sinθ)
            end
            for j = 1:nθ
                grid_v0[:, j, k] .= -U * sin(β) .* sinλc
            end
            grid_div0 .= 0.0
            for i = 1:nλ
                grid_vor0[i, :, k] = 2 * U / radius * (cos(β) * sinθ - sin(β) * cos(λc[i]) * cosθ)
            end
        end
    elseif IC == "random_1"
        @info "Initialize the random initial condition"
        λc = mesh.λc
        sinλc = sin.(λc)
        for k = 1:nd
            for i = 1:nλ
                for j = 1:nθ
                    grid_u0[i, j, k] = cosθ[j]^2
                    grid_v0[i, j, k] = 0.0

                    grid_div0[i,j, k] = 0.0
                    grid_vor0[i,j, k] = 3*sinθ[j]*cosθ[j]/radius
                end
            end
        end
    elseif IC == "random_2"
        @info "Initialize the random initial condition"
        λc = mesh.λc
        sinλc = sin.(λc)
        for k = 1:nd
            for i = 1:nλ
                for j = 1:nθ
                    grid_u0[i, j, k] = 25 * cosθ[j] - 30 * cosθ[j]^3 + 300 * sinθ[j]^2 * cosθ[j]^6   # cosθ[j]^2
                    grid_v0[i, j, k] = 0.0
                    grid_div0[i,j, k] = 0.0
                    grid_vor0[i,j, k] = -(-50 * sinθ[j] + 120 * cosθ[j]^2 * sinθ[j] + 300 * (2cosθ[j]^7 * sinθ[j] - 7cosθ[j]^5 * sinθ[j]^3)) / radius # 3*sinθ[j]*cosθ[j]/radius
                end
            end
        end
    elseif IC == "mix"
        @info "Initialize the mix initial condition"
        U, β = 2 * pi / 256.0, pi / 2.0
        λc = mesh.λc
        sinλc = sin.(λc)
        for k = 1:nd
            for i = 1:nλ
                for j = 1:nθ
                    grid_u0[i, j, 1] = U * (cos(β) * cosθ[j] + sin(β) * cos(λc[i]) * sinθ[j])
                    grid_v0[i, j, 1] = -U * sin(β) * sinλc[i]
                    grid_u0[i, j, 2] = cosθ[j]^2
                    grid_v0[i, j, 2] = 0.0
                    grid_u0[i, j, 3] = 25 * cosθ[j] - 30 * cosθ[j]^3 + 300 * sinθ[j]^2 * cosθ[j]^6   # cosθ[j]^2
                    grid_v0[i, j, 3] = 0.0
                    grid_div0[i, j, 1] = 0.0
                    grid_vor0[i, j, 1] = 2 * U / radius * (cos(β) * sinθ[j] - sin(β) * cos(λc[i]) * cosθ[j])
                    grid_div0[i,j, 2] = 0.0
                    grid_vor0[i,j, 2] = 3*sinθ[j]*cosθ[j]/radius
                    grid_div0[i,j, 3] = 0.0
                    grid_vor0[i,j, 3] = -(-50 * sinθ[j] + 120 * cosθ[j]^2 * sinθ[j] + 300 * (2cosθ[j]^7 * sinθ[j] - 7cosθ[j]^5 * sinθ[j]^3)) / radius # 3*sinθ[j]*cosθ[j]/radius
                end
            end    
        end
    else
        error("Initial condition ", IC, " has not implemented yet")
    end
    return mesh, grid_u0, grid_v0, grid_vor0, grid_div0
end


function test_derivative(IC::String = "rigid_rotation")
    mesh, grid_u0, grid_v0, grid_vor0, grid_div0  = velocity_profile(IC)
    nλ, nθ, nd, num_fourier, num_spherical = mesh.nλ, mesh.nθ, mesh.nd, mesh.num_fourier, mesh.num_spherical 
    vor = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    div = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    grid_vor1 = zeros(Float64, nλ, nθ, nd)
    grid_div1 = zeros(Float64, nλ, nθ, nd)
    Vor_Div_From_Grid_UV!(mesh, grid_u0, grid_v0, vor, div)
    @info "norm(div), norm(vor)", norm(div[:,:,1]), norm(div[:,:,2]),norm(div[:,:,3]), norm(vor)
    Trans_Spherical_To_Grid!(mesh, div,  grid_div1)
    @info "velocity divergence error is ", norm(grid_div0 - grid_div1)
    Trans_Spherical_To_Grid!(mesh, vor,  grid_vor1)
    @info "velocity vorticity error is ", norm(grid_vor0 - grid_vor1)
    grid_u1 = zeros(Float64, nλ, nθ, nd)
    grid_v1 = zeros(Float64, nλ, nθ, nd)
    UV_Grid_From_Vor_Div!(mesh, vor, div, grid_u1, grid_v1)
    @info "u error ", norm(grid_u1 - grid_u0)
    @info "v error ", norm(grid_v1 - grid_v0)
end


function test_advection(IC::String = "rigid_rotation")
    mesh, grid_u0, grid_v0, _, _  = velocity_profile(IC)
    nλ, nθ, nd, num_fourier, num_spherical = mesh.nλ, mesh.nθ, mesh.nd, mesh.num_fourier, mesh.num_spherical 
    cosθ, sinθ = mesh.cosθ, mesh.sinθ
    radius = mesh.radius
    λc = mesh.λc
    grid_hs0 = zeros(Float64, nλ, nθ, nd)
    grid_∇hs0 = zeros(Float64, nλ, nθ, 2, nd)
    for k = 1:nd
        for i = 1:nλ
            for j = 1:nθ
                grid_hs0[i, j, k] = 25 * cosθ[j] - 30 * cosθ[j]^3 + 300 * sinθ[j]^2 * cosθ[j]^6 *sin(λc[i])
                grid_∇hs0[i, j, 1, k] = (300 * sinθ[j]^2 * cosθ[j]^5 *cos(λc[i]))/radius
                grid_∇hs0[i, j, 2, k] = (-25 * sinθ[j] + 90 * cosθ[j]^2*sinθ[j] + 600 * sinθ[j] * cosθ[j]^7*sin(λc[i]) - 1800 * sinθ[j]^3 * cosθ[j]^5*sin(λc[i]))/radius 
            end
        end
    end
    spe_hs0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    
    test_proj_error(mesh, grid_hs0)
    test_proj_error(mesh, grid_∇hs0[:,:,2,:])
    
    Trans_Grid_To_Spherical!(mesh, grid_hs0,  spe_hs0) 
    spe_cos_dλ_hs = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_cos_dθ_hs = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    
    spe_dθ_hs0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    Trans_Grid_To_Spherical!(mesh, grid_∇hs0[:,:,2,:],  spe_dθ_hs0) 
            
    Compute_Gradient_Cos!(mesh, spe_hs0, spe_cos_dλ_hs, spe_cos_dθ_hs)
    #@show norm(spe_cos_dθ_hs - spe_dθ_hs0.*cosθ')
            
    spe_vor0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_div0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_vor1 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_div1 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    grid_vor0 = zeros(Float64, nλ, nθ, nd)
    grid_div0 = zeros(Float64, nλ, nθ, nd)
    grid_vor1 = zeros(Float64, nλ, nθ, nd)
    grid_div1 = zeros(Float64, nλ, nθ, nd)
    
    grid_δhs = zeros(Float64, nλ, nθ, nd)
            
    Vor_Div_From_Grid_UV!(mesh, grid_u0, grid_v0, spe_vor0, spe_div0)
            
    Trans_Spherical_To_Grid!(mesh, spe_div0,  grid_div0)
    #grid_div1 = ∇((u,v)h) =   hdiv(u,v) + (u,v)∇h
    Vor_Div_From_Grid_UV!(mesh, grid_u0.*grid_hs0, grid_v0.*grid_hs0, spe_vor1, spe_div1)
    
    Trans_Spherical_To_Grid!(mesh, spe_div1,  grid_div1)
    #grid_δhs = -(u,v)∇h
    Add_Horizontal_Advection!(mesh, spe_hs0, grid_u0, grid_v0, grid_δhs)
    grid_δhs_ref = -(grid_∇hs0[:, :, 1,:] .*grid_u0 + grid_∇hs0[:, :, 2,:] .* grid_v0)
    
    @show "adv approach error: ", norm(grid_δhs_ref - grid_δhs)
    @show "divergence approach error: ", norm(grid_δhs_ref + grid_div1 + grid_hs0.*grid_div0)
end


if abspath(PROGRAM_FILE) == @__FILE__
    # TODO: cosθ is poorly represented on spherical harmonic basis         
    test_derivative("rigid_rotation")
    test_derivative("random_1")
    test_derivative("random_2")
    test_derivative("mix")

    test_advection("rigid_rotation")
    test_advection("random_1")
    test_advection("random_2")
    test_advection("mix")
end