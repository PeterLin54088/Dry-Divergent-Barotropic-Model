################################################################################
# Operators
export Compute_Legendre!
export Compute_Gaussian!
################################################################################

function Compute_Legendre!(num_fourier, num_spherical, sinθ, nθ)
    """
    Calculate associated Legendre polynomials and its derivative.
    See Spectral Numerical Weather Prediction Models Appendix B.
    """
    ########################################
    ### Main variable
    qnm = zeros(Float64, num_fourier + 1, num_spherical + 2, nθ)
    dqnm = zeros(Float64, num_fourier + 1, num_spherical + 1, nθ)

    ########################################
    ### Temporary variable
    ε = zeros(Float64, num_fourier + 1, num_spherical + 2)
    cosθ = sqrt.(1 .- sinθ .^ 2)
    # 
    for m in 0:num_fourier
        for l in m:(num_spherical + 1)
            ε[m + 1, l + 1] = sqrt((l^2 - m^2) ./ (4 * l^2 - 1))
        end
    end

    ########################################
    ### Legendre polynomials, P_{m,l}
    # P_{0,0} = 1
    qnm[1, 1, :] .= 1.0
    # P_{m,m}
    for m in 1:num_fourier
        qnm[m + 1, m + 1, :] = sqrt((2m + 1) / (2m)) .* cosθ .* qnm[m, m, :]
    end
    # P_{m+1,m}
    for m in 1:(num_fourier + 1)
        qnm[m, m + 1, :] = sqrt(2 * m + 1) * sinθ .* qnm[m, m, :]
    end
    # P_{m,l}
    for m in 0:num_fourier
        for l in (m + 2):(num_spherical + 1)
            qnm[m + 1, l + 1, :] = (sinθ .* qnm[m + 1, l, :] -
                                    ε[m + 1, l] * qnm[m + 1, l - 1, :]) / ε[m + 1, l + 1]
        end
    end

    ########################################
    ### Legendre polynomials with respect to μ, d(P_{m,l})/dμ
    for m in 0:num_fourier
        for l in m:num_spherical
            if l == m
                dqnm[m + 1, l + 1, :] = (-l * ε[m + 1, l + 2] * qnm[m + 1, l + 2, :]) ./
                                        (cosθ .^ 2)
            else
                dqnm[m + 1, l + 1, :] = (-l * ε[m + 1, l + 2] * qnm[m + 1, l + 2, :] +
                                         (l + 1) * ε[m + 1, l + 1] * qnm[m + 1, l, :]) ./
                                        (cosθ .^ 2)
            end
        end
    end
    return qnm[:, 1:(num_spherical + 1), :], dqnm
end

function Compute_Gaussian!(n)
    """
    Calculate weighting function for Gaussian quadrature.
    See Legendre-Gauss Quadrature on Wolfram.
    """
    ########################################
    ### Meta
    itermax = 10000
    tol = 1.0e-15

    ########################################
    ### Main variable
    sinθ = zeros(Float64, n)
    wts = zeros(Float64, n)

    ########################################
    ### Iteration (symmetric)
    n_half = Int64(n / 2)
    for i in 1:n_half
        # init at each latitude 
        dp = 0.0
        z = cos(pi * (i - 0.25) / (n + 0.5)) # first guess
        ### Nest loop ###
        for iter in 1:itermax
            # P_n
            p2 = 0.0
            p1 = 1.0
            for j in 1:n
                p3 = p2 # P_{j-2}
                p2 = p1 # P_{j-1}
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j # P_{j}
            end
            # P'_n
            dp = n * (z * p1 - p2) / (z * z - 1.0)
            # update wts
            z1 = z
            z = z1 - p1 / dp
            if (abs(z - z1) <= tol)
                break
            end
            if iter == itermax
                @error("Compute_Gaussian! does not converge!")
            end
        end
        ### Nest loop ###

        # 
        sinθ[i] = -z
        sinθ[n - i + 1] = z
        wts[i] = 2.0 / ((1.0 - z * z) * dp * dp)
        wts[n - i + 1] = 2.0 / ((1.0 - z * z) * dp * dp)
    end
    return sinθ, wts
end
