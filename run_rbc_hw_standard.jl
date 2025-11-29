# run_rbc_hw_standard.jl
# Simulación de modelo RBC estándar tipo Tabla 3 (Hansen–Wright)
# usando Dynare.jl para obtener la solución lineal.

using Dynare
using Random, LinearAlgebra, Statistics
using SparseArrays
using DataFrames, CSV

# ==================== HP filter rápido ====================

struct HPFilter
    T::Int
    λ::Float64
    Kfact::Cholesky{Float64, Matrix{Float64}}
end

function HPFilter(T::Int; λ::Real=1600.0)
    λ = float(λ)
    I = spdiagm(0 => ones(T))
    D = spzeros(T-2, T)
    @inbounds for t in 1:(T-2)
        D[t, t]   = 1.0
        D[t, t+1] = -2.0
        D[t, t+2] = 1.0
    end
    K = Matrix(I + λ*(D' * D))
    return HPFilter(T, λ, cholesky(Symmetric(K)))
end

@inline function hp_cycle(hp::HPFilter, x::AbstractVector{<:Real})
    trend = hp.Kfact \ Vector{Float64}(x)
    return Vector{Float64}(x) .- trend
end

# ==================== Utilidades ====================

# Extraer matrices de decisión lineal A,B de Dynare
# (equivalente a ghx/ghu)
function extract_AB(context)
    mr = context.results.model_results[1]
    lre = mr.linearrationalexpectations
    A = Matrix(lre.g1_1)           # n×n
    B = Matrix(lre.g1_2)           # n×m
    return A, B
end

# Simulación de desviaciones respecto al steady state
function simulate_devs(A::AbstractMatrix, B::AbstractMatrix;
                      T::Int, burnin::Int, rng::AbstractRNG, σ::Real=1.0)
    n = size(A, 1)
    m = size(B, 2)
    yprev = zeros(n)

    # burn-in
    for _ in 1:burnin
        ε = σ .* randn(rng, m)
        yprev = A*yprev + B*ε
    end

    Ydev = zeros(n, T)
    @inbounds for t in 1:T
        ε = σ .* randn(rng, m)
        yprev = A*yprev + B*ε
        Ydev[:, t] = yprev
    end
    return Ydev
end

# ==================== Steady state según tu desarrollo ====================

"""
Devuelve (ys, names) con el steady state en este orden de variables:
c, i, y, k, h, z
"""
function hw_steady_state()
    β      = 0.99
    δ      = 0.025
    θ      = 0.36
    hbar   = 1/3
    ρ      = 0.95  # solo para info
    σ_ε    = 0.007 # solo para info

    # 1) Euler estacionario: 1 = β (θ (k/h)^{θ-1} + 1-δ)
    #    => (k/h)^{-0.64} = ...  (ya lo hiciste)
    kh_ratio = 37.9893           # tu resultado
    kbar     = kh_ratio * hbar   # 12.6631 aprox

    zbar = 0.0
    ybar = exp(zbar) * kbar^θ * hbar^(1-θ)
    ibar = δ * kbar
    cbar = ybar - ibar

    ys = [cbar, ibar, ybar, kbar, hbar, zbar]
    names = ["c","i","y","k","h","z"]

    return ys, names
end

# ==================== Estadísticas tipo Tabla 3 ====================

function table3_stats_from_levels(hp::HPFilter, Y::AbstractMatrix{<:Real})
    # Orden de variables: c i y k h z
    c = vec(Y[1, :])
    i = vec(Y[2, :])
    y = vec(Y[3, :])
    h = vec(Y[5, :])

    w = y ./ h  # productividad como y/h

    # comprobar positividad para log
    if any(c .<= 0) || any(i .<= 0) || any(y .<= 0) || any(h .<= 0) || any(w .<= 0)
        return nothing
    end

    cy_y = hp_cycle(hp, log.(y))
    cy_c = hp_cycle(hp, log.(c))
    cy_i = hp_cycle(hp, log.(i))
    cy_h = hp_cycle(hp, log.(h))
    cy_w = hp_cycle(hp, log.(w))

    sd_y = 100 * std(cy_y)
    sd_c = 100 * std(cy_c)
    sd_i = 100 * std(cy_i)
    sd_h = 100 * std(cy_h)
    sd_w = 100 * std(cy_w)

    return (
        sd_y        = sd_y,
        sd_c_over_y = sd_c / sd_y,
        sd_i_over_y = sd_i / sd_y,
        sd_h_over_y = sd_h / sd_y,
        sd_w_over_y = sd_w / sd_y,
        sd_h_over_w = sd_h / sd_w,
        cor_h_w     = cor(cy_h, cy_w),
    )
end

# ==================== Cálculo principal ====================

function simulate_and_compute(; modfile::String="rbc_hw_standard.mod",
                                T::Int=179,
                                nsim::Int=100,
                                burnin::Int=200,
                                seed::Int=1234,
                                hp_lambda::Real=1600.0)

    # 1. Resolver el modelo en Dynare
    context = Dynare.@dynare modfile
    A, B    = extract_AB(context)

    # 2. Steady state (vector ys)
    ys, _   = hw_steady_state()

    rng = MersenneTwister(seed)
    hp  = HPFilter(T; λ=hp_lambda)

    draws = DataFrame(
        sd_y        = Float64[],
        sd_c_over_y = Float64[],
        sd_i_over_y = Float64[],
        sd_h_over_y = Float64[],
        sd_w_over_y = Float64[],
        sd_h_over_w = Float64[],
        cor_h_w     = Float64[],
    )

    for s in 1:nsim
        Ydev = simulate_devs(A, B; T=T, burnin=burnin, rng=rng, σ=1.0)
        Y    = ys .+ Ydev
        st   = table3_stats_from_levels(hp, Y)
        st === nothing && continue
        push!(draws, st)
    end

    means = DataFrame(
        sd_y        = mean(draws.sd_y),
        sd_c_over_y = mean(draws.sd_c_over_y),
        sd_i_over_y = mean(draws.sd_i_over_y),
        sd_h_over_y = mean(draws.sd_h_over_y),
        sd_w_over_y = mean(draws.sd_w_over_y),
        sd_h_over_w = mean(draws.sd_h_over_w),
        cor_h_w     = mean(draws.cor_h_w),
        n_used      = nrow(draws),
        T           = T,
        nsim        = nsim,
        burnin      = burnin,
        hp_lambda   = float(hp_lambda),
    )

    CSV.write("rbc_hw_draws.csv", draws)
    CSV.write("rbc_hw_means.csv", means)
