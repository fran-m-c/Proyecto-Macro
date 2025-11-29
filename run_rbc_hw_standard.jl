###############################################################
# run_rbc_hw_standard.jl  (NO RAW STRINGS, NO UNICODE, CLEAN)
###############################################################

using Dynare
using Random, LinearAlgebra, Statistics
using SparseArrays
using DataFrames, CSV


###############################
# HP FILTER
###############################

struct HPFilter
    T::Int
    λ::Float64
    Kfact::Cholesky{Float64, Matrix{Float64}}
end

function HPFilter(T::Int; λ::Real=1600.0)
    λ = float(λ)
    I = spdiagm(0 => ones(T))
    D = spzeros(T-2, T)
    for t in 1:(T-2)
        D[t, t] = 1.0
        D[t, t+1] = -2.0
        D[t, t+2] = 1.0
    end
    K = Matrix(I + λ*(D' * D))
    return HPFilter(T, λ, cholesky(Symmetric(K)))
end

function hp_cycle(hp::HPFilter, x::AbstractVector)
    trend = hp.Kfact \ Vector{Float64}(x)
    return Vector{Float64}(x) .- trend
end


###############################
# Extract A,B from Dynare
###############################

function extract_AB(context)
    mr = context.results.model_results[1]
    lre = mr.linearrationalexpectations
    A = Matrix(lre.g1_1)
    B = Matrix(lre.g1_2)
    return A, B
end


###############################
# Simulate deviations
###############################

function simulate_devs(A,B; T=100, burnin=200, rng=MersenneTwister(0), σ=1.0)
    n = size(A,1)
    m = size(B,2)
    yprev = zeros(n)

    for _ in 1:burnin
        ε = σ .* randn(rng, m)
        yprev = A*yprev + B*ε
    end

    Ydev = zeros(n,T)
    for t in 1:T
        ε = σ .* randn(rng, m)
        yprev = A*yprev + B*ε
        Ydev[:,t] = yprev
    end
    return Ydev
end


###############################
# Steady state from your derivation
###############################

function hw_steady_state()
    h = 1/3
    k = 12.6631
    y = exp(0)*k^0.36*h^(1-0.36)
    i = 0.025*k
    c = y - i
    z = 0.0
    ys = [c,i,y,k,h,z]
    names = ["c","i","y","k","h","z"]
    return ys, names
end


###############################
# Table-3 stats
###############################

function table3_stats_from_levels(hp, Y)
    c = vec(Y[1,:])
    i = vec(Y[2,:])
    y = vec(Y[3,:])
    h = vec(Y[5,:])
    w = y ./ h

    if any(x->x<=0, c) || any(x->x<=0, i) || any(x->x<=0, y) || any(x->x<=0, h) || any(x->x<=0, w)
        return nothing
    end

    cy_y = hp_cycle(hp, log.(y))
    cy_c = hp_cycle(hp, log.(c))
    cy_i = hp_cycle(hp, log.(i))
    cy_h = hp_cycle(hp, log.(h))
    cy_w = hp_cycle(hp, log.(w))

    sd_y = 100*std(cy_y)
    sd_c = 100*std(cy_c)
    sd_i = 100*std(cy_i)
    sd_h = 100*std(cy_h)
    sd_w = 100*std(cy_w)

    return (
        sd_y        = sd_y,
        sd_c_over_y = sd_c/sd_y,
        sd_i_over_y = sd_i/sd_y,
        sd_h_over_y = sd_h/sd_y,
        sd_w_over_y = sd_w/sd_y,
        sd_h_over_w = sd_h/sd_w,
        cor_h_w     = cor(cy_h, cy_w)
    )
end


###############################
# MAIN FUNCTION
###############################

function simulate_and_compute(; 
    modfile="rbc_hw_standard.mod", 
    T=179, 
    nsim=100,
    burnin=200,
    seed=1234,
    hp_lambda=1600.0
)

    # run Dynare
    context = Dynare.@dynare "/content/Proyecto-Macro/rbc_hw_standard.mod"

    A,B = extract_AB(context)
    ys, names = hw_steady_state()

    rng = MersenneTwister(seed)
    hp = HPFilter(T; λ=hp_lambda)

    draws = DataFrame(
        sd_y=Float64[], sd_c_over_y=Float64[], sd_i_over_y=Float64[],
        sd_h_over_y=Float64[], sd_w_over_y=Float64[],
        sd_h_over_w=Float64[], cor_h_w=Float64[]
    )

    for s in 1:nsim
        Ydev = simulate_devs(A,B; T=T, burnin=burnin, rng=rng)
        Y = ys .+ Ydev
        st = table3_stats_from_levels(hp,Y)
        st===nothing && continue
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
        nsim        = nsim
    )

    return draws, means
end


###############################
# MAIN ENTRY
###############################

function main()
    _, means = simulate_and_compute()
    println("=== RBC STANDARD TABLE-3 MEANS ===")
    println(means)
    return means
end

