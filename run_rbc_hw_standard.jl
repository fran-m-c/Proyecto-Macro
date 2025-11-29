using Dynare
using MAT
using LinearAlgebra
using Random, Statistics
using DataFrames, CSV

############################################
# 1. Ejecutar Dynare (modo clásico MATLAB)
############################################

function run_dynare()
    Dynare.@dynare "rbc_hw_standard.mod"
end

############################################
# 2. Cargar matrices desde oo_.mat
############################################

function load_solution()
    file = matread("rbc_hw_standard/results/rbc_hw_standard_results.mat")

    dr = file["oo_"]["dr"]

    A = dr["ghx"]    # state transition
    B = dr["ghu"]    # shock impact

    ys = dr["ys"]    # steady state

    return A, B, ys
end

############################################
# 3. Simulación lineal
############################################

function simulate(A, B, ys; T=200, burn=200)
    rng = MersenneTwister(1234)
    n = size(A,1)
    m = size(B,2)

    y = zeros(n)
    for _ in 1:burn
        eps = randn(rng,m)
        y = A*y + B*eps
    end

    Y = zeros(n,T)
    for t in 1:T
        eps = randn(rng,m)
        y = A*y + B*eps
        Y[:,t] = y
    end

    return Y .+ ys
end

############################################
# 4. HP filter simple
############################################

function hp_cycle(x; λ=1600)
    T = length(x)
    D = zeros(T-2, T)
    for i in 1:T-2
        D[i,i] = 1
        D[i,i+1] = -2
        D[i,i+2] = 1
    end
    K = I + λ*(D' * D)
    trend = K \ x
    return x - trend
end

############################################
# 5. Estadísticos tipo Tabla 3
############################################

function stats(Y)
    c = vec(Y[1,:])
    i = vec(Y[2,:])
    y = vec(Y[3,:])
    h = vec(Y[5,:])

    w = y ./ h

    cy_y = hp_cycle(log.(y))
    cy_c = hp_cycle(log.(c))
    cy_i = hp_cycle(log.(i))
    cy_h = hp_cycle(log.(h))
    cy_w = hp_cycle(log.(w))

    return (
        sd_y        = 100*std(cy_y),
        sd_c_over_y = std(cy_c)/std(cy_y),
        sd_i_over_y = std(cy_i)/std(cy_y),
        sd_h_over_y = std(cy_h)/std(cy_y),
        sd_w_over_y = std(cy_w)/std(cy_y),
        sd_h_over_w = std(cy_h)/std(cy_w),
        cor_h_w     = cor(cy_h, cy_w)
    )
end

############################################
# 6. MAIN
###################
