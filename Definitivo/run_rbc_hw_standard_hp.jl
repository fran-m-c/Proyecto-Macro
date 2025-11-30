module RBC_HW_HP

using Dynare
using Random, LinearAlgebra, Statistics
using DelimitedFiles

# Orden del .mod: c i y k h z
const SS    = [0.918096, 0.3166, 1.2347, 12.6629, 0.3333333333, 0.0]
const IDX_C = 1
const IDX_I = 2
const IDX_Y = 3
const IDX_K = 4
const IDX_H = 5
const IDX_Z = 6

with_silent(f::Function) = redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        f()
    end
end

function run_dynare(modfile::String; quiet::Bool=true)
    expr = Meta.parse("Dynare.@dynare " * repr(modfile))
    return quiet ? with_silent(() -> eval(expr)) : eval(expr)
end

function normalize_ctx(ctx)
    if hasproperty(ctx, :results)
        return ctx
    elseif ctx isa AbstractVector
        for it in ctx
            (it !== nothing && hasproperty(it, :results)) && return it
        end
    end
    error("Dynare devolvió algo sin `.results`. Revisa que el .mod corra bien.")
end

function dynare_AB(modfile::String; quiet::Bool=true)
    ctx0 = run_dynare(modfile; quiet=quiet)
    ctx  = normalize_ctx(ctx0)
    mr   = ctx.results.model_results[1]
    lre  = mr.linearrationalexpectations
    A = Matrix(lre.g1_1)      # nvar × nstate
    B = vec(Matrix(lre.g1_2)) # nvar (1 shock)
    return A, B
end

# ---------- HP filter (HP sobre logs), λ=1600 ----------
struct HP
    Kfact::Cholesky{Float64, Matrix{Float64}}
end

function HP(T::Int; λ::Float64=1600.0)
    Iden = Matrix{Float64}(I, T, T)
    D = zeros(T-2, T)
    @inbounds for t in 1:(T-2)
        D[t,t]   = 1.0
        D[t,t+1] = -2.0
        D[t,t+2] = 1.0
    end
    K = Iden + λ*(D' * D)
    return HP(cholesky(Symmetric(K)))
end

# OJO: acá es \ (una sola)
hp_cycle(hp::HP, x::Vector{Float64}) = x .- (hp.Kfact \ x)

function simulate_one(A::AbstractMatrix, B::AbstractVector, rng::AbstractRNG;
                      T::Int=200, burn::Int=500, sigma::Float64=0.007)

    nstate = size(A,2)
    nstate == 2 || error("Esperaba 2 estados (k,z). Dynare devolvió nstate=" * string(nstate))

    col_z = argmax(abs.(A[IDX_Z, :]))
    col_k = 3 - col_z

    s = zeros(2)

    for _ in 1:burn
        u = sigma * randn(rng)
        xdev = A*s .+ B*u
        s_next = zeros(2)
        s_next[col_k] = xdev[IDX_K]
        s_next[col_z] = xdev[IDX_Z]
        s = s_next
    end

    Ysim = zeros(size(A,1), T)
    for t in 1:T
        u = sigma * randn(rng)
        xdev = A*s .+ B*u
        Ysim[:,t] = SS .+ xdev

        s_next = zeros(2)
        s_next[col_k] = xdev[IDX_K]
        s_next[col_z] = xdev[IDX_Z]
        s = s_next
    end

    return Ysim
end

function moments_table3(hp::HP, Ysim::AbstractMatrix)
    y = vec(Ysim[IDX_Y,:])
    c = vec(Ysim[IDX_C,:])
    i = vec(Ysim[IDX_I,:])
    h = vec(Ysim[IDX_H,:])
    w = y ./ h

    if any(y .<= 0) || any(c .<= 0) || any(i .<= 0) || any(h .<= 0) || any(w .<= 0)
        return nothing
    end

    cy_y = hp_cycle(hp, log.(Float64.(y)))
    cy_c = hp_cycle(hp, log.(Float64.(c)))
    cy_i = hp_cycle(hp, log.(Float64.(i)))
    cy_h = hp_cycle(hp, log.(Float64.(h)))
    cy_w = hp_cycle(hp, log.(Float64.(w)))

    return (
        sd_y        = 100 * std(cy_y),
        sd_c_over_y = std(cy_c) / std(cy_y),
        sd_i_over_y = std(cy_i) / std(cy_y),
        sd_h_over_y = std(cy_h) / std(cy_y),
        sd_w_over_y = std(cy_w) / std(cy_y),
        sd_h_over_w = std(cy_h) / std(cy_w),
        cor_h_w     = cor(cy_h, cy_w)
    )
end

function save_hist_csv(x::Vector{Float64}, statname::String; nbins::Int=40, outdir::String="hist_hw")
    mkpath(outdir)
    xmin, xmax = minimum(x), maximum(x)
    edges = range(xmin, xmax; length=nbins+1)
    counts = zeros(Int, nbins)
    for v in x
        b = clamp(searchsortedlast(edges, v), 1, nbins)
        counts[b] += 1
    end
    centers = (collect(edges[1:end-1]) .+ collect(edges[2:end])) ./ 2
    out = hcat(centers, counts)
    writedlm(joinpath(outdir, "hist_" * statname * ".csv"), out, ',')
    return nothing
end

function tarea(; modfile::String="rbc_hw_standard.mod",
                T::Int=200, nsim::Int=10_000, burn::Int=500,
                sigma::Float64=0.007, seed::Int=1234,
                hp_lambda::Float64=1600.0, quiet_dynare::Bool=true,
                outdir::String="hist_hw",
                save_one_path::Bool=true)

    A, B = dynare_AB(modfile; quiet=quiet_dynare)
    hp   = HP(T; λ=hp_lambda)
    rng  = MersenneTwister(seed)

    if save_one_path
        Y1 = simulate_one(A,B,rng; T=T, burn=burn, sigma=sigma)
        writedlm("one_path_200.csv", permutedims(Y1), ',')
    end

    sd_y        = Float64[]; sd_c_over_y = Float64[]; sd_i_over_y = Float64[]
    sd_h_over_y = Float64[]; sd_w_over_y = Float64[]; sd_h_over_w = Float64[]
    cor_h_w     = Float64[]

    for _ in 1:nsim
        Ysim = simulate_one(A,B,rng; T=T, burn=burn, sigma=sigma)
        m = moments_table3(hp, Ysim)
        m === nothing && continue
        push!(sd_y,        m.sd_y)
        push!(sd_c_over_y, m.sd_c_over_y)
        push!(sd_i_over_y, m.sd_i_over_y)
        push!(sd_h_over_y, m.sd_h_over_y)
        push!(sd_w_over_y, m.sd_w_over_y)
        push!(sd_h_over_w, m.sd_h_over_w)
        push!(cor_h_w,     m.cor_h_w)
    end

    draws = hcat(sd_y, sd_c_over_y, sd_i_over_y, sd_h_over_y, sd_w_over_y, sd_h_over_w, cor_h_w)
    writedlm("table3_draws_10000.csv", draws, ',')

    out = (
        sd_y        = mean(sd_y),
        sd_c_over_y = mean(sd_c_over_y),
        sd_i_over_y = mean(sd_i_over_y),
        sd_h_over_y = mean(sd_h_over_y),
        sd_w_over_y = mean(sd_w_over_y),
        sd_h_over_w = mean(sd_h_over_w),
        cor_h_w     = mean(cor_h_w),
        n_used      = length(sd_y)
    )

    save_hist_csv(sd_y,        "sd_y";        outdir=outdir)
    save_hist_csv(sd_c_over_y, "sd_c_over_y"; outdir=outdir)
    save_hist_csv(sd_i_over_y, "sd_i_over_y"; outdir=outdir)
    save_hist_csv(sd_h_over_y, "sd_h_over_y"; outdir=outdir)
    save_hist_csv(sd_w_over_y, "sd_w_over_y"; outdir=outdir)
    save_hist_csv(sd_h_over_w, "sd_h_over_w"; outdir=outdir)
    save_hist_csv(cor_h_w,     "cor_h_w";     outdir=outdir)

    println(out)
    println("✅ guardé one_path_200.csv, table3_draws_10000.csv y los histogramas CSV en: " * outdir * "/")
    return out
end

end # module
