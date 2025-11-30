module RBC_HW_GOV_HP

using Dynare
using Random, LinearAlgebra, Statistics
using Printf

# Orden del .mod: c i y k h z g
const SS   = [0.6465, 0.3166, 1.2347, 12.6629, 0.3333333333, 0.0, 0.271634]
const IDX_C = 1
const IDX_I = 2
const IDX_Y = 3
const IDX_K = 4
const IDX_H = 5
const IDX_Z = 6
const IDX_G = 7

with_silent(f::Function) = redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        f()
    end
end

# Llamar @dynare con un string variable (evita errores Symbol/modfile)
function run_dynare(modfile::String; quiet::Bool=true)
    expr = Meta.parse("Dynare.@dynare " * repr(modfile))
    return quiet ? with_silent(() -> eval(expr)) : eval(expr)
end

function dynare_AB(modfile::String; quiet::Bool=true)
    ctx = run_dynare(modfile; quiet=quiet)
    mr  = ctx.results.model_results[1]
    lre = mr.linearrationalexpectations
    A = Matrix(lre.g1_1)          # nvar × nstate
    B = Matrix(lre.g1_2)          # nvar × nshock (acá 2)
    return A, B
end

# ---------- HP filter (λ=1600) ----------
struct HP
    Kfact::Cholesky{Float64, Matrix{Float64}}
end

function HP(T::Int; hp_lambda::Float64=1600.0)
    Iden = Matrix(I, T, T)
    D = zeros(T-2, T)
    @inbounds for t in 1:(T-2)
        D[t,t]   = 1.0
        D[t,t+1] = -2.0
        D[t,t+2] = 1.0
    end
    K = Iden + hp_lambda*(D' * D)
    return HP(cholesky(Symmetric(K)))
end

hp_cycle(hp::HP, x::Vector{Float64}) = x .- (hp.Kfact \ x)

# ---------- asignar columnas de estado (k, z, g) ----------
function assign_state_cols(A::AbstractMatrix)
    nstate = size(A,2)
    nstate == 3 || error("Esperaba 3 estados (k,z,g). Dynare devolvió nstate=$(nstate)")

    rows = [IDX_K, IDX_Z, IDX_G]
    perms = ((1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1))

    best = nothing
    bestscore = -Inf
    for p in perms
        score = abs(A[rows[1], p[1]]) + abs(A[rows[2], p[2]]) + abs(A[rows[3], p[3]])
        if score > bestscore
            bestscore = score
            best = p
        end
    end
    col_k, col_z, col_g = best
    return col_k, col_z, col_g
end

# ---------- Simulación ----------
function simulate_one(A::AbstractMatrix, B::AbstractMatrix;
                      T::Int=200, burn::Int=500, seed::Int=1234,
                      sigma_eps::Float64=0.007, sigma_mu::Float64=0.021)

    rng = MersenneTwister(seed)
    col_k, col_z, col_g = assign_state_cols(A)

    s = zeros(3)  # estados internos Dynare (mapeados por columnas)

    for _ in 1:burn
        u = [sigma_eps*randn(rng), sigma_mu*randn(rng)]
        xdev = A*s .+ B*u

        s_next = zeros(3)
        s_next[col_k] = xdev[IDX_K]
        s_next[col_z] = xdev[IDX_Z]
        s_next[col_g] = xdev[IDX_G]
        s = s_next
    end

    Ysim = zeros(size(A,1), T)
    for t in 1:T
        u = [sigma_eps*randn(rng), sigma_mu*randn(rng)]
        xdev = A*s .+ B*u
        Ysim[:,t] = SS .+ xdev

        s_next = zeros(3)
        s_next[col_k] = xdev[IDX_K]
        s_next[col_z] = xdev[IDX_Z]
        s_next[col_g] = xdev[IDX_G]
        s = s_next
    end

    return Ysim
end

# ---------- Momentos tipo Tabla 3 (HP sobre logs) ----------
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

    sd_y = 100 * std(cy_y)

    return (
        sd_y        = sd_y,
        sd_c_over_y = std(cy_c) / std(cy_y),
        sd_i_over_y = std(cy_i) / std(cy_y),
        sd_h_over_y = std(cy_h) / std(cy_y),
        sd_w_over_y = std(cy_w) / std(cy_y),
        sd_h_over_w = std(cy_h) / std(cy_w),
        cor_h_w     = cor(cy_h, cy_w)
    )
end

# ---------- CSV helpers ----------
function write_csv(path::String, header::Vector{String}, rows::Vector{Vector{Float64}})
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join(r, ","))
        end
    end
end

function hist_csv(x::Vector{Float64}; nbins::Int=40)
    xmin = minimum(x); xmax = maximum(x)
    if xmax == xmin
        return [xmin], [length(x)]
    end
    w = (xmax - xmin) / nbins
    counts = zeros(Int, nbins)
    @inbounds for v in x
        b = Int(floor((v - xmin) / w)) + 1
        b = clamp(b, 1, nbins)
        counts[b] += 1
    end
    centers = [xmin + (j-0.5)*w for j in 1:nbins]
    return centers, counts
end

function save_histograms(outdir::String, draws::Vector{NamedTuple}; nbins::Int=40)
    mkpath(outdir)

    specs = (
        ("sd_y",        collect(getfield.(draws, :sd_y))),
        ("sd_c_over_y", collect(getfield.(draws, :sd_c_over_y))),
        ("sd_i_over_y", collect(getfield.(draws, :sd_i_over_y))),
        ("sd_h_over_y", collect(getfield.(draws, :sd_h_over_y))),
        ("sd_w_over_y", collect(getfield.(draws, :sd_w_over_y))),
        ("sd_h_over_w", collect(getfield.(draws, :sd_h_over_w))),
        ("cor_h_w",     collect(getfield.(draws, :cor_h_w))),
    )

    for (nm, vecx) in specs
        centers, counts = hist_csv(vecx; nbins=nbins)
        rows = [ [centers[j], Float64(counts[j])] for j in eachindex(centers) ]
        write_csv(joinpath(outdir, "hist_"*nm*".csv"), ["bin_center","count"], rows)
    end
end

# ---------- TAREA (pasos 4–6) ----------
function tarea(; modfile::String="rbc_hw_gov.mod",
                 T::Int=200, nsim::Int=10_000, burn::Int=500,
                 sigma_eps::Float64=0.007, sigma_mu::Float64=0.021,
                 seed::Int=1234, hp_lambda::Float64=1600.0,
                 quiet_dynare::Bool=true,
                 out_onepath::String="one_path_200.csv",
                 out_draws::String="table3_draws_10000.csv",
                 out_hist_dir::String="hist_gov",
                 nbins::Int=40)

    @printf("=== HW(1992) RBC + Gov Spending ===\n")
    @printf("Running Dynare on: %s\n", modfile)

    A, B = dynare_AB(modfile; quiet=quiet_dynare)
    hp = HP(T; hp_lambda=hp_lambda)

    # Paso 4: un camino
    Y = simulate_one(A, B; T=T, burn=burn, seed=seed, sigma_eps=sigma_eps, sigma_mu=sigma_mu)
    y = vec(Y[IDX_Y,:]); h = vec(Y[IDX_H,:]); w = y ./ h

    open(out_onepath, "w") do io
        println(io, "t,c,i,y,k,h,z,g,w")
        for t in 1:T
            println(io, string(t, ",",
                              Y[IDX_C,t], ",", Y[IDX_I,t], ",", Y[IDX_Y,t], ",",
                              Y[IDX_K,t], ",", Y[IDX_H,t], ",", Y[IDX_Z,t], ",",
                              Y[IDX_G,t], ",", w[t]))
        end
    end

    # Paso 5: 10.000 simulaciones (momentos Tabla 3)
    keep = NamedTuple[]
    for s in 1:nsim
        Ys = simulate_one(A, B; T=T, burn=burn, seed=seed+s,
                          sigma_eps=sigma_eps, sigma_mu=sigma_mu)
        m = moments_table3(hp, Ys)
        m === nothing && continue
        push!(keep, m)
    end

    # guardar draws
    rows = Vector{Vector{Float64}}(undef, length(keep))
    @inbounds for j in eachindex(keep)
        m = keep[j]
        rows[j] = [m.sd_y, m.sd_c_over_y, m.sd_i_over_y, m.sd_h_over_y, m.sd_w_over_y, m.sd_h_over_w, m.cor_h_w]
    end
    write_csv(out_draws,
              ["sd_y","sd_c_over_y","sd_i_over_y","sd_h_over_y","sd_w_over_y","sd_h_over_w","cor_h_w"],
              rows)

    # Paso 6: histogramas (CSV)
    save_histograms(out_hist_dir, keep; nbins=nbins)

    out = (
        sd_y        = mean(getfield.(keep, :sd_y)),
        sd_c_over_y = mean(getfield.(keep, :sd_c_over_y)),
        sd_i_over_y = mean(getfield.(keep, :sd_i_over_y)),
        sd_h_over_y = mean(getfield.(keep, :sd_h_over_y)),
        sd_w_over_y = mean(getfield.(keep, :sd_w_over_y)),
        sd_h_over_w = mean(getfield.(keep, :sd_h_over_w)),
        cor_h_w     = mean(getfield.(keep, :cor_h_w)),
        n_used      = length(keep)
    )

    println(out)
    println("✅ guardé: $(out_onepath), $(out_draws) y histogramas CSV en: $(out_hist_dir)/")
    return out
end

main = tarea

end # module
