using DelimitedFiles
using Plots

function read_hist_csv(path::String)
    M = readdlm(path, ',', Float64; header=true)[1]
    return M[:,1], M[:,2]
end

function main(; outdir::String="hist_gov")
    files = [
        ("hist_sd_y.csv",        "sd_y (%)"),
        ("hist_sd_c_over_y.csv", "sd(c)/sd(y)"),
        ("hist_sd_i_over_y.csv", "sd(i)/sd(y)"),
        ("hist_sd_h_over_y.csv", "sd(h)/sd(y)"),
        ("hist_sd_w_over_y.csv", "sd(w)/sd(y)"),
        ("hist_sd_h_over_w.csv", "sd(h)/sd(w)"),
        ("hist_cor_h_w.csv",     "corr(h,w)"),
    ]

    for (fname, ttl) in files
        path = joinpath(outdir, fname)
        x, c = read_hist_csv(path)
        w = length(x) >= 2 ? (x[2]-x[1]) : 1.0

        p = bar(x, c; bar_width=0.9*w, xlabel=ttl, ylabel="count",
                title=ttl, legend=false)
        display(p)
        savefig(p, joinpath(outdir, replace(fname, ".csv" => ".png")))
    end

    println("✅ Listo: histogramas mostrados + PNG guardados en $(outdir)/")
end
