using DelimitedFiles
using Plots

# Lee un histograma guardado como CSV: (bin_center, count)
function read_hist_csv(path::String)
    M = readdlm(path, ',', Float64)
    x = M[:,1]
    c = M[:,2]
    return x, c
end

# Grafica como barras (conteos)
function plot_hist_csv(path::String; title::String="", xlabel::String="value")
    x, c = read_hist_csv(path)

    # ancho del bin (aprox)
    w = length(x) >= 2 ? (x[2]-x[1]) : 1.0

    bar(x, c;
        bar_width = 0.9*w,
        xlabel = xlabel,
        ylabel = "count",
        title  = title == "" ? basename(path) : title,
        legend = false)
end

function main(; outdir::String="hist_hw", save_png::Bool=true)
    files = [
        ("hist_sd_y.csv",        "sd_y (percent)",            "sd_y"),
        ("hist_sd_c_over_y.csv", "sd(c)/sd(y)",               "sd_c_over_y"),
        ("hist_sd_i_over_y.csv", "sd(i)/sd(y)",               "sd_i_over_y"),
        ("hist_sd_h_over_y.csv", "sd(h)/sd(y)",               "sd_h_over_y"),
        ("hist_sd_w_over_y.csv", "sd(w)/sd(y)",               "sd_w_over_y"),
        ("hist_sd_h_over_w.csv", "sd(h)/sd(w)",               "sd_h_over_w"),
        ("hist_cor_h_w.csv",     "corr(h,w)",                 "cor_h_w"),
    ]

    for (fname, ttl, xlab) in files
        path = joinpath(outdir, fname)
        p = plot_hist_csv(path; title=ttl, xlabel=xlab)
        display(p)
        if save_png
            pngpath = joinpath(outdir, replace(fname, ".csv" => ".png"))
            savefig(p, pngpath)
        end
    end

    println("✅ Histogramas mostrados. (y PNGs guardados en $(outdir)/ si save_png=true)")
end

# Si ejecutas: julia plot_hist_hw.jl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
