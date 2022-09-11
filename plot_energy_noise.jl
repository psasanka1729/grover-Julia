using Plots
using CSV
using DelimitedFiles
file = raw"plot_data_12_6000.txt.txt"
M = readdlm(file)
phi_f = M[:,1];
noise = M[:,2];
s = M[:,3];
gr()
MyTitle = "L = 12, colormap = Entropy, Seed =6000, Page value = 6*log(2)-0.5 = 3.66";
p = plot(phi_f,noise,
    seriestype = :scatter,
    markerstrokecolor = "grey30",
    markerstrokewidth=0.0168,
    markersize=1.0,
    thickness_scaling = 1.4,
    xlims=(0,0.3), 
    ylims=(-3.5,3.5),
    title = MyTitle,
    label = "",
    dpi=600,
    zcolor = s,
    color = :RdBu_4,
    right_margin = 5Plots.mm,
    left_margin = Plots.mm,
    titlefontsize=10,
    guidefontsize=10,
    tickfontsize=13,
    legendfontsize=10,
    )
plot!(size=(900,700))
hline!([[-3.11034138188]],lc=:green,legend=false)
hline!([ [0]],lc=:green,legend=false)
hline!([ [3.11034138188]],lc=:green,legend=false)
xlabel!("Noise")
ylabel!("Energy")
savefig("6000_12.png")
