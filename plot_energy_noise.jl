using Plots
using CSV
using DelimitedFiles
file = raw"plot_data.txt"
M = readdlm(file)
phi_f = M[:,1];
noise = M[:,2];
s = M[:,3];
gr()
MyTitle = "L = 10, colormap = Entropy, Seed =4000, Page value = 5*log(2)-0.5 = 2.96";
p = plot(phi_f,noise,
    seriestype = :scatter,
    markerstrokecolor = "grey30",
    markerstrokewidth=0.0168,
    markersize=1.3,
    thickness_scaling = 1.4,
    xlims=(0,0.5), 
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
hline!([[-3.079082476590806]],lc=:green,legend=false)
hline!([ [0]],lc=:green,legend=false)
hline!([ [3.079082476590806]],lc=:green,legend=false)
xlabel!("Noise")
ylabel!("Energy")
savefig("4000_10.png")