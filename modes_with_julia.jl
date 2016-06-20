using PyCall
using Plots
using DataFrames
gr()

pth = "/home/castaned/Documents/ROTORCmodels/viscalc_suite/"
unshift!(PyVector(pyimport("sys")["path"]), pth)
@pyimport main_modules as mmod
@pyimport legendre_interp as lint

vels = 1:4
zps = Array{Float64}[]
mselect = "1p875"
ells = [6]

for vel in vels
  for m in ells
    amode = mmod.emode(mselect,Int(vel),"6 g4");
    xi_r,xi_t,dt_t,zg,r,zp,sig,cs = mmod.load_MODE_file([amode[:mass]],amode[:vel_idx],amode[:index],amode[:parity],amode[:frequency],true,false);
    if amode[:parity] == "EVEN"
          zp_fine = lint.leg_interp(zp[end-10:end,:],8,"EVEN")
      else
          zp_fine = lint.leg_interp(zp[end-10:end,:],8,"OE")
    end
    push!(zps,zp_fine[1,:])
  end
end

df = DataFrame()
for i in 1:length(zps)
  df[string("v",i)] = zps[i]
end
#using Plots
#gr()
# plotlyjs()
# pl1 = plot(linspace(0,90,length(zps[1][:])),zps[1][:],xlabel="R[zone]",ylabel="zp",legend=false,linewidth=1.5)
# plot!(linspace(0,90,length(zps[5][:])),zps[5][:],xlabel="R[zone]",ylabel="zp",legend=false,linewidth=1.5)
# plot!(linspace(0,90,length(zps[3][:])),zps[3][:],xlabel="R[zone]",ylabel="zp",legend=false,linewidth=1.5)
# PlotlyJS.relayout!(pl1.o, margin = Dict("l"=>100, "b"=>80, "r"=>8, "t"=>20))


# using PGFPlots
# xvec = 1:length(xi_r[:,1])
# plts = [Plots.Linear(1:length(xi_r[:,i]),xi_r[:,i],mark="none")::PGFPlots.Plots.Linear for i in [1,4,8]]
# p =Axis(plts,title = L"M=1.875M$\odot$ V=0 km s$^{-1}$",ylabel=L"$\xi_{r}$",xlabel="R[zone]")
# save("myfile.pdf", p)

# using Gadfly
# pl1=Gadfly.plot([layer(y=xi_r[:,i],x=1:length(xi_r[:,i]), Geom.line,Theme(line_width=0.6mm,default_color=Gadfly.distinguishable_colors(8)[i])) for i in [1,4,8]]...,Guide.xlabel("R[zone]"),Guide.ylabel("ξ_r"))
# draw(PDF("gadflyplot.pdf", 10inch, 6inch), pl1)
