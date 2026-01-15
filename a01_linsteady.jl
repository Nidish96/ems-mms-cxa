# * 
using DifferentialEquations
using NonlinearSolve
using LaTeXStrings
using ForwardDiff
using LinearAlgebra
using ColorSchemes

using GLMakie
using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
includet("./linearsys.jl");

savfigs = true;
if savfigs
    using CairoMakie
    CairoMakie.activate!();
else
    GLMakie.activate!();
end

cfgs = [(zt = 0.25e-2, w0 = 2.0, f = 1.0, Om=1.9),
        (zt = 25e-2, w0 = 2.0, f = 1.0, Om=1.9)];
ci = 2;

# * Compute
pars = cfgs[ci];

typs = [:MMS, :MMS, :EMS, :EMS];
ords = [1, 2, 1, 2];
sols = [];

u0 = [0.0,0.0];
Om0 = 0.05pars.w0;
Om1 = 2pars.w0;
dOm = 0.05pars.w0;
cpars = (save_jacs=true, parm=:arclength, Dsc=:none, nmax=10000);

for (i,ord) in enumerate(ords)
    fun = NonlinearFunction((r,u,p)->sflow!(r, u, (;pars...,
                                                   typ=typs[i],order=ord,Om=p)));

    sol, _, _, _, _ = CONTINUATE(u0, fun, [Om0, Om1], dOm; cpars..., verbosity=0);

    push!(sols, sol);

    display("Done $(typs[i]) O($(ords[i])) with $(length(sol)) points.")    
end

# * Plot
Oms = [s.up[end] for s in sols[1]];
Aans = [ansol((pars..., Om=om)) for om in Oms];
lstl = Linestyle([0, 5,6]);

set_theme!(theme_latexfonts())
fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)",
          ylabel="Forced Response Function (m/N)", yscale=log10);
mszs = [0, 0, [i%3==0 ? 18 : 0 for i in eachindex(sols[3])], 0];
styls = [:solid, lstl, :solid, lstl];
cols = colorschemes[:rainbow][[0.3,0.1,0.8,1.0]];
# cols = [:blue, :blue, :red, :red];

lines!(ax, Oms,  abs.(Aans), color=:black, linewidth=5,
    label="Exact");
for (i,sol) in enumerate(sols)
    if typs[i] == :MMS
        lbl = LaTeXString("$(typs[i]) \$\\mathcal{O}(\\varepsilon^$(ords[i]-1))\$");
    else
        lbl = LaTeXString("$(typs[i]) \$\\mathcal{O}(p^$(ords[i]-1))\$");
    end
    
    lines!(ax, [s.up[end] for s in sol],
           [norm(s.up[1:2]) for s in sol],
           linewidth=3, linestyle=styls[i],
           color=cols[i],
           label=lbl)
end

if ci==1
    axislegend(ax);
end


if Makie.current_backend()==GLMakie
   display(scr, fig);
else
    display(fig)   
end
   
if savfigs
    save("./FIGS/A01_linearsteady_C$(ci).pdf", fig)
end
