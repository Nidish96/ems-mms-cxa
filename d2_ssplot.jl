# * Preamble
using DifferentialEquations
using NonlinearSolve
using LaTeXStrings
using LinearAlgebra
using JLD2
using GLMakie

using Revise

using juliajim.HARMONIC
using juliajim.CONTINUATION
includet("./vdpsys.jl")

analyze = false;
savdats = true;

pars = (c=0.04, w0=2., al=0.1, f=1., Om=0.1);

# * Save figure boolean
savfigs = true;
if savfigs
    using CairoMakie
    CairoMakie.activate!();
else
    GLMakie.activate!();
end

# * Load Data
trdat = load("./DATS/D0_trdat.jld2");
trdat = (; [Symbol(k)=>v for (k,v) in trdat]...);

typs = [:MMS, :MMS, :EMS, :EMS];
ords = [1, 2, 1, 2];

sldat = [load("./DATS/D1_$(typ)_$(ord).jld2") for (typ,ord) in zip(typs,ords)];
sldat = [(; [Symbol(k)=>v for (k,v) in sld]...) for sld in sldat][:];

# * Plot Separately
lstl = Linestyle([0,5,6]);
Om_p = 1.675;

set_theme!(theme_latexfonts())
fsz = 22;
lw = 3;
fig = Figure(fontsize=fsz, size=(1020, 820));
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

axs = [];
for i in 1:2
    for j in 1:2
        ax = Axis(fig[i, j]);
        if j==1
            ax.ylabel="Response Amplitude (m)"
        end
        if i==2
            ax.xlabel="Excitation Frequency (rad/s)";
        end
        push!(axs, ax)
    end
end
for (ti, ax) in enumerate(axs)
    scatter!(ax, trdat.Oms, trdat.Umaxs, color=:black, label="Time Integration")

    lines!(ax, sldat[ti].mainb.Oms./(sldat[ti].mainb.stab.==0),
        sldat[ti].mainb.Umaxs, color=:blue, linewidth=lw,
        label="Periodic, Stable")
    lines!(ax, sldat[ti].mainb.Oms./(sldat[ti].mainb.stab.==2),
        sldat[ti].mainb.Umaxs, color=:red, linewidth=lw,
        label="Periodic, Unstable")    
    lines!(ax, sldat[ti].bifbs[1].Oms, sldat[ti].bifbs[1].Umaxs,
        color=:orange, linewidth=lw, label="QuasiPeriodic")
    lines!(ax, sldat[ti].bifbs[2].Oms, sldat[ti].bifbs[2].Umaxs,
        color=:orange, linewidth=lw)

    if typs[ti] == :MMS
        ax.title = LaTeXString("$(typs[ti]) \$\\mathcal{O}(\\varepsilon^$(ords[ti]-1))\$");
    else
        ax.title = LaTeXString("$(typs[ti]) \$\\mathcal{O}(p^$(ords[ti]-1))\$");
    end


    vlines!(ax, Om_p, color=:gray)
    text!(ax, 1.02Om_p, 0, text=L"\Omega_p", align=(:left, :bottom))
    ylims!(ax, 0, 3)
end
Legend(fig[3,1:2], axs[1], nbanks=4, tellheight=true)
linkaxes!(axs...)

xlims!(axs[4], 1.2, 2.8)

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
   
if savfigs
   save("./FIGS/Bd2_FvdpResps.pdf", fig)
end

# * Transient Analysis at Chosen Omega

ncyc = 60;
nt = 1024;
tper = 2π/Om_p;
t = range(0, tper*ncyc, nt*ncyc);

prob = ODEProblem(rocfun!, zeros(2), [0, tper*ncyc], (;pars..., Om=Om_p));
sol = solve(prob);
sol = sol(t);

typs = [:MMS, :MMS, :MMS, :EMS, :EMS, :EMS];
ords = [1, 2, 3, 1, 2, 3];
ssols = [];
for (typ,ord) in zip(typs, ords)
    sprob = ODEProblem(sflow!, zeros(2), [0, tper*ncyc],
    (;pars..., Om=Om_p,typ=typ,order=ord));
    ssol = solve(sprob);
    push!(ssols, ssol(t));
end

# ** Plot Transient Responses
lstl = Linestyle([0, 5, 6]);

set_theme!(theme_latexfonts())
fsz = 22;
fig1 = Figure(fontsize=fsz, size=(1120,720));
if !isdefined(Main, :scr1) && Makie.current_backend()==GLMakie
   scr1 = GLMakie.Screen();
end

axs = [];
for i in 1:2
    for j in 1:3
        ax = Axis(fig1[i, j]);
        if j==1
            ax.ylabel="Response (m)"
        end
        if i==2
            ax.xlabel="Time (s)";
        end
        push!(axs, ax)
    end
end

for (ax,typ,ord,ssol) in zip(axs,typs,ords,ssols)    
    lines!(ax, sol.t, sol[1,:], color=:gray, label="Time Integration of Second Order System")
    lines!(ax, t, norm.(ssol.u), linewidth=2, color=:blue,
        label="Time integration of Slow Flow System")

    if typ==:MMS
        ax.title = LaTeXString("$(typ) \$\\mathcal{O}(\\varepsilon^$(ord-1))\$");
    else
        ax.title = LaTeXString("$(typ) \$\\mathcal{O}(p^$(ord-1))\$");
    end

end

Legend(fig1[3,1:3], axs[1], nbanks=2, tellheight=true);
linkaxes!(axs...)

# xlims!(ax, tper*ncyc/2, tper*ncyc)

if Makie.current_backend()==GLMakie
   display(scr1, fig1);
else
   display(fig1)
end
   
if savfigs
   save("./FIGS/Bd2_transsols.pdf", fig1)
end
