# * 
using DifferentialEquations
using NonlinearSolve
using LaTeXStrings
using ForwardDiff

using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
includet("./linearsys.jl");

savfigs = false;
if savfigs
    using CairoMakie
    CairoMakie.activate!();
else
    using GLMakie
    GLMakie.activate!();
end

cfgs = [(zt = 0.25e-2, w0 = 2.0, f = 1.0, Om=1.9),
        (zt = 25e-2, w0 = 2.0, f = 1.0, Om=1.9)];
ci = 1;

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

# Highlights for transient analysis
pts_p = ["A", "B", "C"];;
Oms_p = [0.95pars.w0, 1.5pars.w0, 0.1pars.w0];

set_theme!(theme_latexfonts())
fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)",
          ylabel="Forced Response Function (m/N)", yscale=log10);
mszs = [0, 0, [i%3==0 ? 18 : 0 for i in eachindex(sols[3])], 0];
styls = [:solid, :dash, :solid, :dash];
cols = [:blue, :blue, :red, :red];

for (i,sol) in enumerate(sols)
    lines!(ax, [s.up[end] for s in sol],
           [norm(s.up[1:2]) for s in sol],
           linewidth=2, linestyle=styls[i],
           color=cols[i],
           label=LaTeXString("$(typs[i]) \$\\mathcal{O}(\\varepsilon^$(ords[i]-1))\$"))
end
lines!(ax, Oms,  abs.(Aans), color=:black, linewidth=4,
       linestyle=:dot, label="Exact");

vlines!(ax, Oms_p, color=:black)
scatter!(ax, Oms_p, abs.([ansol((;pars...,Om=om)) for om in Oms_p]), color=:white,
         strokewidth=2, strokecolor=:black)
pos = (ax.finallimits.val.origin[2]./ax.finallimits.val.widths[2]*(ci==2 ? 100 : 60));
text!(ax, Oms_p, pos*ones(3),
      text=[LaTeXString("\$\\Omega_$p\$") for p in pts_p]);

Omcr = (sqrt(-(2sqrt(pars.zt^4-3pars.zt^2+1))-2pars.zt^2+3)*pars.w0)/sqrt(5);
scatter!(ax, Omcr, abs.(ansol((;pars...,Om=Omcr))), color=:red,
         markersize=20)
text!(ax, Omcr-0.5, abs.(ansol((;pars...,Om=Omcr)))*1.2, color=:red,
      text=L"$\Omega=\Omega_{cr}$")

if ci==1
    axislegend(ax);
end


if Makie.current_backend()==GLMakie
   display(scr, fig);
else
    display(fig)   
end
   
if savfigs
    save("./FIGS/A0_linearsteady_C$(ci).pdf", fig)
end

# * Nyquist Plot

xll = 0.35pars.w0..1.25pars.w0;

set_theme!(theme_latexfonts())
fsz = 18;
fig1 = Figure(fontsize=fsz, size=(500,600));
if !isdefined(Main, :scr1) && Makie.current_backend()==GLMakie
   scr1 = GLMakie.Screen();
end

ax = Axis(fig1[1, 1], xlabel=L"Coefficient $A_0/F$ (m/N)", ylabel=L"Coefficient $B_0/F$ (m/N)");
lines!(ax, real(Aans), -imag(Aans), color=:black, linewidth=2);
for (i,sol) in enumerate(sols)
    lines!(ax, [s.up[1] for s in sol], [s.up[2] for s in sol],
           color=cols[i], linestyle=styls[i])
end

Omcr = (sqrt(-(2sqrt(pars.zt^4-3pars.zt^2+1))-2pars.zt^2+3)*pars.w0)/sqrt(5);
scatter!(ax, real(ansol((;pars...,Om=Omcr))), -imag(ansol((;pars...,Om=Omcr))),
         color=:red)
algn = (ci==1) ? (:center, :bottom) : (:center, :top);
text!(ax, real(ansol((;pars...,Om=Omcr))), -imag(ansol((;pars...,Om=Omcr))),
      text=L"$\Omega=\Omega_{cr}$", color=:red, align=algn)

scatter!(ax, real(ansol((;pars...,Om=pars.w0))), -imag(ansol((;pars...,Om=pars.w0))),
         color=:black)
algn = (ci==1) ? (:center, :top) : (:left, :top);
text!(ax, real(ansol((;pars...,Om=pars.w0))), -imag(ansol((;pars...,Om=pars.w0))),
      text=L"$\Omega=\omega_0$", align=algn)

ax = Axis(fig1[2, 1], xlabel="Excitation Frequency (rad/s)", ylabel=L"Jacobian Entry $J_{2,1}$");
# lines!(ax, real(Aans), -imag(Aans), color=:black, linewidth=2);
for (i,sol) in enumerate(sols)
    lines!(ax, [s.up[3] for s in sol], [s.J[2,1] for s in sol],
           color=cols[i], linestyle=styls[i],
           label=LaTeXString("$(typs[i]) \$\\mathcal{O}(\\varepsilon^$(ords[i]-1))\$"))
end
axislegend(ax, nbanks=2)
scatter!(ax, pars.w0, 0, color=:black)
text!(ax, pars.w0, 0, text=L"$\Omega=\omega_0$", align=(:left, :top))

Omcr = (sqrt(-(2sqrt(pars.zt^4-3pars.zt^2+1))-2pars.zt^2+3)*pars.w0)/sqrt(5);
scatter!(ax, Omcr, 0, strokewidth=2, color=:red, strokecolor=:red);
text!(ax, Omcr, 0, strokewidth=2, color=:red, text=L"$\Omega=\Omega_{cr}$");
xlims!(ax, (xll.left, xll.right))
ylims!(ax, -3, 3)
# autolimits!(ax)

if Makie.current_backend()==GLMakie
   display(scr1, fig1);
else
   display(fig1)
end
   
if savfigs
   save("./FIGS/A0_Nyquist_C$(ci).pdf", fig1)
end
