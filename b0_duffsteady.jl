# * Preamble
using DifferentialEquations
using NonlinearSolve
using GLMakie
using ColorSchemes
using CairoMakie
using LaTeXStrings
using ForwardDiff
using JLD2
using LinearAlgebra

using Revise
using juliajim.HARMONIC
using juliajim.CONTINUATION
includet("./duffingsys.jl");

analyze = false;
savdats = true;

savfigs = true;
if savfigs
    CairoMakie.activate!();
else
    GLMakie.activate!();
end

cfgs = [(zt = 0.25e-2, w0 = 2.0, al = 0.01, f = [0.25, 0.5, 1.0], Om=1.9),
        (zt = 0.25e-1, w0 = 2.0, al = 0.2, f = [0.4, 0.8, 1.6], Om=1.9)];
ci = 2;
pars = cfgs[ci];

# * HB
h = 0:7;
Nhc = sum((h.==0) + 2(h.!=0));
_, _, _, rinds, iinds = HINDS(1, h);
Nt = 256;

u0 = zeros(Nhc);
u0[[rinds[1], iinds[1]]] .= 1.0;

Om0 = 0.01pars.w0;
Om1 = 2.5pars.w0;
dOm = 0.05;
cpars = (parm=:arclength, minDsc=1e-3, Dsc=:auto);

hbsols = [];

funduff(fi) = NonlinearFunction((du,u,p)->hbresfun!(du, u,
                                                    (;pars..., Om=p, f=pars.f[fi]),
                                                    h, Nt));
if analyze
    for (fi, f) in enumerate(pars.f)
        sol, _, _, _, _ = CONTINUATE(u0, funduff(fi), [Om0, Om1], dOm;
            cpars..., verbosity=0);
        push!(hbsols, sol)

        display("Done $fi/$(length(pars.f)).")
    end
end

# * Slow Flow Approaches
slowfun(fi,typ,ord) = NonlinearFunction((du,u,p)->sflow!(du,u,
    (;pars...,Om=p,f=pars.f[fi], typ=typ,order=ord)));
u0 = ones(2);
cpars = (parm=:arclength, Dsc=:none, nmax=1000, itopt=3);
dOm = 0.05pars.w0;

typs = [:MMS, :MMS, :EMS, :EMS];
ords = [1, 2, 1, 2];
slowsols = [[] for _ in 1:4]
slowh35s = [[] for _ in 1:4]
if analyze
    for (ti, (typ, ord)) in enumerate(zip(typs,ords))
        for (fi, f) in enumerate(pars.f)
            sol, _, _, _, _ = CONTINUATE(u0, slowfun(fi, typ, ord), [Om0, Om1], dOm;
                cpars..., verbosity=0);

            push!(slowsols[ti], sol)
            harms = (up->sflow_harms(up[1:2],
                (;pars...,Om=up[3],f=pars.f[fi],
                 typ=typ,order=ord))).(sol.up);
            push!(slowh35s[ti],
                ((h1,h2)->[real(h1), imag(h1), real(h2), imag(h2)]).(first.(harms),
                    last.(harms)))
                 
            display("($typ, $ord) Done $fi/$(length(pars.f)).")    
        end
    end
end

# * Save Data if Asked
if analyze && savdats
   @save "./DATS/B0_duffss.jld2" hbsols slowsols slowh35s
end
if !analyze
    @load "./DATS/B0_duffss.jld2" hbsols slowsols slowh35s
end

# * Plotting
set_theme!(theme_latexfonts())
fsz = 28;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Excitation Frequency (rad/s)",
          ylabel="Response (m)", yscale=Makie.pseudolog10);
for (fi, (sol,f)) in enumerate(zip(hbsols,pars.f))
    lines!(ax, last.(sol.up), norm.((u->u[1:end-1]).(sol.up)),
           label="F = $(pars.f) N", linewidth=4,
           color=colorschemes[:grays][fi/(length(pars.f)+1)])
end
cols = [:Blues, :Reds, :Greens, :budaS];
for (ti, (typ,ord,slows)) in enumerate(zip(typs,ords,slowsols))
    for (fi,(sol,f)) in enumerate(zip(slows,pars.f))
        lines!(ax, last.(sol.up), norm.((u->u[1:end-1]).(sol.up)),
               label=LaTeXString("$(typ) \$\\mathcal{O}(\\varepsilon^$(ord-1))\$, F = $(f) N"),
               color=colorschemes[cols[ti]][fi/(length(pars.f)+1)]);
    end
end
xlims!(ax, 1.5, 3.5);
axislegend(ax, position=(:left,:top), nbanks=2)

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
   
# * Separately
lstl = Linestyle([0,5,6]);

set_theme!(theme_latexfonts())
fsz = 22;
fig1 = Figure(fontsize=fsz, size=(1020, 820));
if !isdefined(Main, :scr1) && Makie.current_backend()==GLMakie
   scr1 = GLMakie.Screen();
end

axs = [];
for i in 1:2
    for j in 1:2
        local ax = Axis(fig1[i, j]);
        if j==1
            ax.ylabel=L"$H_1$ Response (m)"
        end
        if i==2
            ax.xlabel="Excitation Frequency (rad/s)";
        end
        push!(axs, ax)
    end
end
for (ti, ax) in enumerate(axs)
    for (fi, (sol,f)) in enumerate(zip(hbsols,pars.f))
        l1 = lines!(ax, last.(sol.up), norm.((u->u[[rinds[1], iinds[1]]]).(sol.up)),
            label="F = $f N", linewidth=5,
            color=colorschemes[:grays][fi/(length(pars.f)+1)]);
        if ti==1
            l1.label = "Harmonic Balance";
        end
    end
    
    for (fi, (sol,f)) in enumerate(zip(hbsols,pars.f))
        l2 = lines!(ax, last.(slowsols[ti][fi].up),
            norm.(((u)->u[1:2]).(slowsols[ti][fi].up)),
            label="F = $f N", linewidth=3,
            linestyle=lstl,
            color=colorschemes[:rainbow][1-fi/(length(pars.f)+1)]);
        if ti==1
            l2.label = "Asymptotic Solution";
        end
    end
    xlims!(ax, Om0, Om1)
    ylims!(ax, 0, 1.1round(maximum((u->norm(u[1:2])).(slowsols[1][3].up))))
    ax.yscale = Makie.pseudolog10;
    if typs[ti] == :MMS
        ax.title = LaTeXString("$(typs[ti]) \$\\mathcal{O}(\\varepsilon^$(ords[ti]-1))\$");
    else
        ax.title = LaTeXString("$(typs[ti]) \$\\mathcal{O}(p^$(ords[ti]-1))\$");
    end

    if ti ∈ [1, 3]
        axislegend(ax, merge=true, unique=false, position=:lt)
    end
end

if Makie.current_backend()==GLMakie
   display(scr1, fig1);
else
   display(fig1)
end
   
if savfigs
   save("./FIGS/B0_Duffresp_C$(ci).pdf", fig1)
end

# * Higher Harmonics
set_theme!(theme_latexfonts())
fsz = 22;
fig2 = Figure(fontsize=fsz, size=(1020,820));
if !isdefined(Main, :scr2) && Makie.current_backend()==GLMakie
   scr2 = GLMakie.Screen();
end

axs = [];
for i in 1:2
    for j in 1:2
        local ax = Axis(fig2[i, j]);
        if j==1
            ax.ylabel=L"$H_3$ Response (m)"
        end
        if i==2
            ax.xlabel="Excitation Frequency (rad/s)";
        end
        push!(axs, ax)
    end
end
for (ti, ax) in enumerate(axs)
    for (fi, (sol,f)) in enumerate(zip(hbsols,pars.f))
        l1 = lines!(ax, last.(sol.up), norm.((u->u[[rinds[3], iinds[3]]]).(sol.up)),
                    label="F = $f N", linewidth=5,
                    color=colorschemes[:grays][fi/(length(pars.f)+1)]);
        if ti==2
            l1.label = "Harmonic Balance";
        end
    end
    for (fi, (sol,f)) in enumerate(zip(hbsols,pars.f))
        l2 = lines!(ax, last.(slowsols[ti][fi].up),
                    norm.(((u)->u[1:2]).(slowh35s[ti][fi])),
                    label="F = $f N", linewidth=3,
                    linestyle=lstl,
                    color=colorschemes[:rainbow][1-fi/(length(pars.f)+1)]);
        if ti==2
            l2.label = "Asymptotic Solution";
        end
    end
    xlims!(ax, Om0, Om1)
    # ylims!(ax, 0, maximum((u->norm(u[1:2])).(slowh35s[1][3])))
    ylims!(ax, 0, 1.3maximum((u->norm(u[[rinds[3], iinds[3]]])).(hbsols[end].up)))
    ax.yscale = Makie.pseudolog10;
    if typs[ti] == :MMS
        ax.title = LaTeXString("$(typs[ti]) \$\\mathcal{O}(\\varepsilon^$(ords[ti]-1))\$: \$\\mathcal{O}(\\varepsilon^$(ords[ti]))\$ Inhomogeneous");
    else
        ax.title = LaTeXString("$(typs[ti]) \$\\mathcal{O}(p^$(ords[ti]-1))\$: \$\\mathcal{O}(p^$(ords[ti]))\$ Inhomogeneous");
    end

    if ti ∈ [2, 4]
        axislegend(ax, merge=true, unique=false, position=:rt)
    end
end

if Makie.current_backend()==GLMakie
   display(scr2, fig2);
else
   display(fig2)
end
   
if savfigs
    save("./FIGS/B0_Duffresp_C$(ci)_H3.pdf", fig2)
end
