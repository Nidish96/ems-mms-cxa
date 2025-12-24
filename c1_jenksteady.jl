# * 
using DifferentialEquations
using NonlinearSolve
using LinearAlgebra
using GLMakie
using ColorSchemes
using CairoMakie
using LaTeXStrings
using ForwardDiff
using JLD2: @save, @load

using Revise

using juliajim.HARMONIC
using juliajim.CONTINUATION

includet("./jenkinsys_mod.jl");

savfigs = false;
if savfigs
    CairoMakie.activate!();
else
    GLMakie.activate!();
end

cfgs = [(zt = 0.25e-2, w0 = 2.0, kt=5, fs=1.0, f = [0.1, 0.25, 0.5, 1.0, 1.25], Om=1.9),
        (zt = 0.25e-1, w0 = 2.0, kt=5, fs=1.0, f = [0.4, 0.8, 1.6], Om=1.9)];
ci = 1;
pars = cfgs[ci];

analyze = true;
savdats = true;

# * HB
h = 0:7;
Nhc = sum((h.==0) + 2(h.!=0));
_, _, _, rinds, iinds = HINDS(1, h);
Nt = 512;

u0 = zeros(Nhc);
u0[[rinds[1], iinds[1]]] .= 1e-2;

funjenk(fi) = NonlinearFunction((du,u,p)->hbresfun!(du,u,
                                                    (;pars...,f=pars.f[fi],Om=p),h,Nt));

Om0 = 0.05pars.w0;
Om1 = 2.5pars.w0;
dOm = 0.05pars.w0;
cpars = (parm=:arclength, nmax=1000);

if analyze
    hbsols = Vector{Vector{myNLSoln}}();

    for fi in eachindex(pars.f)

        itopt = (fi==1) ? 3 : :auto;
        
        sol, _, _, _, _ = CONTINUATE(u0, funjenk(fi), [Om0, Om1], dOm; cpars...,
                                     verbosity=0);


        push!(hbsols, sol);
        display("Done $(fi)/$(length(pars.f)).");
    end

    if savdats
        @save "./DATS/C1_hbsols.jld2" hbsols
    end
else
    @load "./DATS/C1_hbsols.jld2" hbsols
    # hbsols = Vector{myNLSoln}.(hbsols);
end

# * Asymptotic Methods
ords = [1, 2, 1, 2, 1, 2];
typs = [:MMS, :MMS, :EMS, :EMS, :EMS0, :EMS0];

fun(ti, fi) = NonlinearFunction((du,u,p)->sflow!(du, u, (;pars...,typ=typs[ti],
                                                         order=ords[ti],
                                                         f=pars.f[fi],Om=p)));
u0 = [1e-2; 1];

dOm = 0.05pars.w0;  # 0.05pars.w0
cpars = (parm=:arclength, Dsc=:none, nmax=4000);

if analyze
    slowsols = [Vector{myNLSoln}[] for _ in typs];
    slowh35s = [Vector{Vector{Float64}}[] for _ in typs];
    for (ti, (typ,ord)) in enumerate(zip(typs,ords))
        for (fi, f) in enumerate(pars.f)
            sol, _, _, _, _ = CONTINUATE(u0, fun(ti, fi), [Om1, Om0], dOm;
                cpars..., verbosity=0);
            sol = sol[end:-1:1];

            harms = (up->sflow_harms(up[1:2],
                (;pars...,Om=up[3],f=f,typ=typ,order=ord))).(sol.up);
            push!(slowh35s[ti],
                ((h1,h2)->[real(h1), imag(h1), real(h2), imag(h2)]).(first.(harms),
                    last.(harms)))

            push!(slowsols[ti], sol);
            display("$(typ), $(ord): Done $(fi)/$(length(pars.f)) " *
                    " with $(length(sol)) points.")
        end
    end
    if savdats
        @save "./DATS/C1_asymsols.jld2" slowsols slowh35s
    end
else
    @load "./DATS/C1_asymsols.jld2" slowsols slowh35s
#    slowsols = [(u->Vector{myNLSoln}(u)).(s) for s in slowsols];
end

# * Plot
lstl = Linestyle([0, 5,6]);

set_theme!(theme_latexfonts())
fsz = 21;
fig = Figure(fontsize=fsz, size=(1020,410));
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

plothb!(ax, sols, cschem, wd=3, lab="") = begin
    for (fi, (sol, f)) in enumerate(zip(sols, pars.f))
        lines!(ax, last.(sol.up), (u->norm(u[[rinds[1], iinds[1]]])/f).(sol.up),
            linewidth=wd, label=(lab=="") ? "F = $(f) N" : lab,
            color=colorschemes[cschem][fi/(length(pars.f)+1+1)])
    end
end;
plotsol!(ax, sols, cschem, wd=3, lab="", lst=:solid) = begin
    for (fi, (sol, f)) in enumerate(zip(sols, pars.f))
        lines!(ax, last.(sol.up), (u->norm(u[1:end-1])/f).(sol.up),
            linewidth=wd, label=(lab=="") ? "F = $(f) N" : lab,
            color=colorschemes[cschem][1-fi/(length(pars.f)+1)],
            linestyle=lst)
    end
end;

plis = [1,2];
axs = [];
for (ti, (typ,ord,slowsol)) in enumerate(zip(typs[plis],ords[plis],slowsols[plis]))
    if typ == :MMS
        ttl = LaTeXString("$(typ) \$\\mathcal{O}(\\varepsilon^$(ord-1))\$");
    else
        ttl = LaTeXString("$(typ) \$\\mathcal{O}(p^$(ord-1))\$");
    end

    ax = Axis(fig[1, ti], xlabel="Excitation Frequency (rad/s)",
              title=ttl);
    if ti==1
        ax.ylabel=L"$H_1$ Response (m/N)"
        lab1 = "";
        lab2 = "";
    else
        lab1 = "Harmonic Balance";
        lab2 = "Asymptotic Solution";
    end
    # hlines!(ax, pars.fs/pars.kt, color=:black)
    plothb!(ax, hbsols, :grays, 5, lab1)
    plotsol!(ax, slowsol, :rainbow, 3, lab2, lstl)

    if ti==1
        axislegend(ax, merge=true, unique=false, position=:lt)
    elseif ti==2
        axislegend(ax, merge=true, unique=false, position=:lt, nbanks=2)
    end
    xlims!(ax, Om0, Om1)
    ylims!(ax, 0, 2.85)

    push!(axs, ax);
end
linkaxes!(axs...)
              
if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end

if savfigs
   save("./FIGS/C1_H1resps_C$(ci).pdf", fig)
end
   
# * Third Harmonic

set_theme!(theme_latexfonts())
fsz = 22;
fig1 = Figure(fontsize=fsz, size=(1020, 410));
if !isdefined(Main, :scr1) && Makie.current_backend()==GLMakie
   scr1 = GLMakie.Screen();
end

plothb3!(ax, sols, cschem, wd=3, lab="") = begin
    for (fi, (sol, f)) in enumerate(zip(sols, pars.f))
        lines!(ax, last.(sol.up), (u->norm(u[[rinds[3], iinds[3]]])).(sol.up),
            linewidth=wd,
            label=(lab=="") ? "F = $(f) N" : lab,
            color=colorschemes[cschem][fi/(length(pars.f)+1+1)])
    end
end;
plotsol3!(ax, sols, h35s, cschem, wd=3, lab="", lst=:solid) = begin
    for (fi, (sol, h35, f)) in enumerate(zip(sols, h35s, pars.f))
        lines!(ax, last.(sol.up), (u->norm(u)).(h35),
            linewidth=wd,
            label=(lab=="") ? "F = $(f) N" : lab,
            color=colorschemes[cschem][1-fi/(length(pars.f)+1)],
            linestyle=lst)
    end
end;

plis = 3:4;
axs = [];
for (ti, (typ, ord, slowsol, slowh35)) in enumerate(zip(typs[plis],ords[plis], slowsols[plis], slowh35s[plis]))
    if typ == :MMS
        ttl = LaTeXString("$(typ) \$\\mathcal{O}(\\varepsilon^$(ord-1)):\\, " *
                          "\\mathcal{O}(\\varepsilon^$(ord))\$ Inhomogeneous");
    else
        ttl = LaTeXString("$(typ) \$\\mathcal{O}(p^$(ord-1)):\\, " *
                          "\\mathcal{O}(p^$(ord))\$ Inhomogeneous");
    end

    ax = Axis(fig1[1, ti], xlabel="Excitation Frequency (rad/s)",
        title=ttl);
    if ti==1
        ax.ylabel=L"$H_3$ Response (m/N)"
        lab1 = "";
        lab2 = "";
    else
        lab1 = "Harmonic Balance";
        lab2 = "Asymptotic Solution";
    end
    # hlines!(ax, pars.fs/pars.kt, color=:black)
    plothb3!(ax, hbsols, :grays, 5, lab1)
    plotsol3!(ax, slowsol, slowh35, :rainbow, 3, lab2, lstl)

    if ti<=2
        axislegend(ax, merge=true, unique=false, position=:rt)
    end
    xlims!(ax, Om0, Om1)

    push!(axs, ax);
end
linkaxes!(axs...);

if Makie.current_backend()==GLMakie
   display(scr1, fig1);
else
   display(fig1)
end
   
if savfigs
   save("./FIGS/C1_H3resps_C$(ci).pdf", fig1)
end

# * Plot the Re-Averaged EMS

pli = 6;

set_theme!(theme_latexfonts())
fsz = 22;
fig2 = Figure(fontsize=fsz, size=(1020, 410));
if !isdefined(Main, :scr2) && Makie.current_backend()==GLMakie
   scr2 = GLMakie.Screen();
end

ax1 = Axis(fig2[1, 1], xlabel="Excitation Frequency (rad/s)",
    ylabel=L"$H_1$ Response (m/N)");
lab1 = "";
lab2 = "";

plothb!(ax1, hbsols, :grays, 5, lab1)
plotsol!(ax1, slowsols[pli], :rainbow, 3, lab2, lstl)

axislegend(ax1, merge=true, unique=false, position=:lt, nbanks=1)

ax2 = Axis(fig2[1, 2], xlabel="Excitation Frequency (rad/s)",
    ylabel=L"$H_3$ Response (m)");
lab1 = "Harmonic Balance";
lab2 = "Asymptotic Solution";
plothb3!(ax2, hbsols, :grays, 5, lab1)
plotsol3!(ax2, slowsols[pli], slowh35s[pli], :rainbow, 3, lab2, lstl)

linkxaxes!(ax1, ax2);
xlims!(ax1, Om0, Om1);
linkyaxes!(ax2, axs[1]);

axislegend(ax2, merge=true, unique=false, position=:rt, nbanks=1)

if Makie.current_backend()==GLMakie
   display(scr2, fig2);
else
   display(fig2)
end

if savfigs
    save("./FIGS/C1_waEMS_C$(ci)_pli$(pli).pdf", fig2)
end
