# * 
using DifferentialEquations
using NonlinearSolve
using GLMakie
using CairoMakie
using LaTeXStrings
using ForwardDiff

using Revise

using juliajim.HARMONIC
using juliajim.CONTINUATION

includet("./linearsys.jl");

savfigs = true;
if savfigs
    CairoMakie.activate!();
else
    GLMakie.activate!();
end

cfgs = [(zt = 0.25e-2, w0 = 2.0, f = 1.0),
        (zt = 25e-2, w0 = 2.0, f = 1.0)];
ci = 1;

# * Compute
pars = cfgs[ci];

typs = [:MMS, :MMS, :EMS, :EMS];
ords = [1, 2, 1, 2];
sols = [];

pts = ["A", "B", "C"];
Oms = [0.95pars.w0, 1.5pars.w0, 0.1pars.w0];
ncyc = (ci==1) ? 500 : 20;
tpl = range(0, 2π*ncyc, 1024ncyc);

iw = 2;

lprob = ODEProblem(rocfun!, zeros(2), [0, 2π/Oms[iw]*ncyc], (;pars..., Om=Oms[iw]));
lsol = solve(lprob);
lsol = lsol(tpl/Oms[iw]);

# * Numerical Slow Flow Solution
# for (ti, typ) in enumerate(typs)
#     sprob = remake(lprob, f=sflow!, p=(;pars..., Om=Oms[iw], typ=typ, order=ords[ti]));
#     sol = solve(sprob, RK4());
#     push!(sols, sol(tpl/Oms[iw]));

#     display("Done $(typ) O($(ords[ti])).")    
# end

# * Analytical Slow Flow Solution
Js = [ForwardDiff.jacobian((du,u)->sflow!(du, u, (;pars...,typ=typ,order=ord,Om=Oms[iw])),
                           ones(2), zeros(2)) for (typ,ord) in zip(typs,ords)];
Fs = [sflow!(zeros(2), zeros(2), (;pars...,typ=typ,order=ord,Om=Oms[iw]))
      for (typ,ord) in zip(typs,ords)];

sols = [(t=tpl/Oms[iw], u=[J\ (exp(J*t)-I(2))*F for t in tpl/Oms[iw]]) for (J,F) in zip(Js,Fs)];

if ci==1
    ncyc2 = floor(Int, ncyc/10);
    tpl2 = range(0.0, 2π*ncyc2, 16192ncyc2);
    sols[4] = (t=tpl2/Oms[iw], u=[Js[4]\ (exp(Js[4]*t)-I(2))*Fs[4]
                                 for t in tpl2/Oms[iw]]);
end

# Reconstruction
ureco = [[u[1]cos(Oms[iw]*t)+u[2]sin(Oms[iw]*t) for (u,t) in zip(sol.u,sol.t)] for sol in sols];

# * Plot

xl1 = 2π/Oms[iw]*[10, 200];
xl2 = 2π/Oms[iw]*(ncyc .- [100, 0])
ylmin = 0;
ylmax = max(maximum(abs.(getindex.(lsol.u,1))),
            maximum([maximum(norm.(s.u)) for s in sols[[1,3]]]));    
if ci==1
    if iw==2
        xl1 = [0, 2π/Oms[iw]*20];
        ylmin = -0.2ylmax;
    elseif iw==3
        xl1 = [0, 2π/Oms[iw]*2];
        xl2 = 2π/Oms[iw]*(ncyc .- [20, 0])
        ylmin = -1.1ylmax;
    end
else
    xl1 = 2π/Oms[iw]*[0, 5];
    if iw==2
        xl1 = [0, 2π/Oms[iw]*5];
        ylmin = -0.2ylmax;
    elseif iw==3
        xl1 = [0, 2π/Oms[iw]*2];
        xl2 = 2π/Oms[iw]*(ncyc .- [20, 0])
        ylmin = -1.1ylmax;
    end
end

set_theme!(theme_latexfonts())
fsz = 18;
if (ci==1 && iw==1) || (ci==2 && iw==1)
    fsz = 27;
end
fig = Figure(fontsize=fsz, size=(800, 600));
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end
cols = [:blue, :blue, :red, :red];
styls = [:solid, :dash, :solid, :dot];
lws = [2, 3, 2, 4];

if iw==3 || iw==2
    ax = Axis(fig[1, 1:4], ylabel="Response (m/N)");
else
    ax = Axis(fig[1, 3:8], ylabel="Response (m/N)");
end
plotall!(axs) = begin 
    lines!(axs, lsol.t, getindex.(lsol.u, 1), color=:black,
           linewidth=0.25, label="Reference");
    for (ti, sol) in enumerate(sols)
        lines!(axs, sol.t, norm.(sol.u), linewidth=lws[ti],
               color=cols[ti], linestyle=styls[ti],
               label=LaTeXString("$(typs[ti]) \$\\mathcal{O}(\\varepsilon^$(ords[ti]-1))\$"));
    end
end
box1!(axs, w=2, ylm=ylmax) = poly!(axs, Point2f.([(xl1[1],ylmin), (xl1[2],ylmin),
                                                  (xl1[2],1.1ylm), (xl1[1],1.1ylm)]),
    color=(:red,0), strokecolor=:green1, strokewidth=w);
box2!(axs, w=2) = poly!(axs, Point2f.([(xl2[1],ylmin), (xl2[2],ylmin),
                                       (xl2[2],1.1ylmax), (xl2[1],1.1ylmax)]),
    color=(:red,0), strokecolor=:red, strokewidth=w);

plotall!(ax)
box1!(ax)
# if iw!=3 && iw!=2
#     box2!(ax)
# end

xlims!(ax, 0, 2π/Oms[iw]*ncyc);
ylims!(ax, -1.1ylmax, 1.1ylmax);

if iw==1 # show frequency content of CXA
    Legend(fig[1,1:2], ax);
# elseif iw==2
#     ax = Axis(fig[1,1:2], xlabel="Excitation Frequency (rad/s)", ylabel="Response (m/N)",
#               yscale=log10)
#     Nom = 100;
#     Oms_pl = range(0.05pars.w0, 2pars.w0, Nom);
#     Aan = [ansol((;pars...,Om=om)) for om in Oms_pl];
#     Aans = [ansol((;pars...,Om=om)) for om in Oms];

#     lines!(ax, Oms_pl, abs.(Aan), color=:black)
#     scatter!(ax, Oms, abs.(Aans), color=:white, strokecolor=:black, strokewidth=2)
#     text!(ax, Oms.+[-0.5,0,0], abs.(Aans), text=pts)    
end

if iw==1
    ax = Axis(fig[2, 1:8], xlabel="Time (s)", ylabel="Response (m/N)");
else
    ax = Axis(fig[2, 1:4], xlabel="Time (s)", ylabel="Response (m/N)");
end
plotall!(ax);
xlims!(ax, xl1[1], xl1[2]);

box1!(ax, 6)
ylims!(ax, ylmin, 1.1ylmax)

if iw==3 || iw==2
    ax = Axis(fig[1, 5:8]);

    cols_e = [:blue, :green1];
    for (ur,ord,col) in zip(ureco[1:2],ords[1:2],cols_e)
        lines!(ax, lsol.t, ur, color=col, linewidth=2,
               label=LaTeXString("MMS \$\\mathcal{O}(\\varepsilon^{$(ord-1)})\$"))
    end
    lines!(ax, lsol.t, getindex.(lsol.u, 1), color=:black, linewidth=1,
           label="Reference");    
    axislegend(ax, nbanks=2)
    xlims!(ax, xl1[1], xl1[2]);
    ylims!(ax, -ylmax, 1.75ylmax)

    ax = Axis(fig[2, 5:8], xlabel="Time (s)");

    cols_e = [:red, :gold]
    for (sol,ur,ord,col) in zip(sols[3:4],ureco[3:4],ords[3:4],cols_e)
        lines!(ax, sol.t, ur, color=col, linewidth=2,
               label=LaTeXString("EMS \$\\mathcal{O}(\\varepsilon^{$(ord-1)})\$"))
    end
    lines!(ax, lsol.t, getindex.(lsol.u, 1), color=:black, linewidth=1,
           label=L"2$^{nd}$ Order ODE");    
    axislegend(ax, nbanks=2)
    xlims!(ax, xl1[1], xl1[2]);
    ylims!(ax, -ylmax, 1.75ylmax)

    if iw==3
        colp = :magenta;
        poly!(ax, [(xl1[1], -0.2ylmax), (xl1[2]/150, -0.2ylmax), (xl1[2]/150, 1.1ylmax), (xl1[1], 1.1ylmax)], color=(:red,0), strokewidth=2, strokecolor=colp);

        axin = Axis(fig[2, 5:8],
                    width=Relative(0.4),
                    height=Relative(0.3),
                    halign=0.97,
                    valign=0.125,
                    backgroundcolor=:white,
                    xticklabelsize=0.7fsz,
                    yticklabelsize=0.7fsz);
        for (sol,ur,ord,col) in zip(sols[3:4],ureco[3:4],ords[3:4],cols_e)
            lines!(axin, sol.t, ur, color=col)
        end
        lines!(axin, lsol.t, getindex.(lsol.u, 1), color=:black, linewidth=0.5);
        xlims!(axin, xl1[1], xl1[2]/150);
        ylims!(axin, -0.2ylmax, 1.1ylmax);
        translate!(axin.blockscene, 0, 0, 1000);
        axin.yticks = LinearTicks(3);
        poly!(axin, [(xl1[1], -0.2ylmax), (xl1[2]/150, -0.2ylmax), (xl1[2]/150, 1.1ylmax), (xl1[1], 1.1ylmax)], color=(:red,0), strokewidth=4, strokecolor=colp);
    end
# else
#     ax = Axis(fig[2, 5:8], xlabel="Time (s)");
    
#     plotall!(ax)
#     box2!(ax, 6)
#     xlims!(ax, xl2[1], xl2[2]);
#     ylims!(ax, ylmin, 1.1ylmax)
end

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
    display(fig)   
end
   
if savfigs
    save("./FIGS/A1_linearSS_C$(ci)_W$(iw).pdf", fig)
end
