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
pars = cfgs[ci];

typs = [:MMS, :MMS, :EMS, :EMS];
ords = [1, 2, 1, 2];

# * Compute Forced Response Through Continuation

u0 = [0.0,0.0];
Om0 = 0.05pars.w0;
Om1 = 2pars.w0;
dOm = 0.05pars.w0;
cpars = (save_jacs=true, parm=:arclength, minDsc=1e-1, nmax=10000);

fun = NonlinearFunction((r,u,p)->sflow!(r, u, (;pars...,
                                               typ=typs[3],order=1,Om=p)));
sol = first(CONTINUATE(u0, fun, [Om0, Om1], dOm; cpars...));

# * Transient analysis
ncyc = (ci==1) ? 500 : 20;
ncyc = 10;
tpl = range(0, 2π*ncyc, 1024ncyc);

lprob(om) = ODEProblem(rocfun!, zeros(2), [0, 2π/om*ncyc], (;pars..., Om=om));

corrs = [];

Om = sol.up[1][3];

for (wi, Om) in enumerate(last.(sol.up))
    lsol = solve(lprob(Om), RK4());
    Aan = ansol((;pars..., Om=Om));

    lsol = lsol(tpl/Om);
    utr = first.(lsol.u)-real(Aan*exp.(im*Om*lsol.t));

    # Slow flow matrices
    Js = [ForwardDiff.jacobian((du,u)->sflow!(du, u, (;pars...,typ=typ,order=ord,Om=Om)),
                               ones(2), zeros(2)) for (typ,ord) in zip(typs,ords)];
    Fs = [sflow!(zeros(2), zeros(2), (;pars...,typ=typ,order=ord,Om=Om))
          for (typ,ord) in zip(typs,ords)];

    sols = [[J\ (exp(J*t)-I(2))*F for t in tpl/Om] for (J,F) in zip(Js,Fs)];
    solss = [-J\F for (J,F) in zip(Js,Fs)];

    ureco = [[u[1]cos(Om*t)+u[2]sin(Om*t) for (u,t) in zip(sol,tpl/Om)] for sol in sols];
    utrr = [u-(uss[1]cos.(tpl)+uss[2]sin.(tpl)) for (u,uss) in zip(ureco,solss)];

    push!(corrs, [utr'uu/norm(uu)/norm(utr) for uu in utrr])

    display("Done $(wi)/$(length(sol.up))")
end

# * Plot
set_theme!(theme_latexfonts())
fsz = 22;
lwd = 2.5;
fig = Figure(fontsize=fsz, size=(610, 610));
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], ylabel="Response (m/N)", yscale=log10);
lines!(ax, last.(sol.up), (u->norm(u[1:2])).(sol.up), color=:black, linewidth=lwd)

ax = Axis(fig[2, 1], ylabel="MMS Correlation");
lines!(ax, last.(sol.up), getindex.(corrs, 1), color=:blue,
       label=L"$\mathcal{O}(\varepsilon^0)$", linewidth=lwd)
lines!(ax, last.(sol.up), getindex.(corrs, 2), color=:red,
       label=L"$\mathcal{O}(\varepsilon^1)$", linewidth=lwd)
axislegend(ax, position=:rb)
ylims!(ax, -0.1, 1.1)

ax = Axis(fig[3, 1], xlabel="Excitation Frequency (rad/s)", ylabel="EMS Correlation");
lines!(ax, last.(sol.up), abs.(getindex.(corrs, 3)), color=:blue, linewidth=lwd)
lines!(ax, last.(sol.up), abs.(getindex.(corrs, 4)), color=:red, linewidth=lwd)
ylims!(ax, -0.1, 1.1)
      
if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
   
if savfigs
   save("./FIGS/A3_tcorr_C$(ci).pdf", fig)
end
