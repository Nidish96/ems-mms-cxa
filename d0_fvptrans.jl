# * Preamble
using DifferentialEquations
using NonlinearSolve
using LaTeXStrings
using LinearAlgebra
using JLD2: @save, @load

using Revise

using juliajim.HARMONIC
using juliajim.CONTINUATION
includet("./vdpsys.jl")

savfigs = false;
if savfigs
    using CairoMakie
    CairoMakie.activate!();
else
    using GLMakie
    GLMakie.activate!();
end

savdats = true;
analyze = true;

pars = (c=0.04, w0=2., al=0.1, f=1., Om=0.1);

# * Transient Simulation
Nom = 20;
Oms = range(0.1, 4, Nom);
ncyc = 80;
nt = 1024;

iw = 1;
prob = ODEProblem(rocfun!, zeros(2), [0, 2π/Oms[iw]*ncyc], (;pars..., Om=Oms[iw]));
sol = solve(prob);

tl = range(0, 2π/Oms[iw]*ncyc, ncyc*nt÷2);
soll = sol(tl);

set_theme!(theme_latexfonts())
fsz = 18;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="A0");
lines!(ax, soll.t, getindex.(soll.u, 1))
# lines!(ax, soll)
# lines!(ax, getindex.(soll.u, 1), getindex.(soll.u, 2))

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end

# * Forced Response
if analyze
    Nom = 160;
    Oms = range(0.1, 4, Nom);
    ncyc = 480;
    nt = 512;

    Umaxs = zeros(Nom);
    for (iw,Om) in enumerate(Oms)
        lprob = ODEProblem(rocfun!, zeros(2), [0, 2π/Om*ncyc], (;pars..., Om=Om));
        lsol = solve(lprob);

        ltl = range(2π/Om*ncyc÷2, 2π/Om*ncyc, ncyc*nt÷2);
        lsoll = lsol(ltl);

        Umaxs[iw] = maximum(abs.(getindex.(lsoll.u,1)));
    end
    if savdats
        @save "./DATS/D0_trdat.jld2" Oms Umaxs
    end
else
    @load "./DATS/D0_trdat.jld2" Oms Umaxs
end

# Plotting
set_theme!(theme_latexfonts())
fsz = 18;
fig1 = Figure(fontsize=fsz);
if !isdefined(Main, :scr1) && Makie.current_backend()==GLMakie
   scr1 = GLMakie.Screen();
end

ax = Axis(fig1[1, 1], xlabel=L"Excitation Frequency $\Omega$",
    ylabel="Peak Response (m)");
scatterlines!(ax, Oms, Umaxs)

if Makie.current_backend()==GLMakie
   display(scr1, fig1);
else
   display(fig1)
end
   


