# * 
using DifferentialEquations
using NonlinearSolve
using GLMakie
using LaTeXStrings
using ForwardDiff

using Revise
includet("/home/nidish/Documents/Research/b_Programming/juliajim/src/HARMONIC.jl")
includet("/home/nidish/Documents/Research/b_Programming/juliajim/src/CONTINUATION.jl")
includet("./linearsys.jl");

# * Second Order System
pars = (zt = 0.1e-2, w0 = 2.0, f = 1.0, Om=1.9);

ncyc = 500;
Tmax = 2π/pars.Om*ncyc;
lprob = ODEProblem(rocfun!, [0.0, 0.0], [0.0, Tmax], pars);
sol = solve(lprob, RK4());

Aan = ansol(pars);
uan = abs(Aan)cos.(pars.Om*sol.t .+ angle(Aan));
utr = [u[1] for u in sol.u]-uan;

# ** MMS
mmspars = (;pars..., typ=:MMS, order=1);
sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,Tmax], mmspars);
sol_mms1 = solve(sprob, RK4());

sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,Tmax], (;mmspars...,order=2));
sol_mms2 = solve(sprob, RK4());

# *** Steady State
ssprob = SteadyStateProblem(sflow!, [0.0,0.0], mmspars);
ss_mms1 = solve(ssprob);

ssprob = SteadyStateProblem(sflow!, [0.0,0.0], (;mmspars...,order=2));
ss_mms2 = solve(ssprob);

# ** EMS
emspars = (;pars..., typ=:EMS, order=1);
sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,Tmax], emspars);
sol_ems1 = solve(sprob, RK4());

sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,Tmax], (;emspars...,order=2));
sol_ems2 = solve(sprob, RK4());

sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,Tmax], (;emspars...,order=3));
sol_ems3 = solve(sprob, RK4());

# *** Steady State
ssprob = SteadyStateProblem(sflow!, [0.0,0.0], emspars);
ss_ems1 = solve(ssprob);

ssprob = SteadyStateProblem(sflow!, [0.0,0.0], (;emspars...,order=2));
ss_ems2 = solve(ssprob);

ssprob = SteadyStateProblem(sflow!, [0.0,0.0], (;emspars...,order=3));
ss_ems3 = solve(ssprob);

# * Summaries, Interpolated on ODE Grid
u_mms1 = [u[1]cos(pars.Om*sol.t[i])+u[2]sin(pars.Om*sol.t[i]) for (i, u) in enumerate(sol_mms1(sol.t))];
u_mms2 = [u[1]cos(pars.Om*sol.t[i])+u[2]sin(pars.Om*sol.t[i]) for (i, u) in enumerate(sol_mms2(sol.t))];

utr_mms1 = u_mms1 - (ss_mms1[1]cos.(pars.Om*sol.t)+ss_mms1[2]sin.(pars.Om*sol.t));
utr_mms2 = u_mms2 - (ss_mms2[1]cos.(pars.Om*sol.t)+ss_mms2[2]sin.(pars.Om*sol.t));

u_ems1 = [u[1]cos(pars.Om*sol.t[i])+u[2]sin(pars.Om*sol.t[i]) for (i, u) in enumerate(sol_ems1(sol.t))];
u_ems2 = [u[1]cos(pars.Om*sol.t[i])+u[2]sin(pars.Om*sol.t[i]) for (i, u) in enumerate(sol_ems2(sol.t))];
u_ems3 = [u[1]cos(pars.Om*sol.t[i])+u[2]sin(pars.Om*sol.t[i]) for (i, u) in enumerate(sol_ems3(sol.t))];

utr_ems1 = u_ems1 - (ss_ems1[1]cos.(pars.Om*sol.t)+ss_ems1[2]sin.(pars.Om*sol.t));
utr_ems2 = u_ems2 - (ss_ems2[1]cos.(pars.Om*sol.t)+ss_ems2[2]sin.(pars.Om*sol.t));
utr_ems3 = u_ems3 - (ss_ems3[1]cos.(pars.Om*sol.t)+ss_ems3[2]sin.(pars.Om*sol.t));

# * Plotting

set_theme!(theme_latexfonts())
fsz = 24;
fig = Figure(fontsize=fsz);
if !(@isdefined scr) && (@isdefined GLMakie)
    scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1:2], xlabel="Time (s)", ylabel="Response", title="Transient Response");

lines!(ax, sol.t, [u[1] for u in sol.u], label="ODE Solution")

lines!(ax, sol_mms1.t, norm.(sol_mms1.u),
       label=L"$\mathcal{O}(1)$ MMS")
lines!(ax, sol_mms2.t, norm.(sol_mms2.u),
       label=L"$\mathcal{O}(\epsilon)$ MMS")

lines!(ax, sol_ems1.t, norm.(sol_ems1.u), 
       label=L"$\mathcal{O}(1)$ EMS")
lines!(ax, sol_ems2.t, norm.(sol_ems2.u), 
       label=L"$\mathcal{O}(\epsilon)$ EMS")
lines!(ax, sol_ems3.t, norm.(sol_ems3.u), 
       label=L"$\mathcal{O}(\epsilon^2)$ EMS")

hlines!(ax, abs(Aan), label="Analytical Steady State")
hlines!(ax, norm(ss_mms1), label=L"$\mathcal{O}(1)$ MMS-SS")
hlines!(ax, norm(ss_mms2), label=L"$\mathcal{O}(\epsilon)$ MMS-SS")
hlines!(ax, norm(ss_ems1), label=L"$\mathcal{O}(1)$ EMS-SS")
hlines!(ax, norm(ss_ems2), label=L"$\mathcal{O}(\epsilon)$ EMS-SS")
hlines!(ax, norm(ss_ems3), label=L"$\mathcal{O}(\epsilon^2)$ EMS-SS")

axislegend(ax, "Curves Key", nbanks=2)

ax = Axis(fig[2, 1], title="Transients Amplitudes");

lines!(ax, sol.t, utr, label="ODE Solution")

lines!(ax, sol.t, [norm(u-ss_mms1.u) for u in sol_mms1(sol.t)])
lines!(ax, sol.t, [norm(u-ss_mms2.u) for u in sol_mms2(sol.t)])

lines!(ax, sol.t, [norm(u-ss_ems1.u) for u in sol_ems1(sol.t)])
lines!(ax, sol.t, [norm(u-ss_ems2.u) for u in sol_ems2(sol.t)])
lines!(ax, sol.t, [norm(u-ss_ems3.u) for u in sol_ems3(sol.t)])

ax = Axis(fig[2, 2], title="Errors");
lines!(ax, sol.t, [u[1] for u in sol.u], label="ODE Solution")
lines!(ax, sol.t, [u[1] for u in sol.u]-u_mms1, label=L"$\mathcal{O}(1)$ MMS Error")
lines!(ax, sol.t, [u[1] for u in sol.u]-u_mms2, label=L"$\mathcal{O}(\epsilon)$ MMS Error")

lines!(ax, sol.t, [u[1] for u in sol.u]-u_ems1, label=L"$\mathcal{O}(1)$ EMS Error")
lines!(ax, sol.t, [u[1] for u in sol.u]-u_ems2, label=L"$\mathcal{O}(\epsilon)$ EMS Error")
lines!(ax, sol.t, [u[1] for u in sol.u]-u_ems3, label=L"$\mathcal{O}(\epsilon^2)$ EMS Error")

axislegend(ax)

if isdefined(Main, :GLMakie)
    display(scr, fig);
else
    fig   
end
