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

# * Plot
pars = (zt = 0.1e-2, w0 = 2.0, f = 1.0, Om=0.1);

t_an = range(0, 2π/pars.Om, 128);
Aan = ansol(pars);
uan = real(Aan*exp.(im*pars.Om*t_an));
udan = real(im*pars.Om*Aan*exp.(im*pars.Om*t_an));

xlims = -0.5..2.5;
xdlims = -0.5..0.5;

Alims = -1.5abs(Aan)..1.5abs(Aan);
Blims = -1.5abs(Aan)..1.5abs(Aan);

du = zeros(2);

set_theme!(theme_latexfonts())
fsz = 24;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie)
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1:3], xlabel=L"State $A_0$", ylabel=L"State $B_0$",
          title=L"$\mathcal{O}(1)$ MMS");
streamplot!(ax, x->Point2(sflow!(du, x, (;pars...,typ=:MMS,order=1), 0.)),
            Alims, Blims)
scatter!(ax, real(Aan), -imag(Aan), marker=:xcross, color=:red, markersize=20);
scatter!(ax, Point2(solve(SteadyStateProblem(sflow!, [0.,0.],
                                             (;pars...,typ=:MMS,order=1)))),
         color=:blue, markersize=16);

ax = Axis(fig[1, 4:6], xlabel=L"State $A_0$", ylabel=L"State $B_0$",
          title=L"$\mathcal{O}(\epsilon)$ MMS");
streamplot!(ax, x->Point2(sflow!(du, x, (;pars...,typ=:MMS,order=2), 0.)),
            Alims, Blims)
scatter!(ax, real(Aan), -imag(Aan), marker=:xcross, color=:red, markersize=20);
scatter!(ax, Point2(solve(SteadyStateProblem(sflow!, [0.,0.],
                                             (;pars...,typ=:MMS,order=2)))),
         color=:blue, markersize=16);

ax = Axis(fig[2, 1:2], xlabel=L"State $A_0$", ylabel=L"State $B_0$",
          title=L"$\mathcal{O}(1)$ EMS");
streamplot!(ax, x->Point2(sflow!(du, x, (;pars...,typ=:EMS,order=1), 0.)),
            Alims, Blims)
scatter!(ax, real(Aan), -imag(Aan), marker=:xcross, color=:red, markersize=20);
scatter!(ax, Point2(solve(SteadyStateProblem(sflow!, [0.,0.],
                                             (;pars...,typ=:EMS,order=1)))),
         color=:blue, markersize=16);

ax = Axis(fig[2, 3:4], xlabel=L"State $A_0$", ylabel=L"State $B_0$",
          title=L"$\mathcal{O}(\epsilon)$ EMS");
streamplot!(ax, x->Point2(sflow!(du, x, (;pars...,typ=:EMS,order=2), 0.)),
            Alims, Blims)
scatter!(ax, real(Aan), -imag(Aan), marker=:xcross, color=:red, markersize=20);
scatter!(ax, Point2(solve(SteadyStateProblem(sflow!, [0.,0.],
                                             (;pars...,typ=:EMS,order=2)))),
         color=:blue, markersize=16);

ax = Axis(fig[2, 5:6], xlabel=L"State $A_0$", ylabel=L"State $B_0$",
          title=L"$\mathcal{O}(\epsilon^2)$ EMS");
streamplot!(ax, x->Point2(sflow!(du, x, (;pars...,typ=:EMS,order=3), 0.)),
            Alims, Blims)
scatter!(ax, real(Aan), -imag(Aan), marker=:xcross, color=:red, markersize=20);
scatter!(ax, Point2(solve(SteadyStateProblem(sflow!, [0.,0.],
                                             (;pars...,typ=:EMS,order=3)))),
         color=:blue, markersize=16);

if isdefined(Main, :GLMakie)
   display(scr, fig);
else
    fig   
end
   
# * Assess Stability
spars = (;pars..., typ=:EMS, order=3, Om=2.0);

ss = solve(SteadyStateProblem(sflow!, [0.,0.], spars));
J = ForwardDiff.jacobian(u -> sflow!(similar(u), u, spars, 0.0), ss.u);
display(eigvals(J))
