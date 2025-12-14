# * 
using DifferentialEquations
using NonlinearSolve
using GLMakie
using LaTeXStrings

using Revise
includet("/home/nidish/Documents/Research/b_Programming/juliajim/src/HARMONIC.jl")
includet("/home/nidish/Documents/Research/b_Programming/juliajim/src/CONTINUATION.jl")
includet("./linearsys.jl");

# * Simulations
pars = (zt = 0.1e-2, w0 = 2.0, f = 1.0, Om=1.9);

ncyc = 100;
Nom = 200;
thv = tanh.(range(-0.5, 0.5, Nom));
Oms = (thv.-thv[1])/(thv[end]-thv[1])*(2pars.w0-0.05pars.w0) .+ 0.05pars.w0;
taumax = 2π*ncyc;

Aan = [ansol((;pars..., Om=om)) for om in Oms];

trerrs = zeros(Nom, 5);
sserrs = zeros(Nom, 5);
sssols = zeros(Complex, Nom, 5);

for (iw, Om) in enumerate(Oms)
    # Original ODE Solution
    prob = ODEProblem(rocfun!, [0., 0.], [0., taumax/Om], (;pars...,Om=Om));
    local sol = solve(prob, RK4());

    local uan = abs(Aan[iw])cos.(Om*sol.t .+ angle(Aan[iw]));
    local utr = [u[1] for u in sol.u]-uan;

    # Slow Flow Solutions

    # MMS
    for oi in 1:2
    	local sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,taumax/Om], (;pars...,Om=Om,typ=:MMS,order=oi));
        sfsol = solve(sprob, RK4());

        local ssprob = SteadyStateProblem(sflow!, [0.0,0.0], (;pars...,Om=Om,typ=:MMS,order=oi));
        local ss = solve(ssprob);

        ut = [u[1]cos(Om*sol.t[i])+u[2]sin(Om*sol.t[i]) for (i, u) in enumerate(sfsol(sol.t))];
        utr_sf = ut - (ss[1]cos.(Om*sol.t)+ss[2]sin.(Om*sol.t));

        trerrs[iw, oi] = utr'utr_sf/norm(utr)/norm(utr_sf);
        sserrs[iw, oi] = norm(ss.u-[real(Aan[iw]), imag(Aan[iw])]);
        sssols[iw, oi] = ss[1]-im*ss[2];
    end

    # EMS
    for oi in 1:3
    	local sprob = ODEProblem(sflow!, [0.0,0.0], [0.0,taumax/Om], (;pars...,Om=Om,typ=:EMS,order=oi));
        sfsol = solve(sprob, RK4());

        ssprob = SteadyStateProblem(sflow!, [0.0,0.0], (;pars...,Om=Om,typ=:EMS,order=oi));
        local ss = solve(ssprob);

        ut = [u[1]cos(Om*sol.t[i])+u[2]sin(Om*sol.t[i]) for (i, u) in enumerate(sfsol(sol.t))];
        utr_sf = ut - (ss[1]cos.(Om*sol.t)+ss[2]sin.(Om*sol.t));

        trerrs[iw, 2+oi] = utr'utr_sf/norm(utr)/norm(utr_sf);
        sserrs[iw, 2+oi] = norm(ss.u-[real(Aan[iw]), -imag(Aan[iw])]);
        sssols[iw, 2+oi] = ss[1]-im*ss[2];
    end
end

# * Plotting

meths = [L"$\mathcal{O}(1)$ MMS", L"$\mathcal{O}(\epsilon)$ MMS",
         L"$\mathcal{O}(1)$ EMS", L"$\mathcal{O}(\epsilon)$ EMS",
         L"$\mathcal{O}(\epsilon^2)$ EMS"];

stls = [(color=:blue, linestyle=:solid),
         (color=:blue, linestyle=:dash),
         (color=:red, linestyle=:solid),
         (color=:red, linestyle=:dash),
         (color=:red, linestyle=:dot)];

set_theme!(theme_latexfonts())
fsz = 28;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && isdefined(Main, :GLMakie)
   scr = GLMakie.Screen();
end

ax = Axis(fig[1:2, 1], ylabel="Steady State Response",
          yscale=log10);
scatterlines!(ax, Oms, abs.(Aan), color=:black, label="Analytical")
for i in eachindex(meths)
    lines!(ax, Oms, abs.(sssols[:,i]), label=meths[i]; stls[i]...)
end
axislegend()

ax = Axis(fig[3, 1],
          ylabel="Relative SS Error", yscale=log10);
for i in eachindex(meths)
    lines!(ax, Oms, abs.(sserrs[:,i]./sssols[:,i]), label=meths[i]; stls[i]...)
end

ax = Axis(fig[4, 1], xlabel="Excitation Frequency (rad/s)",
          ylabel="Relative Transient Error", yscale=log10);
for i in eachindex(meths)
    lines!(ax, Oms, abs.(1 .-trerrs[:,i]), label=meths[i]; stls[i]...)
end

if isdefined(Main, :GLMakie)
   display(scr, fig);
else
    fig   
end
