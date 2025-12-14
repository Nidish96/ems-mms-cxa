# * Preamble
using DifferentialEquations
using NonlinearSolve
using LaTeXStrings
using LinearAlgebra
using JLD2

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
analyze = false;
savdats = true;

pars = (c=0.04, w0=2., al=0.1, f=1., Om=0.1);

# * Analyze or Load Data
slowfun(typ,ord) = NonlinearFunction((du,u,p)-> sflow!(du,u,(;pars...,Om=p,
                                                             typ=typ,order=ord)));

typ = :EMS;
ord = 3;

Om0 = 0.1;
Om1 = 4.0;
dOm = 0.08;

if analyze
    u0 = ones(2);
    cpars = (parm=:arclength, nmax=1000, Dsc=:none, save_jacs=true);

    sol, _, _, _, _ = CONTINUATE(u0, slowfun(typ,ord), [Om0, Om1], dOm;
        cpars...);
    stab = [sum(real(eigvals(j)).>0) for j in sol.J];

    # ** Higher Harmonics
    h35s = [[sflow_harms(up[1:2], (;pars...,typ=typ,order=ord,Om=up[3]))...]
           for up in sol.up];
    uh = [[v[1]-im*v[2]; h35...] for (v,h35) in zip(sol.u, h35s)];

    Nt = 128;
    tau = range(0, 2π, Nt+1)[1:Nt];
    et = [exp.(-1im*n*tau) for n in 1:2:5];
    Umaxs = [maximum(abs.(sum(u.*et))) for u in uh];

    mainb = (Oms=sol.p, Umaxs=Umaxs, stab=stab, uh=uh);

    # * Bifurcation Detection
    bifis = findall(stab[1:end-1].!=stab[2:end]);
    dxis = sign.(stab[bifis.+1]-stab[bifis]);  # Direction of bifurcation
    bifis[1] -= 1;
    bifis[2] += 2;

    abhbs  = [];
    wselbs = [];
    Omsb   = [];

    # Parameters for HB on the Slow Flow
    h = 0:3;
    N = 64;

    Nhc = NHC(h);
    _,_, zinds,rinds,iinds = HINDS(2, h);

    cL = I(2Nhc)[:, setdiff(1:2Nhc, iinds[1])];  # Phase constraint
    qamp = 0.1;  # Step size for branch switch

    cpars = (parm=:arclength, nmax=1000, Dsc=:none);
    for (bifi,dxi) in zip(bifis,dxis)
        eVals, eVecs = eigen(sol.J[bifi]);
        ei = findall(imag(eVals).>0);

        eVecs = eVecs.*exp.(-1im*eVecs[1, ei]);  # Phase normalize
        eVec = normalize(eVecs[:, ei]);
        eVal = eVals[ei];

        lu0 = zeros(2Nhc);
        lu0[zinds] = sol.up[bifi][1:2];
        lu0[[rinds[1:2];iinds[1:2]]] = qamp*[real(eVec); imag(eVec)];

        uC = zeros(2Nhc);
        uC[zinds] = sol.up[bifi][1:2];

        w0 = imag(eVal);

        fun = NonlinearFunction((R,uw,p)-> hbslow!(R, [cL*uw[1:end-1];uw[end]],
            (;pars..., Om=p,typ=typ,order=ord), h,N; U0=[uC;w0]));
        prob = NonlinearProblem(fun, [cL'lu0;w0], mainb.Oms[bifi]);
        solb0 = solve(prob, show_trace=Val(false));

        # Continuation
        Om00 = mainb.Oms[bifi];
        Om11 = dxi<0 ? Om0 : Om1;
        dOm1 = dxi<0 ? 0.5 : 0.25;
        solb, _, _, _, _ = CONTINUATE(solb0.u, fun, [Om00, Om11], dOm1;
            cpars...);

        # Also get previous points right "behind" (without deflation)
        fun = NonlinearFunction((R,uw,p)-> hbslow!(R, [cL*uw[1:end-1];uw[end]],
            (;pars..., Om=p,typ=typ,order=ord), h,N));

        solbp = [];
        for i in 1:2
            lprob = NonlinearProblem(fun, solb0.u, mainb.Oms[bifi-i*dxi]);
            lsol = solve(lprob, show_trace=Val(false));
            push!(solbp, myNLSoln([lsol.u;mainb.Oms[bifi-i*dxi]]));
        end

        solb = [solbp[end:-1:1]..., solb...];
        
        push!(abhbs, [cL*up[1:end-2] for up in solb.up]);
        push!(wselbs, [up[end-1] for up in solb.up]);
        push!(Omsb, [up[end] for up in solb.up]);

        display("Done $(dxi)")
    end

    # * Get Peak Responses in Time Domain
    N = 64;

    Umaxs = [];
    for bi = 1:length(abhbs)
        Umax = zeros(length(abhbs[bi]));
        for iw = 1:length(abhbs[bi])
            uh1t = eachrow(AFT(reshape(abhbs[bi][iw], 2,Nhc)', h,N, :f2t));
            hfun = (u)-> sflow_harms(u, (;pars...,typ=typ,order=ord,Om=Omsb[bi][iw]));
            h35t = hfun.(uh1t);
            uht = [[v[1]-im*v[2]; h35...] for (v,h35) in zip(uh1t, h35t)];
            Umaxt = [maximum(abs.(sum(u.*et))) for u in uht];
            Umax[iw] = maximum(Umaxt);
        end
        push!(Umaxs, Umax);
    end

    bifbs = [(Oms=om, Umaxs=umax) for (om,umax) in zip(Omsb,Umaxs)];

    if savdats
        @save "./DATS/D1_$(typ)_$(ord).jld2" mainb bifbs
    end
else
    @load "./DATS/D1_$(typ)_$(ord).jld2" mainb bifbs
end

# * Load Transient Data
trdat = load("./DATS/D0_trdat.jld2");
trdat = (; [Symbol(k)=>v for (k,v) in trdat]...);

# * Plot
set_theme!(theme_latexfonts())
fsz = 22;
fig = Figure(fontsize=fsz);
if !isdefined(Main, :scr) && Makie.current_backend()==GLMakie
   scr = GLMakie.Screen();
end

ax = Axis(fig[1, 1], xlabel=L"Excitation Frequency (rad/s)",
    ylabel="Peak Absolute Response (m)");
scatterlines!(ax, mainb.Oms./(mainb.stab.==0), mainb.Umaxs)
scatterlines!(ax, mainb.Oms./(mainb.stab.==2), mainb.Umaxs)

for bifb in bifbs
    scatterlines!(ax, bifb.Oms, bifb.Umaxs)
end

lines!(ax, trdat.Oms, trdat.Umaxs, color=:black);

xlims!(Om0, Om1);

if Makie.current_backend()==GLMakie
   display(scr, fig);
else
   display(fig)
end
   
