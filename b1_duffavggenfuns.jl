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

analyze = true;
savdats = true;

savfigs = false;
if savfigs
    CairoMakie.activate!();
else
    GLMakie.activate!();
end

cfgs = [(zt = 0.25e-2, w0 = 2.0, al = 0.01, f = [0.25, 0.5, 1.0], Om=1.9),
        (zt = 0.25e-1, w0 = 2.0, al = 0.2, f = [0.4, 0.8, 1.6], Om=1.9)];
ci = 2;
pars = cfgs[ci];

# * Load Data
typs = [:MMS, :MMS, :EMS, :EMS];
ords = [1, 2, 1, 2];

@load "./DATS/B0_duffss.jld2" slowsols slowh35s

