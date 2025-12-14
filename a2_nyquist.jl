# * 
using DifferentialEquations
using NonlinearSolve
using GLMakie
using CairoMakie
using LaTeXStrings
using ForwardDiff

using Revise
includet("/home/nidish/Documents/Research/b_Programming/juliajim/src/HARMONIC.jl")
includet("/home/nidish/Documents/Research/b_Programming/juliajim/src/CONTINUATION.jl")
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
