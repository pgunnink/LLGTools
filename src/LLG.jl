module LLG
# using DifferentialEquations
using StaticArrays: reshape
using Base: Float64
using ModelingToolkit
using Random, Distributions
import LinearAlgebra: cross, similar, det, ⋅
import Symbolics: build_function
using GLMakie


## diffeq stuff:


import Symbolics: istree, operation, arguments


export LatticeDescription, generateLatticeShape, Atom, positionsLattice, mappingAtoms, coupleNN, coupleInternal, isless

export plot_atoms, plot_atoms!, plot_connections, plot_connections!, absPositionsParticle
export dmi_interaction_honeycomb, nn_coupling_honeycomb, nnn_coupling_honeycomb
export initBoltzmann


# export custom alg
# export LLGProblem, SIB

# export exchangeCoupling, DMICoupling, setupEOM

include("lattice.jl")
# include("diffeq.jl")
include("plotting.jl")
include("graphene.jl")
include("diffeq_toolkit.jl")
# include("custom_cache.jl")
# include("custom_alg.jl")
# Write your package code here.
end
