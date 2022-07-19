using LLG
using ModelingToolkit
using BenchmarkTools


l = LatticeDescription([[1, 0], [0.5, 0.5 * sqrt(3)]], Dict(Atom(1) => [0, 0], Atom(2) => [0, 1 / sqrt(3)]))
a = 20
squ = pos -> (-a / 2 < pos[1] < a / 2) && (-a / 2 < pos[2] < a / 2)
unit_cells = generateLatticeShape(l, squ)
pos, mapping = LLG.positionsLattice(unit_cells, l)

couplingNNN_matr = coupleNN(unit_cells, dmi_interaction_honeycomb)
couplingNN_matr = coupleNN(unit_cells, LLG.nn_coupling_honeycomb)
N = size(couplingNN_matr)[1]
vars = LLG.VarsCart(N)
J = 1
Dfield = [0, 0, 0.05 ]
H = 1.
# EJ = LLG.exchangeEnergy(J * couplingNN_matr, vars)
# EDMI = LLG.DMIEnergy(couplingNNN_matr, Dfield, vars)
# EB = LLG.magneticEnergy([0,0,H], vars)
# E  = EJ + EDMI + EB
# Ef = LLG.generateEnergy(E, vars)

HJ = LLG.exchangeEffectiveField(J * couplingNN_matr, vars)
HDMI = LLG.DMIEffectiveField(couplingNNN_matr, Dfield, vars)
HH = LLG.magneticEffectiveField([0,0,H], vars)

Heff = HJ + HDMI + HH
torque_field = zeros(Num, N, 3)
for i in 1:N
    torque_field[i,3] = .2
end
torque = LLG.spinTorqueField(torque_field, vars)
# Heff = LLG.effectiveFieldFromEnergy(E, vars)
α = 1e-2

##
f = LLG.generateA!(Heff, torque, vars, α);
fpar = LLG.generateA!(Heff, torque, vars, α; parallel = true);


m = zeros(Float64, N, 3)
u0 = similar(m)

f(u0, m, (), 0.)
fpar(u0, m, (), 0.)

##

@benchmark f($u0, $m, (), 0.)


@benchmark fpar($u0, $m, (), 0.)

