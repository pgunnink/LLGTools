


struct Vars
    t::Num
    θ::Vector{Num}
    ϕ::Vector{Num}
end

Vars(N) = Vars((@variables t, θ[1:N](t), ϕ[1:N](t))...)

struct RHS
    θ::Vector{Num}
    ϕ::Vector{Num}
end

Base.:+(x::RHS, y::RHS) = RHS(x.θ .+ y.θ, x.ϕ .+ y.ϕ)


function effectiveField(field::Array, vars::Vars, α)
    # field is Nx3 dimensional, given in cartesian coordinates
    N = length(vars.ϕ)
    ϕres = Vector{Num}(undef, N)
    θres = Vector{Num}(undef, N)
    for i in 1:N
        hx, hy, hz = field[i,:]
        θ = vars.θ[i]
        ϕ = vars.ϕ[i]
        ϕres[i] = hz + cos(ϕ) * (-(hx * cot(θ)) + hy * α * csc(θ)) - hy * cot(θ) * sin(ϕ) - hx * α * csc(θ) * sin(ϕ)
        θres[i] = (hy + hx * α * cos(θ)) * cos(ϕ) - hz * α * sin(θ) + (-hx + hy * α * cos(θ)) * sin(ϕ)
    end
    RHS(θres, ϕres)
end


 
function spinTorque(torque::Array, vars::Vars)
    # takes as input an array
    N = length(vars.ϕ)
    ϕres = Vector{Num}(undef, N)
    θres = Vector{Num}(undef, N)
    for i in 1:N
        ϕ = vars.ϕ[i]
        θ = vars.θ[i]
        ϕres[i] = csc(θ) * (torque[i,1] * sin(ϕ) - torque[i,2] * cos(ϕ))
        θres[i] = torque[i,3] * sin(θ) - cos(θ) * (torque[i,1] * cos(ϕ) + torque[i,2] * sin(ϕ))
    end
    return RHS(θres, ϕres)
end

function spinTorque(torquefunc::Function, vars::Vars)
    # takes as input an function with as parameter Particle and should return a 3D array
    mapping = mappingAtoms(unit_cells)
    N = length(mapping)
    res = zeros(N, 3)
    for (i, p) in mapping
        res[i, :] = torquefunc(p)
    end
    spintorque(res, vars)
end


function exchangeCoupling(coupling::Matrix, vars::Vars)
    # coupling is a NxN array which describes which sites couple to which
    # returns two arrays with θ and ϕ right hand sides of the EOM
    ϕres = Vector{Num}(undef, size(coupling)[1])
    θres = Vector{Num}(undef, size(coupling)[1])
    for i in 1:size(coupling)[1]
        tempθ = 0
        tempϕ = 0
        ϕi = vars.ϕ[i]
        θi = vars.θ[i]
        for j in 1:size(coupling)[2]
            ϕj = vars.ϕ[j]
            θj = vars.θ[j]
            tempθ += -coupling[i,j] * sin(θj) * sin(θi - θj)
            tempϕ += coupling[i,j] * (
                cos(θj) - cos(ϕi - ϕj) * cot(θi) * sin(θj) 
            )
        end
        θres[i] = tempθ
        ϕres[i] = tempϕ
    end
    RHS(θres, ϕres)
end

function DMICoupling(coupling::Array, DMI::Array, vars::Vars)
    # coupling is a NxN array which describes which sites couple to which
    # DMI is a 3D vector that describes the DMI field
    # returns two arrays with θ and ϕ right hand sides of the EOM
    ϕres = Vector{Num}(undef, size(coupling)[1])
    θres = Vector{Num}(undef, size(coupling)[1])
    Dx, Dy, Dz = DMI
    for i in 1:size(coupling)[1]
        tempθ = 0
        tempϕ = 0
        ϕi = vars.ϕ[i]
        θi = vars.θ[i]
        for j in 1:size(coupling)[2]
            ϕj = vars.ϕ[j]
            θj = vars.θ[j]

            tempθ += coupling[i,j] * (
                Dz * cos(ϕi - ϕj) * sin(θj) - cos(θj) * (Dx * cos(ϕi) + Dy * sin(ϕi))
            )
            tempϕ += coupling[i,j] * (
                -(cos(θj) * cot(θi) * (Dy * cos(ϕi) - Dx * sin(ϕi))) + 
                -  sin(θj) * (-(cos(ϕj) * (Dy + Dz * cot(θi) * sin(ϕi))) + (Dx + Dz * cos(ϕi) * cot(θi)) * sin(ϕj))
                )
        end
        θres[i] = tempθ
        ϕres[i] = tempϕ
    end
    RHS(θres, ϕres)
end


function setupEOM(vars::Vars, rhs::RHS)
    D = Differential(vars.t)
    [D.(vars.θ) .~ rhs.θ; D.(vars.ϕ) .~ rhs.ϕ]
end




function initBoltzmann(vars::Vars, λ)
    # boltzmann is proportional to θ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    N = length(vars.θ)
    dist = Exponential(λ / 2)
    init_θ = rand(dist, N)
    init_ϕ = rand(Float64, N) * 2π
    return [[vars.θ[i] => init_θ[i] for i in 1:N]; [vars.ϕ[i] => init_ϕ[i] for i in 1:N]]
end