function exchangeEffectiveField(coupling::Matrix, u)
    # coupling is a NxN array which describes which sites couple to which
    # returns an array with effective field
    Heff = zeros(Float64, size(coupling)[1], 3)
    
    for i in 1:size(coupling)[1]
        for j in 1:size(coupling)[1]
            @. Heff[i,:] += coupling[i,j] .* u[j,:]
        end
    end
    Heff
end


function exchangeEffectiveField!(Heff, coupling::Matrix, u)
    # coupling is a NxN array which describes which sites couple to which
    # returns an array with effective field    
    for i in 1:size(coupling)[1]
        for j in 1:size(coupling)[1]
            uj = @view u[j,:]
            cij = @view coupling[i,j]
            @. Heff[i,:] += cij .* uj
        end
    end
end

function DMIEffectiveField(coupling::Matrix, Dfield::Vector, u)
    # coupling is a NxN array which describes which sites couple to which
    # returns an array with effective field
    Heff = zeros(Float64, size(coupling)[1], 3)
    
    for i in 1:size(coupling)[1]
        for j in 1:size(coupling)[1]
            Heff[i,:] .-= coupling[i,j] * cross(Dfield, u[j,:])
        end
    end
    Heff
end


function magneticEffectiveField(Hfield::Vector, N)
    Heff = zeros(Float64, N, 3)
    for i in 1:N
        @. Heff[i,:] += Hfield
    end
    return Heff
end

function magneticEffectiveField!(Heff, Hfield::Vector, N)
    for i in 1:N
        @. Heff[i,:] += Hfield
    end
end

function spinTorqueEffectiveField(Js, u)
    Heff = zeros(Float64, size(Js)[1], 3)
    for i in 1:size(Js)[1]
        Heff[i,:] = cross(Js[i,:], u[i,:])
    end
    return Heff
end

# function generateEOM(Heff, α, vars::VarsCart)
#     res = zeros(Num, length(vars.x), 3)
#     for i in 1:length(vars.x)
#         S = [vars.x[i], vars.y[i], vars.z[i]]
#         res[i, :] = - cross(S, Heff[i,:]) - α * cross(S, cross(S, Heff[i,:]))
#     end
#     D = Differential(vars.t)
#     return [D.(vars.x) .~ res[:, 1]; D.(vars.y) .~ res[:, 2]; D.(vars.z) .~ res[:, 3]]
# end

function initRandomSpin(N)
    # boltzmann is proportional to θ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    ϕ = rand(N) * 2π
    θ = rand(N) * π
    mx = sin.(θ) .* cos.(ϕ)
    my = sin.(θ) .* sin.(ϕ)
    mz = cos.(θ)
    m0 = zeros(N, 3)
    m0[:,1] = mx
    m0[:,2] = my
    m0[:,3] = mz
    return m0
end

function initBoltzmann(N, λ)
    # boltzmann is proportional to θ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    dist = Exponential(λ / 2)
    θ = rand(dist, N)
    ϕ = rand(Float64, N) * 2π
    mx = sin.(θ) .* cos.(ϕ)
    my = sin.(θ) .* sin.(ϕ)
    mz = cos.(θ)
    m0 = zeros(N, 3)
    m0[:,1] = mx
    m0[:,2] = my
    m0[:,3] = mz
    return m0

end

