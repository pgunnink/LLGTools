struct VarsCart
    t::Num
    x::Vector{Num}
    y::Vector{Num}
    z::Vector{Num}
end

function VarsCart(N)
    t = (@variables t)[1]
    x = (@variables x[1:N])[1] |> Symbolics.scalarize
    y = (@variables y[1:N])[1] |> Symbolics.scalarize
    z = (@variables z[1:N])[1] |> Symbolics.scalarize
    VarsCart(t, x, y, z)
end
mi(i, vars::VarsCart) = [vars.x[i]; vars.y[i]; vars.z[i]]

# struct RHS
#     x::Vector{Num}
#     y::Vector{Num}
#     z::Vector{Num}
# end

# Base.:+(x::RHS, y::RHS) = RHS(x.x .+ y.x, x.y .+ y.y, x.z .+ y.z)

# we use the conventions of SPIRIT: dtm = - m × Beff and A = -Beff

function exchangeEffectiveField(coupling::Matrix, vars::VarsCart)
    # coupling is a NxN array which describes which sites couple to which
    # returns an array with effective field
    Heff = zeros(Num, size(coupling)[1], 3)

    for i in 1:size(coupling)[1]
        for j in 1:size(coupling)[1]
            Heff[i, :] .+= coupling[i, j] .* [vars.x[j], vars.y[j], vars.z[j]]
        end
    end
    Heff
end

function exchangeXYEffectiveField(coupling::Matrix, vars::VarsCart)
    # coupling is a NxN array which describes which sites couple to which
    # returns an array with effective field
    Heff = zeros(Num, size(coupling)[1], 3)

    for i in 1:size(coupling)[1]
        for j in 1:size(coupling)[1]
            Heff[i, :] .+= coupling[i, j] .* [vars.x[j], vars.y[j], 0]
        end
    end
    Heff
end

function spinTorqueField(field, vars::VarsCart)
    # spin torque is defined as m × (m × J)
    Js = zeros(Num, length(vars.x)[1], 3)
    for i in 1:size(Js)[1]
        Js[i, :] .+= cross([vars.x[i], vars.y[i], vars.z[i]], field[i, :])
    end
    Js
end

function DMIEffectiveField(coupling::Matrix, Dfield::Vector, vars::VarsCart)
    # coupling is a NxN array which describes which sites couple to which
    # returns an array with effective field
    Heff = zeros(Num, size(coupling)[1], 3)
    for i in 1:size(coupling)[1]
        for j in 1:size(coupling)[1]
            Heff[i, :] .-= coupling[i, j] * cross(Dfield, mi(j, vars))
        end
    end
    Heff
end

function anistropyEnergy(anistropy::Vector, vars::VarsCart)
    E = Num(0)
    for i in 1:length(vars.x)
        E -= sum(anistropy .* mi(i, vars) .^ 2)
    end
    return E
end

function anistropyFieldZ(K, vars::VarsCart)
    Heff = zeros(Num, length(vars.x), 3)
    for i in 1:length(vars.x)
        Heff[i, 3] = -2K * vars.z[i]
    end
    Heff
end

function exchangeEnergy(coupling::Matrix, vars::VarsCart)
    # we count every pair twice, so multiply by 0.5
    E = Num(0)
    for i in 1:length(vars.x)
        for j in 1:length(vars.x)
            E += -0.5coupling[i, j] * mi(i, vars) ⋅ mi(j, vars)
        end
    end
    return E
end

function magneticEnergy(B::Vector, vars::VarsCart)
    E = Num(0)
    for i in 1:length(vars.x)
        E += -B ⋅ mi(i, vars)
    end
    return E
end

function DMIEnergy(coupling::Matrix, Dfield::Vector, vars::VarsCart)
    # we count every pair twice, so multiply by 0.5
    E = Num(0)
    for i in 1:length(vars.x)
        for j in 1:length(vars.x)
            E += -0.5coupling[i, j] * Dfield ⋅ cross(mi(i, vars), mi(j, vars))
        end
    end
    return E
end

function DipolarEnergy(R, vars::VarsCart)
    E = Num(0)
    for i in 1:length(vars.x)
        Si = mi(i, vars)
        for j in 1:length(vars.x)
            if i == j
                continue
            end
            Sj = mi(j, vars)
            Rij = R[i] .- R[j]
            Rijnorm = Rij .^ 2 |> sum |> sqrt
            Rijhat = Rij ./ Rijnorm
            E -= 0.5 * Rijnorm^(-3) * (3 * (Si ⋅ Rijhat) * (Sj ⋅ Rijhat) - Si ⋅ Sj)
        end
    end
    return E
end

function effectiveFieldFromEnergy(E, vars::VarsCart)
    Heff = zeros(Num, length(vars.x), 3)
    Heff[:, 1] = -Symbolics.gradient(E, vars.x)
    Heff[:, 2] = -Symbolics.gradient(E, vars.y)
    Heff[:, 3] = -Symbolics.gradient(E, vars.z)

    return Heff

end

# custom functions to make sure the substitute works:

istree(c::Complex{Num}) = true
operation(::Complex{Num}) = Complex{Num}
arguments(c::Complex{Num}) = [real(c), imag(c)]


function HamiltonianFromHeff(Heff, vars)
    N = length(vars.x)
    δm = (@variables δm[1:N])[1] |> Symbolics.scalarize
    δmconj = (@variables δmconj[1:N])[1] |> Symbolics.scalarize
    repl_rules = Dict([])
    for i in 1:N
        repl_rules[vars.x[i]] = (δm[i] + δmconj[i]) / 2
        repl_rules[vars.y[i]] = (δm[i] - δmconj[i]) * -1im / 2
        repl_rules[vars.z[i]] = 1
    end
    # Heffδ = similar(Heff)
    # for i in 1:length(Heff)
    #     Heffδ[i] = substitute(Heff[i], repl_rules)
    # end
    eom_δm = zeros(Complex{Num}, N)
    eom_δm_imag = zeros(Num, N)

    eom_δmconj = zeros(Complex{Num}, N)
    eom_δmconj_imag = zeros(Num, N)

    for i in 1:N
        eom_m = cross(LLG.mi(i, vars), Heff[i, :])
        eom_δm[i] = substitute((eom_m[1] + 1im * eom_m[2]), repl_rules)
        # eom_δm_imag[i] = substitute((eom_m[1] + 1im * eom_m[2]) |> imag, repl_rules) 
        eom_δmconj[i] = substitute((eom_m[1] - 1im * eom_m[2]), repl_rules)
        # eom_δmconj_imag[i] = substitute((eom_m[1] - 1im * eom_m[2]) |> imag, repl_rules) 
    end
    # A = zeros(Complex{Num}, N, N)
    # B = zeros(Complex{Num}, N, N)
    # for i in 1:N
    #     for j in 1:N
    #         A[i,j] = Symbolics.linear_expansion(real(eom[i,1]), δm[j])[1] + im * Symbolics.linear_expansion(imag(eom[i,1]), δm[j])[1]
    #         B[i,j] = Symbolics.linear_expansion(real(eom[i, 2]), δmconj[j])[1] + im * Symbolics.linear_expansion(imag(eom[i, 2]), δmconj[j])[1]
    #     end
    # end
    # A, B
    i, j = (2, 3)
    eom_m = cross(LLG.mi(i, vars), Heff[i, :])[1]
    # @info eom_m
    # @time Symbolics.linear_expansion(eom_m, vars.z[j])[1]
    @info Symbolics.value(eom_δm[i])
    # @time Symbolics.linear_expansion(eom_δm[i], δm[j])[1]
end

function magneticEffectiveField(Hfield::Vector, vars::VarsCart)
    Heff = zeros(Num, length(vars.x), 3)
    for i in 1:length(vars.x)
        Heff[i, :] = Hfield
    end
    return Heff
end

function generateEnergy(E, vars::VarsCart)
    eval(build_function(E, [vars.x; vars.y; vars.z]))
end

function generateEnergyFromHeff(Heff, vars::VarsCart)
    E = Num(0)
    for i in 1:length(vars.x)
        E -= Heff[:, i] ⋅ mi(i, vars)
    end
    return eval(build_function(E, [vars.x; vars.y; vars.z]))
end

function generateEOM(Heff, α, vars::VarsCart)
    res = zeros(Num, length(vars.x), 3)
    for i in 1:length(vars.x)
        S = [vars.x[i], vars.y[i], vars.z[i]]
        res[i, :] = -cross(S, Heff[:, i]) - α * cross(S, cross(S, Heff[:, i]))
    end
    D = Differential(vars.t)
    return [D.(vars.x) .~ res[:, 1]; D.(vars.y) .~ res[:, 2]; D.(vars.z) .~ res[:, 3]]
end


function generateA(Heff, vars::VarsCart, p=[])
    eval(build_function(-Heff, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true})[1])
end
function generateA!(Heff, vars::VarsCart, p=[]; parallel=false)
    if parallel
        eval(build_function(-Heff, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true}, parallel=Symbolics.MultithreadedForm())[2])
    else
        eval(build_function(-Heff, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true})[2])
    end
end

function generateA!(Heff, vars::VarsCart, α::Number, p=[]; parallel=false)
    totalH = -Heff
    for i in 1:length(vars.x)
        totalH[i, :] .-= α * cross([vars.x[i], vars.y[i], vars.z[i]], Heff[i, :])
    end
    if parallel
        eval(build_function(totalH, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true}, parallel=Symbolics.MultithreadedForm())[2])
    else
        eval(build_function(totalH, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true})[2])
    end
end

function generateA!(Heff, auxfield, vars::VarsCart, α::Number, p=[]; parallel=false)
    totalH = -Heff
    for i in 1:length(vars.x)
        totalH[i, :] .-= α * cross([vars.x[i], vars.y[i], vars.z[i]], Heff[i, :])
    end
    if parallel
        build_function(totalH + auxfield, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{false}, parallel=Symbolics.MultithreadedForm())[2]
    else
        build_function(totalH + auxfield, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{false})[2]
    end
end




function generateA!(Heff, auxfield, vars::VarsCart, α::Vector{Number}, p=[]; parallel=false)
    totalH = -Heff
    for i in 1:length(vars.x)
        totalH[i, :] .-= α[i] * cross([vars.x[i], vars.y[i], vars.z[i]], Heff[i, :])
    end
    if parallel
        eval(build_function(totalH + auxfield, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true}, parallel=Symbolics.MultithreadedForm())[2])
    else
        eval(build_function(totalH + auxfield, [vars.x; vars.y; vars.z], p, vars.t, expression=Val{true})[2])
    end
end

function noiseLLG(vars::VarsCart)
    # eval(build_function(0, [vars.x; vars.y; vars.z], vars.t)[1])
    (u, t) -> 0
end

function generateEOM(Heff, Js::Array, α, vars::VarsCart)
    res = zeros(Num, length(vars.x), 3)
    for i in 1:length(vars.x)
        S = [vars.x[i], vars.y[i], vars.z[i]]
        res[i, :] = -cross(S, Heff[i, :]) - α * cross(S, cross(S, Heff[i, :])) + cross(S, cross(S, Js[i, :]))
    end
    D = Differential(vars.t)
    return [D.(vars.x) .~ res[:, 1]; D.(vars.y) .~ res[:, 2]; D.(vars.z) .~ res[:, 3]]
end


function initRandomSpinToolkit(vars::VarsCart, λ)
    # boltzmann is proportional to θ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    N = length(vars.x)
    ϕ = rand(N) * 2π
    θ = rand(N) * π
    mx = sin.(θ) .* cos.(ϕ)
    my = sin.(θ) .* sin.(ϕ)
    mz = cos.(θ)
    return [
        [vars.x[i] => mx[i] for i in 1:N]
        [vars.y[i] => my[i] for i in 1:N]
        [vars.z[i] => mz[i] for i in 1:N]
    ]
end

function initRandomSpin(N::Int, λ)
    ϕ = rand(N) * 2π
    θ = rand(N) * π
    mx = sin.(θ) .* cos.(ϕ)
    my = sin.(θ) .* sin.(ϕ)
    mz = cos.(θ)
    return [
    mx my mz
]
end

function initRandomSpin(vars::VarsCart, λ)
    # boltzmann is proportional to θ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    N = length(vars.x)
    initRandomSpin(N, λ)
end

function initBoltzmannSpin(vars::VarsCart, λ)
    # permutation allows for easy switching the axis around
    # boltzmann is proportional to λ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    N = length(vars.x)
    ϕ = rand(N) * 2π
    dist = Exponential(λ / 2)
    p = rand(dist, N)
    θ = acos.(1 .- 2p)
    mx = sin.(θ) .* cos.(ϕ)
    my = sin.(θ) .* sin.(ϕ)
    mz = cos.(θ)
    m = zeros(Float64, N, 3)
    m[:, 1] .= mx
    m[:, 2] .= my
    m[:, 3] .= mz
    return m
end

function initBoltzmannSpinToolkit(vars::VarsCart, λ)
    # boltzmann is proportional to θ (for now)
    # the angle ϕ is chosen randomly between 0 and 2pi
    N = length(vars.x)
    ϕ = rand(N) * 2π
    dist = Exponential(λ / 2)
    p = rand(dist, N)
    θ = acos.(1 .- 2p)
    mx = sin.(θ) .* cos.(ϕ)
    my = sin.(θ) .* sin.(ϕ)
    mz = cos.(θ)
    return [
        [vars.x[i] => mx[i] for i in 1:N]
        [vars.y[i] => my[i] for i in 1:N]
        [vars.z[i] => mz[i] for i in 1:N]
    ]
end


