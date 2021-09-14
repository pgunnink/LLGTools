struct LLGProblem{uType,tType,P,NP,K,ND} <: AbstractSDEProblem{uType,tType,true,ND}
    f::Function
    g::Function
    α::Number
    vars::VarsCart
    u0::uType
    tspan::tType
    kbT::Number
    p::P
    noise_rate_prototype::ND
    noise::NP
    kwargs::K
    seed::UInt64
    
    @add_kwonly function LLGProblem(A, α, vars, u0, tspan, p = NullParameters() ; 
    kbT = 0., noise_rate_prototype = nothing,
        noise = nothing, seed = UInt64(0),
        kwargs...)
        g = LLG.noiseLLG(vars)
        f = convert(SciMLBase.SDEFunction{true}, A, g)
        _tspan = promote_tspan(tspan)
        new{typeof(u0),typeof(_tspan),typeof(p),typeof(noise),typeof(kwargs),typeof(nothing)}(f, g, α, vars, u0, tspan, kbT, p, noise_rate_prototype, noise, kwargs, seed)
    end
end


abstract type LLGAlgorithm <: StochasticDiffEqAlgorithm end
struct SIB <: LLGAlgorithm end

alg_order(alg::SIB) = 1
alg_compatible(prob::LLGProblem,alg::LLGAlgorithm) = true

@cache struct SIBCache{uType} <: StochasticDiffEqMutableCache
    up::uType
    Acache::uType
    N::Int64
    kbT::Float64
    α::Float64
    detM::Vector{Float64}
    detMx::Vector{Float64}
    detMy::Vector{Float64}
    detMz::Vector{Float64}
end
    
function alg_cache(alg::SIB, prob::LLGProblem, u, dW, dZ, p, rate_prototype, noise_rate_prototype, jump_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, f, t, dt, iip)
    Acache = similar(u)
    up = similar(u)
    N = size(u)[1]
    detM = zeros(Float64, N)
    detMx = zeros(Float64, N)
    detMy = zeros(Float64, N)
    detMz = zeros(Float64, N)
    SIBCache(up, Acache, N, prob.kbT, prob.α, detM, detMx, detMy, detMz)
end


function DiffEqBase.__solve(prob::LLGProblem,
    alg::LLGAlgorithm,
    timeseries = [],ts = [],ks = nothing, # needed for variable rate
recompile::Type{Val{recompile_flag}} = Val{true};
    kwargs...) where recompile_flag
    integrator = DiffEqBase.__init(prob, alg, timeseries, ts, recompile;kwargs...)
    solve!(integrator)
    integrator.sol
end

@fastmath @muladd function perform_step!(integrator, cache::SIBCache)
    @unpack t, dt, u, f, W = integrator
    # @unpack M, Mx, My, Mz, a, Acache, N, up  = cache
    @unpack Acache, N, up, α, kbT, detM, detMx, detMy, detMz = cache

    # do the predictor step:
    f(Acache, u, t)
    @. Acache = .5dt * (Acache + W.dW * √(2α * kbT))
    
    fdet(x11, x21, x31, x12, x22, x32, x13, x23, x33) = x11 * (x22 * x33 - x23 * x32) - x12 * (x21 * x33 - x23 * x31) + x13 * (x21 * x32 - x22 * x31)
    
    combined_cross(uloc, Ax, Ay, Az) = (uloc[1] + uloc[2] * Az - Ay * uloc[3], uloc[2] + uloc[3] * Ax - uloc[1] * Az, uloc[3] + uloc[1] * Ay - uloc[2] * Ax)

    @inbounds for i in 1:N
        uloc = @view u[i,:]
        Ax, Ay, Az = @view Acache[i,:]
        a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
        detMx[i] =  fdet(a1, a2, a3, -Az, 1, Ax, Ay, -Ax, 1)
        detMy[i] = fdet(1, Az, -Ay, a1, a2, a3, Ay, -Ax, 1)
        detMz[i] = fdet(1, Az, -Ay, -Az, 1, Ax, a1, a2, a3)
        detM[i] = 1 / fdet(1, Az, -Ay, -Az, 1, Ax, Ay, -Ax, 1)
    end
    @. up[:,1] = detMx * detM
    @. up[:,2] = detMy * detM
    @. up[:,3] = detMz * detM

    # # and the final step:
    @. up = .5 * ( up + u)
    f(Acache,  up, t + .5dt) 
    @. Acache = .5dt * (Acache + W.dW * √(2α * kbT))

    @inbounds for i in 1:N
        uloc = @view u[i,:]
        Ax, Ay, Az = @view Acache[i,:]
        a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
        detMx[i] =  fdet(a1, a2, a3, -Az, 1, Ax, Ay, -Ax, 1)
        detMy[i] = fdet(1, Az, -Ay, a1, a2, a3, Ay, -Ax, 1)
        detMz[i] = fdet(1, Az, -Ay, -Az, 1, Ax, a1, a2, a3)
        detM[i] = 1 / fdet(1, Az, -Ay, -Az, 1, Ax, Ay, -Ax, 1)
    end
    @. u[:,1] = detMx * detM
    @. u[:,2] = detMy * detM
    @. u[:,3] = detMz * detM
    # @inbounds for i in 1:N
    #     # Ax, Ay, Az = @view Acache[:,i]
    #     # uloc = @view u[:,i]
    #     # # a = SVector{3}(uloc + cross(uloc, @view Acache[:,i]))
    #     # a .= uloc .+ cross(uloc, @view Acache[:,i])
    #     # # up[:,i] .= calc_up(Ax, Ay, Az, a)
    #     calc_up!(@view(up[:,i]), @view( Acache[:,i]), @view(u[:,i]))
# end
end
# @muladd function perform_step!(integrator, cache::SIBCache)
#     @unpack t, dt, u, f, W = integrator
#     # @unpack M, Mx, My, Mz, a, Acache, N, up  = cache
#     @unpack Acache, N, up, α, kbT = cache

#     # do the predictor step:
#     f(Acache, u, t)
#     Acache .+= W.dW * √(2α * kbT)
#     Acache .*= .5dt
    
#     fdet(x11, x21, x31, x12, x22, x32, x13, x23, x33) = x11 * (x22 * x33 - x23 * x32) - x12 * (x21 * x33 - x23 * x31) + x13 * (x21 * x32 - x22 * x31)
    
#     combined_cross(uloc, Ax, Ay, Az) = (uloc[1] + uloc[2] * Az - Az * uloc[2], uloc[2] + uloc[3] * Ay - uloc[1] * Az, uloc[3] + uloc[1] * Ay - uloc[2] * Ax)
#     function calc_up!(res, Acachei, uloc)
#         Ax, Ay, Az = Acachei
#         # a = uloc .+ cross(uloc, Acachei)
#         a1, a2, a3 = combined_cross(uloc, Ax, Ay, Az)
#         detMx =  fdet(a1, -Az, Ay, a2, 1, -Ax, a3, Ax, 1)
        
#         detMy = fdet(1, a1, Ay, Az,  a2, -Ax, -Ay, a3, 1)
#         detMz = fdet(1, -Az, a1, Az, 1, a2, -Ay, Ax, a3)
#         detM = fdet(1, -Az, Ay, Az, 1, -Ax, -Ay, Ax, 1)

#         # detMx = a[1] * (1 + Ax^2) + Az * (a[2] * Ax - a[3]) + Ay * (a[2] * Ax - a[3])
#         # detMy = det(@SMatrix  [  1  a[1] Ay 
#         #         Az  a[2] -Ax 
#         #         -Ay a[3] 1
#         #         ] )
#         # detMz = det(@SMatrix  [  1 -Az a[1] 
#         #     Az 1 a[2] 
#         #     -Ay Ax a[3]
#         #     ] )
#         # detM = det(@SMatrix [1 -Az Ay 
#         # Az 1 -Ax 
#         # -Ay Ax 1])

#         res .= [detMx / detM, detMy / detM, detMz / detM]
#     end

#     @inbounds for i in 1:N
#         calc_up!(@view(up[:,i]), @view( Acache[:,i]), @view(u[:,i]))
#     end
        
#     # and the final step:
#     f(Acache, .5 * (up + u), t + dt / 2) 
#     Acache += W.dW * √(2α * kbT)
#     Acache *= dt / 2 
#     @inbounds for i in 1:N
#         # Ax, Ay, Az = @view Acache[:,i]
#         # uloc = @view u[:,i]
#         # # a = SVector{3}(uloc + cross(uloc, @view Acache[:,i]))
#         # a .= uloc .+ cross(uloc, @view Acache[:,i])
#         # # up[:,i] .= calc_up(Ax, Ay, Az, a)
#         calc_up!(@view(u[:,i]), @view( Acache[:,i]), @view(u[:,i]))
#     end
# end
# @muladd function perform_step!(integrator, cache::SIBCache)
#     @unpack t, dt, u, f, W = integrator
#     # @unpack M, Mx, My, Mz, a, Acache, N, up  = cache
#     @unpack Acache, N, up, α, kbT  = cache

#     # do the predictor step:
#     f(Acache, u, t)
#     Acache += W.dW * √(2α * kbT)
#     Acache *= dt / 2 
    
#     @inbounds for i in 1:N
#         Ax, Ay, Az = @view Acache[i,:]
#         uloc = @view u[i,:]
#         a = SVector{3}(uloc + cross(uloc, Acache[i,:]))
#         M = @SMatrix [1 -Az Ay 
#             Az 1 -Ax 
#         -Ay Ax 1]
        
#         Mx = @SMatrix [  a[1] -Az Ay 
#         a[2] 1 -Ax 
#         a[3] Ax 1
#         ]
#         My = @SMatrix  [  1  a[1] Ay 
#                 Az  a[2] -Ax 
#                 -Ay a[3] 1
#                 ]
#         Mz = @SMatrix  [  1 -Az a[1] 
#                 Az 1 a[2] 
#                 -Ay Ax a[3]
#                 ]
        
#         up[i,:] .= [det(Mx) / det(M), det(My) / det(M), det(Mz) / det(M)]
#         end
        
#         # # and the final step:
#     f(Acache, .5 * (up + u), t + dt / 2) 
#     Acache += W.dW * √(2α * kbT)

#     Acache *= dt / 2 

#     @inbounds for i in 1:N
#         Ax, Ay, Az = @view Acache[i,:]
#         uloc = @view u[i,:]
#         a = SVector{3}(uloc + cross(uloc, Acache[i,:]))
#         M = @SMatrix [1 -Az Ay 
#             Az 1 -Ax 
#         -Ay Ax 1]
        
#         Mx = @SMatrix [  a[1] -Az Ay 
#         a[2] 1 -Ax 
#         a[3] Ax 1
#         ]
#         My = @SMatrix  [  1  a[1] Ay 
#                 Az  a[2] -Ax 
#                 -Ay a[3] 1
#                 ]
#         Mz = @SMatrix  [  1 -Az a[1] 
#                 Az 1 a[2] 
#                 -Ay Ax a[3]
#                 ]
        
#         u[i,:] .= [det(Mx) / det(M), det(My) / det(M), det(Mz) / det(M)]
#     end
# end

