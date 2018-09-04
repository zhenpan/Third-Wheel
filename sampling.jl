using SpecialFunctions

function Samp(N_smp::Int, E_σ::Float64, q_min::Float64, q_max::Float64)
    ϕ = 2π*rand(N_smp)
    θ = acos.(2*rand(N_smp)-1.)
    ψ = 2π*rand(N_smp)
    EA3_smp = hcat(ϕ, θ, ψ)

    U_uni  = rand(N_smp)
    E3_smp = zeros(N_smp)
    for i = 1:N_smp
        E3_smp[i] = (erfinv(U_uni[i]))^2 * E_σ
    end

    U_uni  = rand(N_smp)
    q_smp  = q_min + (q_max-q_min)*U_uni
    return EA3_smp, E3_smp, q_smp
end
