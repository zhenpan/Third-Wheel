using SpecialFunctions

function Samp(N_smp::Int, E_σ::Float64, q_min::Float64, q_max::Float64)
    U_uni  = rand(N_smp)
    E3_smp = zeros(N_smp)
    for i = 1:N_smp
        E3_smp[i] = (erfinv(U_uni[i]))^2 * E_σ
    end

    U_uni  = rand(N_smp)
    q_smp  = q_min + (q_max-q_min)*U_uni
    return E3_smp, q_smp
end
