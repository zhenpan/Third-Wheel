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

function scatter!(ith, inbn, outsd)
    M0   = 2π*rand()
    inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
    E3   = outsd.E3 + δE3

    Nsct = 1
    Npos = δE3 > 0. ? 1:0
    Nneg = δE3 > 0. ? 0:1
    println(@sprintf "(ith,Ns,Np)=(%4i, %3i, %3i), E3=%.2e, δE3=%.2e, E12=%.2e, M0/2π=%.2f, q/a=%.4e, e=%.4e, e3=%.4e, a=%.3e, χ=%.3f" ith Nsct Npos outsd.E3 δE3 inbn.E M0/(2π) outsd.q/inbn.a norm(inbn.e_vec) outsd.e3 inbn.a inbn.χ)

    while (E3 < 0.) & (Nsct < 200)
        outsd = Outsd_updater!(outsd, inbn, δL3_vec, δE3)
        M0    = 2π*rand()
        inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
        E3 += δE3

        Nsct += 1
        Npos += δE3 > 0. ? 1:0
        Nneg += δE3 > 0. ? 0:1
        println(@sprintf "(ith,Ns,Np)= (%4i %3i, %3i), E3= %.2e, δE3=%.2e, E12=%.2e, M0/2π=%.2f, q/a=%.4e, e=%.4e, e3=%.4e, a=%.3e,χ=%.3f" ith Nsct Npos outsd.E3 δE3 inbn.E M0/(2π) outsd.q/inbn.a norm(inbn.e_vec) outsd.e3 inbn.a inbn.χ)
    end

    cap = Nsct >= 200 ? 1:0
    return inbn, Nsct, cap
end


function fake_scatter!(ith, inbn, outsd)
    M0   = 2π*rand()
    inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
    E3   = outsd.E3 + δE3

    Nsct = 1
    Npos = δE3 > 0. ? 1:0
    Nneg = δE3 > 0. ? 0:1
    println(@sprintf "(ith,Ns,Np)=(%4i, %3i, %3i), E3=%.2e, δE3=%.2e, E12=%.2e, M0/2π=%.2f, q/a=%.4e, e=%.4e, e3=%.4e, a=%.3e, χ=%.3f" ith Nsct Npos E3 δE3 inbn.E M0/(2π) outsd.q/inbn.a norm(inbn.e_vec) outsd.e3 inbn.a inbn.χ)

    while (E3 < 0.) & (Nsct < 200)
        #outsd = Outsd_updater!(outsd, inbn, δL3_vec, δE3)
        M0    = 2π*rand()
        inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
        E3 += δE3

        Nsct += 1
        Npos += δE3 > 0. ? 1:0
        Nneg += δE3 > 0. ? 0:1
        println(@sprintf "(ith,Ns,Np)= (%4i %3i, %3i), E3= %.2e, δE3=%.2e, E12=%.2e, M0/2π=%.2f, q/a=%.4e, e=%.4e, e3=%.4e, a=%.3e,χ=%.3f" ith Nsct Npos E3 δE3 inbn.E M0/(2π) outsd.q/inbn.a norm(inbn.e_vec) outsd.e3 inbn.a inbn.χ)
    end

    cap = Nsct >= 200 ? 1:0
    return inbn, cap
end
