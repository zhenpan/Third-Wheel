include("BS.jl")
include("sampling.jl")

inbn  = Inner_binary(m1 = 10., m2 = 10., a = 1.e8, e_mag = 0.1)

N_smp = 10000
m3    = 1.
E_σ   =  0.5*(m3*(inbn.m1+inbn.m2)/(m3+inbn.m1+inbn.m2))*(1.e-4)^2
q_min =  10.*inbn.a
q_max =  11.*inbn.a
E3_smp, q_smp = Samp(N_smp, E_σ, q_min, q_max)

for i = 1:N_smp
    println(@sprintf "%i-th scattering" i)
    outsd = Third_wheel(EA3= 2π*rand(3), E3=E3_smp[i], q = q_smp[i], m1 = inbn.m1, m2 = inbn.m2, m3 = m3)
    inbn  = scatter!(inbn, outsd)
end

function scatter!(inbn, outsd)
    M0   = 2π*rand()
    inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
    E3   = outsd.E3 + δE3

    Nsct = 1
    Npos = δE3 > 0. ? 1:0
    Nneg = δE3 > 0. ? 0:1
    println(@sprintf "(Nsct, Npos)= (%3i, %3i), E3= %.2e, δE3=%.2e, E12=%.2e, q/a=%.4e, ecc=%.4e, a=%.3e, χ=%.3f" Nsct Npos E3 δE3 inbn.E outsd.q/inbn.a norm(inbn.e_vec) inbn.a inbn.χ)

    while (E3 < 0.) & (Nsct < 200)                               ##Third boday captured, multiple scatterings
        outsd = Outsd_updater!(outsd, inbn, δL3_vec, δE3)
        M0    = 2π*rand()
        inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
        E3 += δE3

        Nsct += 1
        Npos += δE3 > 0. ? 1:0
        Nneg += δE3 > 0. ? 0:1
        println(@sprintf "(Nsct, Npos)= (%3i, %3i), E3= %.2e, δE3=%.2e, E12=%.2e,q/a=%.4e, ecc=%.4e, a=%.3e,χ=%.3f" Nsct Npos E3 δE3 inbn.E outsd.q/inbn.a norm(inbn.e_vec) inbn.a inbn.χ)
    end

    return inbn
end
