include("BS.jl")
include("sampling.jl")

inbn  = Inner_binary(m1 = 10., m2 = 10., a = 1.e8, e_mag = 1.e-5)

N_smp = 5000
m3    = 1.
E_σ   =  0.5*(m3*(inbn.m1+inbn.m2)/(m3+inbn.m1+inbn.m2))*(1.e-4)^2
q_min =  2.*inbn.a
q_max =  20.*inbn.a
EA3_smp, E3_smp, q_smp = Samp(N_smp, E_σ, q_min, q_max)

N_cap = 0

for i = 1:N_smp
    outsd = Third_wheel(EA3= EA3_smp[i,:], E3=E3_smp[i], q = q_smp[i], m1 = inbn.m1, m2 = inbn.m2, m3 = m3)
    inbn, cap  = scatter!(i, inbn, outsd)
    N_cap += cap
end

function scatter!(ith, inbn, outsd)
    M0   = 2π*rand()
    inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
    E3   = outsd.E3 + δE3

    Nsct = 1
    Npos = δE3 > 0. ? 1:0
    Nneg = δE3 > 0. ? 0:1
    println(@sprintf "(ith,Nsct,Np)=(%4i, %3i, %3i), E3=%.2e, δE3=%.2e, E12=%.2e, M0/2π=%.2f, q/a=%.4e, ecc=%.4e, a=%.3e, χ=%.3f" ith Nsct Npos E3 δE3 inbn.E M0/(2π) outsd.q/inbn.a norm(inbn.e_vec) inbn.a inbn.χ)

    while (E3 < 0.) & (Nsct < 200)                               
        outsd = Outsd_updater!(outsd, inbn, δL3_vec, δE3)
        M0    = 2π*rand()
        inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
        E3 += δE3

        Nsct += 1
        Npos += δE3 > 0. ? 1:0
        Nneg += δE3 > 0. ? 0:1
        println(@sprintf "(ith,Nsct,Npos)= (%4i %3i, %3i), E3= %.2e, δE3=%.2e, E12=%.2e, M0/2π=%.2f, q/a=%.4e, ecc=%.4e, a=%.3e,χ=%.3f" ith Nsct Npos E3 δE3 inbn.E M0/(2π) outsd.q/inbn.a norm(inbn.e_vec) inbn.a inbn.χ)
    end

    cap = Nsct >= 200 ? 1:0
    return inbn, cap
end
