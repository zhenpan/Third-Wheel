include("BS.jl")
include("sampling.jl")
using PyPlot

inbn  = Inner_binary(m1 = 10., m2 = 10., a = 1.e8, e_mag = 0.1)

N_smp = 5000
m3    = 1.
E_σ   =  0.5*(m3*(inbn.m1+inbn.m2)/(m3+inbn.m1+inbn.m2))*(1.e-4)^2
q_min =  5.*inbn.a
q_max =  10.*inbn.a
EA3_smp, E3_smp, q_smp = Samp(N_smp, E_σ, q_min, q_max)

Nsct= zeros(N_smp)
Ein = zeros(N_smp)
ein = zeros(N_smp)

N_cap = 0

for i = 1:N_smp
    outsd = Third_wheel(EA3= EA3_smp[i,:], E3=E3_smp[i], q = q_smp[i], m1 = inbn.m1, m2 = inbn.m2, m3 = m3)
    inbn, Ns, cap  = scatter!(i, inbn, outsd, Nlim = 1000)

    N_cap += cap
    Nsct[i]= Ns
    Ein[i] = inbn.E
    ein[i] = norm(inbn.e_vec)
end


plot(ein)
plot(Nsct/1000.)
