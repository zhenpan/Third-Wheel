include("BS.jl")
include("sampling.jl")
using PyPlot

inbn  = Inner_binary(m1 = 10., m2 = 10., a = 1.e8, e_mag = 0.1)

N_smp = 5000
m3    = 1.
E_σ   =  0.5*(m3*(inbn.m1+inbn.m2)/(m3+inbn.m1+inbn.m2))*(1.e-4)^2
q_min =  3.*inbn.a
q_max =  10.*inbn.a
EA3_smp, E3_smp, q_smp = Samp(N_smp, E_σ, q_min, q_max)

Nsct= zeros(N_smp)
Ein = zeros(N_smp)
ein = zeros(N_smp)
χin = zeros(N_smp)
e3  = zeros(N_smp)

N_cap = 0

for i = 1:N_smp
    outsd = Third_wheel(EA3= EA3_smp[i,:], E3=E3_smp[i], q = q_smp[i], m1 = inbn.m1, m2 = inbn.m2, m3 = m3)
    inbn, Ns, cap  = scatter!(i, inbn, outsd, Nlim = 1000)

    N_cap += cap
    Nsct[i]= Ns
    Ein[i] = inbn.E
    ein[i] = norm(inbn.e_vec)
    χin[i] = inbn.χ
    e3[i]  = outsd.e3

end


plot(ein, lw=2, label=L"$e$")
plot(χin, lw=2, label=L"$χ$")
#plot(e3, lw=2, label=L"$e_3$")
plot(Nsct/1000., lw =0.6, "--", label=L"$N_{\rm cp}/10^3$")
legend()
