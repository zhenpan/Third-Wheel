function Samp(N_smp::Int)
    ϕ = 2π*rand(N_smp)
    θ = acos.(2*rand(N_smp)-1.)
    ψ = 2π*rand(N_smp)

    A_hat = [cos.(ψ).*cos.(ϕ)-cos.(θ).*sin.(ϕ).*sin.(ψ), -sin.(ψ).*cos.(ϕ)-cos.(θ).*sin.(ϕ).*cos.(ψ),  sin.(θ).*sin.(ϕ)]
    B_hat = [cos.(ψ).*sin.(ϕ)+cos.(θ).*cos.(ϕ).*sin.(ψ), -sin.(ψ).*sin.(ϕ)+cos.(θ).*cos.(ϕ).*cos.(ψ), -sin.(θ).*cos.(ϕ)]
    C_hat = [sin.(ψ).*sin.(θ), cos.(ψ).*sin.(θ), cos.(θ)]

    return A_hat, B_hat, C_hat
end


using PyPlot

A_hat, B_hat, C_hat = Samp(1000)

tt = A_hat[1].*A_hat[2] + B_hat[1].*B_hat[2]
histogram(tt,nbins=20)

function rand_cos(Nlim)
    Ns = 1
    E3 = cos(2π*rand())
    while (E3 < 0.) & (Ns < Nlim)
        E3 += cos(2π*rand())
        Ns +=1
    end
    return Ns
end

function rand_uni()
    Ns = 1
    E3 = rand([-1,1])
    while (E3 < 0.) & (Ns < 5000)
        E3 += rand([-1,1])
        Ns +=1
    end
    return Ns
end

N_smp = 500
N_sct = zeros(N_smp)

for i = 1:N_smp
    N_sct[i] = rand_cos(100000)
end

plot(N_sct)


function rand_walk(Nlim)
    j_pls = zeros(Nlim+2)
    j_mns = zeros(Nlim+2)

    j_pls[1] = 0.5
    j_mns[1] = 0.
    leak     = 0.5

    for t = 1:Nlim
        jt_pls   = zeros(Nlim+2)
        jt_mns   = zeros(Nlim+2)

        jt_pls[1] = 0.
        jt_mns[1] = 0.5*(j_mns[2] + j_pls[1])
        for i = 2:Nlim+1
            jt_pls[i] = 0.5*(j_pls[i-1] + j_mns[i])
            jt_mns[i] = 0.5*(j_mns[i+1] + j_pls[i])
        end

        j_pls = deepcopy(jt_pls)
        j_mns = deepcopy(jt_mns)
        leak += j_mns[1]
    end
    return leak, sum(j_pls[2:end]+j_mns[2:end])
end

Nlim = [1e2, 1e3, 1e4, 1e5]
Capt = zeros(4)

for i = 1:4
    leak, Capt[i] = rand_walk(Int64(Nlim[i]))
end
