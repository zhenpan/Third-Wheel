function Samp(N_smp::Int)
    ϕ = 2π*rand(N_smp)
    θ = acos.(2*rand(N_smp)-1.)
    ψ = 2π*rand(N_smp)

    A_hat = [cos.(ψ).*cos.(ϕ)-cos.(θ).*sin.(ϕ).*sin.(ψ), -sin.(ψ).*cos.(ϕ)-cos.(θ).*sin.(ϕ).*cos.(ψ),  sin.(θ).*sin.(ϕ)]
    B_hat = [cos.(ψ).*sin.(ϕ)+cos.(θ).*cos.(ϕ).*sin.(ψ), -sin.(ψ).*sin.(ϕ)+cos.(θ).*cos.(ϕ).*cos.(ψ), -sin.(θ).*cos.(ϕ)]
    C_hat = [sin.(ψ).*sin.(θ), cos.(ψ).*sin.(θ), cos.(θ)]

    return A_hat, B_hat, C_hat
end

A_hat, B_hat, C_hat = Samp(1000)

tt = A_hat[1].*A_hat[2] + B_hat[1].*B_hat[2]
histogram(tt,nbins=20)

function rand_walk()
    Ns = 1
    E3 = cos(2π*rand())
    while (E3 < 0.) & (Ns < 1000)
        E3 += cos(2π*rand())
        Ns +=1
    end
    return Ns
end

N_smp = 500
N_sct = zeros(N_smp)

for i = 1:N_smp
    N_sct[i] = rand_walk()
end
