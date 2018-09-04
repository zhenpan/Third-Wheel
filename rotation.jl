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
