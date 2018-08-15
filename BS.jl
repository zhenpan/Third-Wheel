using SpecialFunctions

immutable Inner_binary
    m1::Float64
    m2::Float64
    a_hat::Array{Float64}
    b_hat::Array{Float64}
    S1::Array{Float64}
    S2::Array{Float64}
    L::Array{Float64}
    E::Float64
    a::Float64
    ecc::Float64
end

immutable Third_wheel
    A_hat::Array{Float64}
    B_hat::Array{Float64}
    Z_hat::Array{Float64}
    L3::Array{Float64}
end

function Inner_binary(;m1 = 1., m2 = 1., a = 1000., ecc = 0.3)
    E, Lmag = geo2phy(m1, m2, a, ecc)

    a_hat = [1., 0., 0.]
    b_hat = [0., 1., 0.]
    S1 = 0.1*(m1^2)*[0., 0., 1.]
    S2 = 0.1*(m2^2)*[0., 0., 1.]
    L  = Lmag*[0., 0., 1.]

    inbn = Inner_binary(m1, m2, a_hat, b_hat, S1, S2, L, a, ecc)
    return inbn
end

function Third_wheel(EA3::Array{Float64}, q::Float64, m1::Float64, m2::Float64, m3::Float64)
    ϕ = EA3[1]
    θ = EA3[2]
    Ψ = EA3[3]

    A_hat = [cos(ψ)*cos(ϕ)-cos(θ)*sin(ϕ)*sin(ψ), -sin(ψ)*cos(ϕ)-cos(θ)*sin(ϕ)*cos(ψ),  sin(θ)*sin(ϕ)]
    B_hat = [cos(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*sin(ψ), -sin(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*cos(ψ), -sin(θ)*cos(ϕ)]
    Z_hat = [sin(ψ)*sin(θ), cos(ψ)*sin(θ), cos(θ)]

    M12   = m1+m2
    μ3    = m3*(m1+m2)/(m1+m2+m3)
    magL3 = sqrt(q*μ3*μ3*M12)
    L3 = magL3*Z_hat

    outsd = Third_wheel(A_hat, B_hat, Z_hat, L3)
    return outsd
end

function f12(inbn::Inner_binary, outsd::Third_wheel)
    A_a = dot(outsd.A_hat, inbn.a_hat)
    A_b = dot(outsd.A_hat, inbn.b_hat)
    B_a = dot(outsd.B_hat, inbn.a_hat)
    B_b = dot(outsd.B_hat, inbn.b_hat)

    ecc = inbn.ecc
    e_1 = -4*ecc*(besselj(0, ecc)- besselj(2,ecc)) + (3-ecc^2)*( besselj(-1, ecc)- besselj(3,ecc) )
    e_2 = 0.
    e_3 = 2*ecc*(besselj(0, ecc)- besselj(2,ecc)) - (3-2*ecc^2)*( besselj(-1, ecc)- besselj(3,ecc) )
    e_4 = 6*sqrt(1-ecc^2)*( besselj(-1, ecc)+besselj(3,ecc)- ecc*( besselj(0, ecc)+ besselj(2,ecc)) )

    f1 = e_1*A_a*B_a + e_3*A_b*B_b + 0.5*(e_2 + e_4)*(A_a*B_b+ B_a*A_b)
    f2 = e_1*(A_a^2 -B_a^2) + (e_2+e_4)*(A_a*A_b-B_a*B_b) + e_3*(A_b^2-B_b^2)
    f1 = -8*f1
    f2 = 4*f2
    return f1, f2
end

function δphy(inbn::Inner_binary, outsd::Third_wheel, M0::Float64)
    f1, f2 = f12(inbn, outsd)

    K = sqrt(2*(m1+m2)/(m1+m2+m3)*q^3/a^3)

    δE  = (m1*m2*m3)/(m1+m2)*(-5*a^2/q^3)*sqrt(pi)/120.*(K^2.5)*exp(-2K/3.)*(f1*cos(M0) + f2*sin(M0))
    δS1 =
    δS2 =
    δL  = -3π*m3*a/(8*m1*m2)*sqrt(2./q^3/M123) f(L, L3)
    return y, δS1, δS2, δS3
end

function geo2phy(m1::Float64, m2::Float64, a::Float64, ecc::Float64)
    μ = m1*m2/(m1+m2)
    E = -0.5*m1*m2/a
    L = sqrt( (ecc^1-1)/(2E)*(μ*(m1*m2)^2) )
    return E, L
end

function phy2geo(m1::Float64, m2::Float64, E::Float64, L::Float64)
    μ = m1*m2/(m1+m2)
    a = m1*m2/(2*abs(E))
    ecc = sqrt( 1+2*E*L^2/(μ*(m1*m2)^2) )
    return a, ecc
end
