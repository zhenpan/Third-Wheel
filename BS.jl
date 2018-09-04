using SpecialFunctions

immutable Inner_binary
    m1::Float64
    m2::Float64
    a_hat::Array{Float64}
    b_hat::Array{Float64}
    c_hat::Array{Float64}
    S1::Array{Float64}
    S2::Array{Float64}
    L_vec::Array{Float64}
    e_vec::Array{Float64}
    E::Float64
    a::Float64
    χ::Float64
end

immutable Third_wheel
    m3::Float64
    A_hat::Array{Float64}
    B_hat::Array{Float64}
    C_hat::Array{Float64}
    L3_vec::Array{Float64}
    A_vec::Array{Float64}
    E3::Float64
    q::Float64
end

function Inner_binary(;m1 = 1., m2 = 1., a = 1000., e_mag = 0.3)
    a_hat = [1., 0., 0.]
    b_hat = [0., 1., 0.]
    c_hat = [0., 0., 1.]
    S1 = 0.1*(m1^2)*[0., 0., 1.]
    S2 = 0.1*(m2^2)*[0., 0., 1.]

    E, L_mag = geo2phy(m1, m2, a, e_mag)
    L_vec = L_mag * c_hat
    e_vec = e_mag * a_hat

    χ = dot( (S1/m1 + S2/m2)/(m1+m2), c_hat)

    return Inner_binary(m1, m2, a_hat, b_hat, c_hat, S1, S2, L_vec, e_vec, E, a, χ)
end


function Third_wheel(;EA3=[0.,0.,0.], E3=0., q = 1.e4, m1 = 1., m2 = 1., m3 = 1.)
    ϕ = EA3[1]
    θ = EA3[2]
    ψ = EA3[3]

    A_hat = [cos(ψ)*cos(ϕ)-cos(θ)*sin(ϕ)*sin(ψ), -sin(ψ)*cos(ϕ)-cos(θ)*sin(ϕ)*cos(ψ),  sin(θ)*sin(ϕ)]
    B_hat = [cos(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*sin(ψ), -sin(ψ)*sin(ϕ)+cos(θ)*cos(ϕ)*cos(ψ), -sin(θ)*cos(ϕ)]
    C_hat = [sin(ψ)*sin(θ), cos(ψ)*sin(θ), cos(θ)]

    M12   = m1+m2
    M123  = m1+m2+m3
    μ3    = m3*M12/M123
    L3_vec= sqrt(2q*μ3*m3*M12)*C_hat
    A_vec = μ3 * M123 * A_hat

    return Third_wheel(m3, A_hat, B_hat, C_hat, L3_vec, A_vec, E3, q)
end

function Inbn_updater!(inbn::Inner_binary, outsd::Third_wheel, M0::Float64)
    m1  = inbn.m1
    m2  = inbn.m2
    m3  = outsd.m3
    a   = inbn.a
    q   = outsd.q

    M12 = m1+m2
    M123= m1+m2+m3

    f1_α, f2_α, f1_β, f2_β = f12(inbn, outsd)
    K   = sqrt( 2*(M12/M123)*(q^3/a^3) )
    δE  = (m1*m2*m3)/(m1+m2)*(a^2/q^3)*(-sqrt(π)/24. )*(  K^2.5 )*exp(-2K/3.)*(f1_α*cos(M0)  + f2_α*sin(M0))
    δE += (m1*m2*m3)/(m1+m2)*(a^2/q^3)*(-sqrt(8π)/24.)*((2K)^2.5)*exp(-4K/3.)*(f1_β*cos(2M0) + f2_β*sin(2M0))

    A_hat = outsd.A_hat
    B_hat = outsd.B_hat
    C_hat = outsd.C_hat
    a_hat = inbn.a_hat
    b_hat = inbn.b_hat
    c_hat = inbn.c_hat

    A_a  = dot(A_hat, a_hat)
    A_b  = dot(A_hat, b_hat)
    A_c  = dot(A_hat, c_hat)
    B_a  = dot(B_hat, a_hat)
    B_b  = dot(B_hat, b_hat)
    B_c  = dot(B_hat, c_hat)

    e_vec = inbn.e_vec
    e_mag = norm(e_vec)
    e_A   = dot(e_vec, A_hat)
    e_B   = dot(e_vec, B_hat)
    e_C   = dot(e_vec, C_hat)

    𝚥_mag = sqrt(1-e_mag^2)
    𝚥_vec = 𝚥_mag*inbn.c_hat
    𝚥_A   = dot(𝚥_vec, A_hat)
    𝚥_B   = dot(𝚥_vec, B_hat)
    𝚥_C   = dot(𝚥_vec, C_hat)

    scale = (3π/4)*(a/q)^1.5*m3/sqrt(2*M12*M123)
    δe1   = scale * cross(𝚥_A*A_hat+𝚥_B*B_hat-4𝚥_vec+5(A_a^2 + B_a^2)*𝚥_vec, e_vec)
    e_vec = rot_app(e_vec, δe1)
    e_vec = ( 1.- scale * (5𝚥_mag)*(A_a*A_b+B_a*B_b) )*e_vec

    δ𝚥1   = scale*( cross(𝚥_A*A_hat+𝚥_B*B_hat, 𝚥_vec) -5e_mag^2*(A_a*A_c+B_a*B_c)*b_hat )
    𝚥_vec = rot_app(𝚥_vec, δ𝚥1)
    𝚥_vec = ( 1. + scale*(5e_mag^2)*(A_a*A_b + B_a*B_b)/𝚥_mag )*𝚥_vec

    e_mag = norm(e_vec)                                   # dot(e, 𝚥) = 0, j=√1-e^2
    𝚥_vec = 𝚥_vec-dot(𝚥_vec, e_vec)/e_mag^2*e_vec
    𝚥_vec = sqrt(1-e_mag^2)/norm(𝚥_vec) * 𝚥_vec

    a_hat = normalize(e_vec)
    c_hat = normalize(𝚥_vec)
    b_hat = cross(c_hat, a_hat)

    E = inbn.E + δE
    a = m1*m2/(2*abs(E))
    L_mag = m1*m2*sqrt(a/M12)*sqrt(1-e_mag^2)
    L_vec = L_mag*c_hat

    δS1 = [0., 0., 0.]
    δS2 = [0., 0., 0.]
    S1  = inbn.S1 + δS1
    S2  = inbn.S2 + δS2
    χ   = dot( (S1/m1 + S2/m2)/(m1+m2), c_hat)

    δL3_vec = -(L_vec - inbn.L_vec)
    δE3     = -δE

    return Inner_binary(m1, m2, a_hat, b_hat, c_hat, S1, S2, L_vec, e_vec, E, a, χ), δL3_vec, δE3
end



function Outsd_updater!(outsd::Third_wheel, inbn::Inner_binary, δL3_vec::Array{Float64}, δE3::Float64)
    e_vec = inbn.e_vec
    e_mag = norm(e_vec)
    e_A   = dot(e_vec, outsd.A_hat)
    e_B   = dot(e_vec, outsd.B_hat)

    𝚥_vec = sqrt(1-e_mag^2)*inbn.c_hat
    𝚥_A   = dot(𝚥_vec, outsd.A_hat)
    𝚥_B   = dot(𝚥_vec, outsd.B_hat)

    m1 = inbn.m1
    m2 = inbn.m2
    m3 = outsd.m3
    a  = inbn.a
    q  = outsd.q

    scale = -1.5*m1*m2*m3/(m1+m2)*(a^2/q^2)
    δtmp  = (π/4)*(𝚥_B*𝚥_vec-5e_B*e_vec)
    δA1   = (π/16)*( (5𝚥_B^2-𝚥_A^2)- 5*(5e_B^2-e_A^2))*outsd.B_hat + dot(δtmp, outsd.B_hat)*outsd.B_hat + dot(δtmp, outsd.C_hat)*outsd.C_hat
    δA1   = scale*δA1
    δA_mag= scale*( (3π/8)*(𝚥_A*𝚥_B-5e_A*e_B) + dot(δtmp, outsd.A_hat) )
    A_vec = rot_app(outsd.A_vec, δA1)
    A_vec = (1 + δA_mag/norm(A_vec))*A_vec

    L3_vec = outsd.L3_vec + δL3_vec
    L3_mag = norm(L3_vec)
    A_vec  = A_vec - dot(A_vec, L3_vec)/L3_mag^2*L3_vec
    A_vec  = m3*(m1+m2)/norm(A_vec)*A_vec

    A_hat  = normalize(A_vec)
    C_hat  = normalize(L3_vec)
    B_hat  = cross(C_hat, A_hat)

    E3     = outsd.E3 + δE3
    μ3     = m3*(m1+m2)/(m1+m2+m3)
    q      = 0.5*L3_mag^2/(μ3*m3*(m1+m2))

    return Third_wheel(m3, A_hat, B_hat, C_hat, L3_vec, A_vec, E3, q)
end

function f12(inbn::Inner_binary, outsd::Third_wheel)
    A_a = dot(outsd.A_hat, inbn.a_hat)
    A_b = dot(outsd.A_hat, inbn.b_hat)
    B_a = dot(outsd.B_hat, inbn.a_hat)
    B_b = dot(outsd.B_hat, inbn.b_hat)

    ecc = norm(inbn.e_vec)
    α_1 = -4*ecc*(besselj(0, ecc)- besselj(2,ecc)) + (3-ecc^2)*( besselj(-1, ecc)- besselj(3,ecc) )
    β_1 = -4*ecc*(besselj(1, 2ecc)-besselj(3,2ecc))+ (3-ecc^2)*( besselj(0, 2ecc)- besselj(4,2ecc) )
    α_2 = 0.
    β_2 = 0.
    α_3 = 2*ecc*(besselj(0, ecc)- besselj(2,ecc)) - (3-2*ecc^2)*( besselj(-1, ecc)- besselj(3,ecc) )
    β_3 = 2*ecc*(besselj(1, 2ecc)-besselj(3,2ecc))- (3-2*ecc^2)*( besselj(0, 2ecc)- besselj(4,2ecc) )
    α_4 = 6*sqrt(1-ecc^2)*( besselj(-1, ecc)+besselj(3, ecc)- ecc*( besselj(0,  ecc)+ besselj(2,ecc)) )
    β_4 = 6*sqrt(1-ecc^2)*( besselj(0, 2ecc)+besselj(4,2ecc)- ecc*( besselj(1, 2ecc)+ besselj(3,2ecc)) )

    f1_α = -8*( α_1*A_a*B_a + α_3*A_b*B_b + 0.5*(α_2 + α_4)*(A_a*B_b+ B_a*A_b) )
    f2_α =  4*( α_1*(A_a^2 -B_a^2) + (α_2+α_4)*(A_a*A_b-B_a*B_b) + α_3*(A_b^2-B_b^2) )

    f1_β = -8*( β_1*A_a*B_a + β_3*A_b*B_b + 0.5*(β_2 + β_4)*(A_a*B_b+ B_a*A_b) )
    f2_β =  4*( β_1*(A_a^2 -B_a^2) + (β_2+β_4)*(A_a*A_b-B_a*B_b) + β_3*(A_b^2-B_b^2) )

    return f1_α, f2_α, f2_β, f2_β
end

function rot_app(Vec::Array{Float64}, δVec::Array{Float64})
    return norm(Vec)/norm(Vec+δVec)*(Vec+δVec)
end

function geo2phy(m1::Float64, m2::Float64, a::Float64, ecc::Float64)
    μ = m1*m2/(m1+m2)
    E = -0.5*m1*m2/a
    L = sqrt((ecc^2-1)*μ/(2E)) *(m1*m2)
    return E, L
end

function phy2geo(m1::Float64, m2::Float64, E::Float64, L::Float64)
    μ = m1*m2/(m1+m2)
    a = m1*m2/(2*abs(E))
    ecc = sqrt( 1+2*E*L^2/(μ*(m1*m2)^2) )
    return a, ecc
end
