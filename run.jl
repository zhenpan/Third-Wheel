include("BS.jl")

inbn  = Inner_binary(m1 = 1., m2 = 1., a = 1000., e_mag = 1.e-3) #m1 is fixed as unity, m2<=m1, no constraint on m3


for iwhl = 1:10
    EA3, q, m3 =  2π*rand(3), 1.e4, 1.        #Euler Angles of the 3rd body, q and m3
    outsd = Third_wheel(EA3, q, inbn.m1, inbn.m2, m3)

    M0   = 2π*rand()
    inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)

    isct = 1

    println(@sprintf "(iwhl, isct)= (%2i, %2i), E3= %.2e, δE3=%.2e, q/a=%.2e, ecc=%.4f"  iwhl isct outsd.E3+δE3 δE3 outsd.q/inbn.a norm(inbn.e_vec))
    while (outsd.E3 + δE3 < 0.) & (isct < 100)                ##Third boday captured, multiple scatterings
        isct += 1
        outsd = Outsd_updater!(outsd, inbn, δL3_vec, δE3)
        M0    = 2π*rand()
        inbn, δL3_vec, δE3 = Inbn_updater!(inbn, outsd, M0)
        println(@sprintf "(iwhl, isct)= (%2i, %2i), E3= %.2e, δE3=%.2e, q/a=%.2e, ecc=%.4f"  iwhl isct outsd.E3+δE3 δE3 outsd.q/inbn.a norm(inbn.e_vec))
        end
end
