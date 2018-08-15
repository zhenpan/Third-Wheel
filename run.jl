include("BS.jl")

m1 = 1.
m2 = 1.
m3 = 1.   #m1 is fixed as unity, m2<=m1, no constraint on m3

EA3  = 2*π*rand(3)          #Euler Angles of the 3rd body

inbn = Inner_binary(m1, m2, a = 1000., ecc = 0.3)
outsd = Third_wheel(m1, m2, m3, EA3, q  = 1.e4)
δE, δS1, δS2, δL = δphy(inbn, outsd, M0)
