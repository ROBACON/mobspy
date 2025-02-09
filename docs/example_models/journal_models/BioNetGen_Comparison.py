from mobspy import *

R0, L0, A0, kon, koff, kAon, kAoff, kAp, kAdp = ModelParameters(100, 500, 100, 0.01,
                                                                0.1, 0.01, 0.1, 0.01, 0.1)

r_link = BaseSpecies()
r_link.r_0, r_link.r_1
R, L, A = New(r_link)

R.r_0 + L.r_0 >> R.r_1 + L.r_1 [kon, koff]

R.a_0 + A.r_0 >> R.a_1 + A.r_1 [kAon, kAoff]

(R.a_1+ L + A.n_p).r_1 >> (R.a_1 + L + A.y_p).r_1 [lambda R: kAp*R**3]

A.y_p >> A.n_p [kAdp]

R(R0), L(L0), A(A0)
S = Simulation(R, L, A)
print(S.compile())


