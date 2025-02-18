from mobspy import *

# We start with all base species
Link1, Link2, Phos_2 = BaseSpecies()
Link1.nl_1, Link1.l_1, Link2.nl_2, Link2.l_2
Phos_2.not_phos, Phos_2.phos_1, Phos_2.phos_2

# A and B can be linked
A = Link1 * Link2
B = New(Link1)
# C can be phosporolized twice
C = Phos_2 * Link2

# Rev stands for reversible reaction
Rev[A.nl_1 + B.nl_1 >> A.l_1 + B.l_1][1e-4, lambda r: 0.1 * r]

# A binds to C
A.nl_2.l_1 + C.nl_2.not_phos >> A.l_2.l_1 + C.l_2.not_phos[1e-4]

# C is phosphoralized releases A
C.l_2.not_phos + A.l_2.l_1 >> C.nl_2.phos_1 + A.nl_2.l_1[lambda r: 1 * r]

# C phosporalized binds to unbound A
A.l_2.nl_1 + C.nl_2.phos_1 >> A.l_2.nl_1 + C.l_2.phos_1[1e-4]

# Final C site is modified
A.l_2.nl_1 + C.l_2.phos_1 >> A.nl_2.nl_1 + C.nl_2.phos_2[lambda r: r]

# initial conditions
A(1000), B(1000), C(10000)
S = Simulation(A | B | C)
print(S.compile())
