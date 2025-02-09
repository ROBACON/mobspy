from mobspy import *

Promoter, Start_Positions, Tet, Mortal = BaseSpecies()
Promoter.inactive, Promoter.active
Ribo, RNA_Poly = New(Start_Positions)
P2, Ptet = New(Promoter)
Mrna_P2, Mrna_Ptet, GFP, RFP, CFP, GFP_F_RFP = New(Mortal)

# Death here does not consider a compound for degradation.
Mortal >> Zero [1]

# Promoter activation - only Ptet is activated. P2 is always active
Rev[Ptet.inactive + Tet >> Ptet.active][1, 1]

# Read expresses the reading of RNA and DNA - both transcription and translation
# For translation the RNA Polymerase is the reader, for protein it is the Ribossome
def Read(Pro, R, strand):
    sp = 'started_' + strand[0][0] # sp stands for start position
    Start_Positions.c(sp)
    rate = [lambda r1, r2: 2 if r1.active else 1, 1]
    # From free position to bounded to a site
    Rev[Pro + R.c('free_' + str(R)) >> Pro + R.c(sp).c('at_' + strand[0][0])][rate]
    next_location = strand[1:]
    # Movement of the reader
    for (location, Product), (next_l, _) in zip(strand, next_location):
        R.c(sp).c('at_' + location) >> R.c(sp).c('at_' + next_l) + Product [1]
    # Remove reader from final location
    R.c(sp).c('at_' + next_location[-1][0]) >> R.c(sp).c('free_' + str(R)) [1]

# Zero produces nothing
Read(Ptet, RNA_Poly, [('ptet_dna', Mrna_Ptet), ('ptet_end', Zero)])
Read(P2, RNA_Poly, [('p2_dna', Mrna_P2), ('p2_end', Zero)])
Read(Mrna_Ptet, Ribo, [('gfp', GFP), ('rfp', GFP_F_RFP), ('rfp_end', Zero)])
Read(Mrna_P2, Ribo, [('cfp',  CFP), ('cfp_end',  Zero)])
Read(Mrna_Ptet, Ribo, [('rfp', RFP), ('rfp_end', Zero)])

# P2 is set to active as it is a constitutive promoter
model = set_counts({RNA_Poly: 100, Ribo: 100, GFP: 0, RFP: 0, CFP: 0, Ptet: 1, P2.active: 1,
                    Tet: 100, Mrna_Ptet: 0, Mrna_P2: 0, GFP_F_RFP: 0})
S = Simulation(model)
print(S.compile())










