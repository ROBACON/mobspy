"""Tests for all example models in docs/example_models/.

Each test compiles (and where necessary runs) an example model to verify
that the DSL, compiler, and simulation pipeline work end-to-end.
Tests that run simulations use save_data=False, plot_data=False.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from mobspy import *
from mobspy.modules.meta_class import ListSpecies

# ---------------------------------------------------------------------------
# Application models
# ---------------------------------------------------------------------------


class TestAB2CD:
    def test_compile(self):
        A, B, C, D = BaseSpecies()
        A(200) + B(100) >> 2 * C + D[0.1]
        S = Simulation(A | B | C | D)
        S.duration = 5
        S.volume = 10
        result = S.compile()
        assert result is not None

    def test_run(self):
        A, B, C, D = BaseSpecies()
        A(200) + B(100) >> 2 * C + D[0.1]
        S = Simulation(A | B | C | D)
        S.save_data = False
        S.plot_data = False
        S.duration = 5
        S.volume = 10
        S.simulation_method = "stochastic"
        S.run()
        assert S.results is not None


class TestSimpleInfection:
    def test_compile(self):
        Age, Mortal, Infectable, Virus = BaseSpecies()
        Reproducer = New(Age)
        V1, V2 = New(Virus)

        Age.young >> Age.old[1]
        Reproducer >> 2 * Reproducer.young[0.1]

        def infection_rate(r1, r2):
            factor = 0.01
            factor = 2 * factor if r1.old else 1 * factor
            factor = 2 * factor if r2.is_a(V2) else 1 * factor
            return factor

        Infectable.not_infected + Virus >> Infectable.infected[infection_rate]
        Mortal >> Zero[lambda r1: 2 if r1.infected else 0.01]
        Cell = Infectable * Mortal * Age * Reproducer

        Cell(100)
        V1(20), V2(25)
        S = Simulation(Cell | V1 | V2)
        S.duration = 5
        result = S.compile()
        assert result is not None


class TestSimpleRepressor:
    def test_run(self):
        Chemical, Promoter, Protein = BaseSpecies()
        Rev[Chemical + Promoter.inactive >> Promoter.active][2, 1]
        Promoter >> Promoter + Protein[lambda promoter: 1 if promoter.inactive else 0]
        Protein >> Zero[0.1]
        Promoter.inactive(100), Chemical(1000)

        S = Simulation(Chemical | Promoter | Protein)
        S.save_data = False
        S.plot_data = False
        S.duration = 200
        S.run()
        assert S.results is not None


class TestDonorReceptor:
    def test_compile(self):
        Resource, Mortal, Infectible = BaseSpecies()
        AA, Glu = New(Resource)
        AA(100), Glu(100)

        Mortal >> Zero[0.1]
        Donor, Phage = New(Mortal)
        Donor(100)

        dup_rate = lambda _, resource: 0.2 if resource.is_a(AA) else 0.1
        Donor + Resource >> 2 * Donor[dup_rate]
        Donor + Resource >> Donor + Resource + Phage[0.1]
        Infectible.low_inf >> Infectible.high_inf[0.1]
        Receptor = Mortal * Infectible
        inf_rate = lambda receptor: 0.2 if receptor.high_inf else 0.1
        Receptor.not_infected + Phage >> Receptor.early_infection[inf_rate]
        Receptor.early_infection >> Receptor.late_infection[0.1]
        Receptor + Resource >> Receptor.low_inf + Receptor[dup_rate]

        S = Simulation(Donor | Receptor | Phage | AA | Glu)
        result = S.compile()
        assert result is not None


class TestOscillator:
    def test_run(self):
        Promoter, Chemical = BaseSpecies()
        Promoter.inactive, Promoter.active
        DNAPromoter = New(Promoter)
        DNAPromoter.PLcl(100), DNAPromoter.PLacl(100), DNAPromoter.PTetR(100)
        Chemical.TetR(1)

        list_of_chemicals = ["TetR", "Lcl", "Lacl"]
        list_of_promoters = ["PLcl", "PLacl", "PTetR"]
        for che, pro in zip(list_of_chemicals, list_of_promoters, strict=False):
            Rev[
                DNAPromoter.inactive.c(pro) + Chemical.c(che)
                >> DNAPromoter.active.c(pro)
            ][1, 1]
        Chemical >> Zero[1]

        repressed = ["TetR", "Lcl", "Lacl"]
        repressors = ["Lacl", "TetR", "Lcl"]
        hill = lambda che: f"10/(1 + ({che})^3)"
        for rpsor, rpsed in zip(repressed, repressors, strict=False):
            Chemical.c(rpsor) >> Chemical.c(rpsed) + Chemical.c(rpsor)[hill]

        S = Simulation(DNAPromoter | Chemical)
        S.duration = 300
        S.plot_data = False
        S.save_data = False
        S.run()
        assert S.results is not None


class TestCRISPROscillator:
    def test_compile(self):
        Promoter, dCas, CasBinding = BaseSpecies()
        Promoter.active, Promoter.inactive, CasBinding.no_cas, CasBinding.cas

        DNAPro = New(Promoter)
        gRNA = New(CasBinding)

        G = ListSpecies(3, gRNA)
        P = ListSpecies(3, Promoter)

        rev_rt = (1.8e-3 / (u.nanomolar * u.second), 2.3e-2 / u.minute)
        Rev[gRNA.no_cas + dCas >> gRNA.cas][rev_rt]
        gRNA.no_cas >> Zero[0.0069 / u.second]

        Promoter.active >> 2 * Promoter.active[2.3e-2 / u.minute]
        Promoter >> Zero[2.3e-2 / u.minute]

        for Prom, Grna in zip(P, G, strict=False):
            act_rt = lambda dna: 5 / u.minute if dna.active else 0
            Prom >> Grna.no_cas + Prom[act_rt]

        gRNA_rep_List = [G[-1], G[0], G[1]]
        for Prom, Grna in zip(P, gRNA_rep_List, strict=False):
            dna_rt1 = 1.2e-2 * u.liter / (u.nanomoles * u.second)
            dna_rt2 = 2.3e-2 / u.minute
            Prom.active + Grna.cas >> Prom.inactive[dna_rt1]
            Prom.inactive >> 2 * Prom.active + Grna.cas[dna_rt2]

        Rev[Zero >> dCas][1 / u.minute, 2.3e-2 / u.minute]

        for i in range(3):
            P[i].active(1)
        G[0].no_cas(0), G[1].no_cas(40 * u.nanomolar), G[2].no_cas(1 * u.nanomolar)
        dCas(43 * u.nanomolar)

        S = Simulation(G | P | dCas)
        S.volume = 1 * u.femtoliter
        S.duration = 650 * u.hours
        result = S.compile()
        assert result is not None


class TestANDGate:
    def test_run(self):
        A, B, AB, Promoter, Protein = BaseSpecies()
        A + B >> AB[0.5]
        AB + Promoter.inactive >> Promoter.active[0.5]
        Promoter >> Promoter + Protein[lambda promoter: 1 if promoter.active else 0]
        Protein >> Zero[2]

        Promoter(100)
        A(50), B(50)
        S = Simulation(A | B | AB | Promoter | Protein)
        S.duration = 10
        S.level = 0
        S.save_data = False
        S.plot_data = False
        S.run()
        protein_final = S.results[Protein][0][-1]
        assert protein_final > 0


class TestNORGate:
    def test_run(self):
        Repressor, Promoter, Protein = BaseSpecies()
        Repressor + Promoter.active >> Promoter.inactive[0.5]
        A, B = New(Repressor)
        Promoter >> Promoter + Protein[lambda promoter: 1 if promoter.active else 0]
        Protein >> Zero[2]

        Promoter(100)
        A(0), B(0)
        S = Simulation(A | B | Promoter | Protein)
        S.duration = 10
        S.level = 0
        S.save_data = False
        S.plot_data = False
        S.run()
        protein_final = S.fres["Protein"][-1]
        # NOR gate with 0,0 inputs should produce protein
        assert protein_final > 0


class TestForTheTrees:
    def test_compile(self):
        Ager, Mortal, Colored, Location = BaseSpecies()
        Colored.green, Colored.yellow, Colored.brown
        Location.dense, Location.sparse
        Ager.young >> Ager.old[1 / 10 / u.year]
        Mortal >> Zero[lambda r1: 0.1 / u.year if r1.old else 0]
        Tree = Ager * Colored * Mortal * Location

        Tree.old >> Tree + Tree.green.young[0.1 / u.year]
        (
            Tree.dense.old + Tree.dense.young
            >> Tree.dense.old[1e-10 * u.decimeter**2 / u.year]
        )

        colors = ["green", "yellow", "brown"]
        for color, next_color in zip(colors, colors[1:] + colors[:1], strict=False):
            Tree.c(color) >> Tree.c(next_color)[10 / u.year]

        Tree.dense(50), Tree.dense.old(50), Tree.sparse(50), Tree.sparse.old(50)
        S = Simulation(Tree)
        S.duration = 100 * u.years
        S.volume = 1 * u.meter**2
        result = S.compile()
        assert result is not None


class TestPositivePhageFeedbackLoop:
    def test_compile(self):
        Activatable, Age, Infected, Mortal, Signal, Phage = BaseSpecies()
        Signal, Phage = New(Mortal)
        Reproducer = New(Age)

        Mortal >> Zero[0.1 / u.h]
        Signal + Activatable.inactive >> Activatable.active[20 / u.h]
        Activatable.active >> Activatable.inactive + Signal[0.01 / u.h]
        Age.young >> Age.old[1 / u.h]
        Reproducer >> 2 * Reproducer.young[0.1 / u.h]

        infection_rate = lambda r1, r2: 2 / u.h if r1.old else 1 / u.h
        Infected.not_infected + Phage >> Infected.infected[infection_rate]

        Cell = Activatable * Reproducer * Infected * Mortal
        Cell.active >> Cell.active + Signal[10 / u.h]
        Cell.active.infected >> Cell.active.infected + Phage[10 / u.h]

        initial_signal = ModelParameters(0)
        counts = {
            Cell.not_infected: 50,
            Cell.infected: 50,
            Phage: 0,
            Signal: initial_signal,
        }
        model = set_counts(counts)
        S = Simulation(model)
        S.duration = 0.3 * u.h
        result = S.compile()
        assert result is not None


class TestRandomWalk:
    def test_compile(self):
        Mesh = BaseSpecies()
        n = 3  # smaller grid for faster test
        for i in range(n):
            for j in range(n):
                coordinate = "p_" + str(i) + "_" + str(j)
                if i + 1 < n:
                    Mesh.c(coordinate) >> Mesh.c(f"p_{i + 1}_{j}")[0.1]
                if i - 1 > -1:
                    Mesh.c(coordinate) >> Mesh.c(f"p_{i - 1}_{j}")[0.1]
                if j - 1 > -1:
                    Mesh.c(coordinate) >> Mesh.c(f"p_{i}_{j - 1}")[0.1]
                if j + 1 < n:
                    Mesh.c(coordinate) >> Mesh.c(f"p_{i}_{j + 1}")[0.1]

        Bacteria, Phage = New(Mesh)
        (
            Bacteria.not_infected + Phage
            >> Bacteria.infected[lambda r1, r2: 1000000 if Mesh(r1) == Mesh(r2) else 0]
        )
        Bacteria.p_0_0(1)
        Phage.c(f"p_{n - 1}_{n - 1}")(1)

        S = Simulation(Bacteria | Phage)
        S.duration = 100
        result = S.compile()
        assert result is not None


class TestSimpleRuleBasedANDGate:
    def test_compile(self):
        A, B, C, Pa, Pb = BaseSpecies()
        C >> Zero[1]

        def Promoter_Rule(Promoter, Ligand, Protein, K):
            (
                Promoter + Ligand
                >> Promoter
                + Ligand
                + Protein[lambda p, loc: (p / u.h) * loc**4 / (K**4 + loc**4)]
            )

        Promoter_Rule(Pa, A, C, 5)
        Promoter_Rule(Pb, B, C, 8)

        initial_a = ModelParameters(0)
        initial_b = ModelParameters(0)

        model = set_counts({A: initial_a, B: initial_b, C: 0, Pa: 1, Pb: 1})
        S = Simulation(model)
        S.duration = 30 * u.h
        result = S.compile()
        assert result is not None


# ---------------------------------------------------------------------------
# Journal models
# ---------------------------------------------------------------------------


class TestAntimonyComparison:
    def test_generate_antimony(self):
        Replicator, Mortal = BaseSpecies()
        Replicator >> 2 * Replicator[lambda r: 2 * (100 - r) * r]
        Mortal >> Zero[1]

        A, B, C = New(Replicator * Mortal)
        A + B >> C[1]

        S = Simulation(A | B | C)
        antimony = S.generate_antimony()
        assert antimony is not None
        assert len(antimony) > 0


class TestBioCRNpyler1:
    def test_compile(self):
        A, B, C, D = BaseSpecies()
        A >> 2 * B[3]
        B >> C + D[1.4]
        S = Simulation(A | B | C | D)
        result = S.compile()
        assert result is not None


class TestBioCRNpyler2:
    def test_compile(self):
        Promoter, Start_Positions, Tet, Mortal = BaseSpecies()
        Promoter.inactive, Promoter.active
        Ribo, RNA_Poly = New(Start_Positions)
        P2, Ptet = New(Promoter)
        Mrna_P2, Mrna_Ptet, GFP, RFP, CFP, GFP_F_RFP = New(Mortal)

        Mortal >> Zero[1]
        Rev[Ptet.inactive + Tet >> Ptet.active][1, 1]

        def Read(Pro, R, strand):
            sp = "started_" + strand[0][0]
            Start_Positions.c(sp)
            rate = [lambda r1, r2: 2 if r1.active else 1, 1]
            Rev[Pro + R.c("free_" + str(R)) >> Pro + R.c(sp).c("at_" + strand[0][0])][
                rate
            ]
            next_location = strand[1:]
            for (location, Product), (next_l, _) in zip(
                strand, next_location, strict=False
            ):
                R.c(sp).c("at_" + location) >> R.c(sp).c("at_" + next_l) + Product[1]
            R.c(sp).c("at_" + next_location[-1][0]) >> R.c(sp).c("free_" + str(R))[1]

        Read(Ptet, RNA_Poly, [("ptet_dna", Mrna_Ptet), ("ptet_end", Zero)])
        Read(P2, RNA_Poly, [("p2_dna", Mrna_P2), ("p2_end", Zero)])
        Read(Mrna_Ptet, Ribo, [("gfp", GFP), ("rfp", GFP_F_RFP), ("rfp_end", Zero)])
        Read(Mrna_P2, Ribo, [("cfp", CFP), ("cfp_end", Zero)])
        Read(Mrna_Ptet, Ribo, [("rfp", RFP), ("rfp_end", Zero)])

        model = set_counts(
            {
                RNA_Poly: 100,
                Ribo: 100,
                GFP: 0,
                RFP: 0,
                CFP: 0,
                Ptet: 1,
                P2.active: 1,
                Tet: 100,
                Mrna_Ptet: 0,
                Mrna_P2: 0,
                GFP_F_RFP: 0,
            }
        )
        S = Simulation(model)
        result = S.compile()
        assert result is not None


class TestBioNetGenComparison:
    def test_compile(self):
        R0, L0, A0, kon, koff, kAon, kAoff, kAp, kAdp = ModelParameters(
            100, 500, 100, 0.01, 0.1, 0.01, 0.1, 0.01, 0.1
        )

        r_link, C = BaseSpecies()
        r_link.r_0, r_link.r_1
        R, L, A = New(r_link)

        R.r_0 + L.r_0 >> R.r_1 + L.r_1[kon, lambda r: koff * r]
        R.a_0 + A.r_0 >> R.a_1 + A.r_1[kAon, lambda r: kAoff * r]

        C.assign(A.r_1.n_p - R.a_1.r_0)
        C + A.n_p.r_1 >> C + A.y_p.r_1[lambda c: kAp * c]
        A.y_p >> A.n_p[kAdp]

        R(R0), L(L0), A(A0)
        S = Simulation(R | L | A | C)
        result = S.compile()
        assert result is not None


class TestBistableGutInflammation:
    def test_run(self):
        Dummy, Dilutable, Promoter, TTR, Phosphorable = BaseSpecies()
        TF, Reporter = New(Dilutable)
        TtrR, TtrS = New(Phosphorable)
        Cro, CI = New(TF)
        Pcro, PcI = New(Promoter)
        Promoter.bound, Promoter.unbound
        Phosphorable.dephospho, Phosphorable.phospho

        TtrS.dephospho + TTR >> TtrS.phospho + TTR[1 / u.min]
        TtrS.phospho + TtrR.dephospho >> TtrR.phospho + TtrS.dephospho[1 / u.min]
        TtrR.phospho >> TtrR.dephospho[5 * 1e-3 * 1 / u.second]
        TtrR.phospho >> Cro + TtrR.phospho[50 * 0.02 * 1 / u.min]

        Dummy >> CI + Dummy[30 * 0.02 * 1 / u.min]

        Dilutable >> Zero[0.02 / u.min]
        TF >> Zero[lambda r1: 0 if r1.is_a(CI) else 1.6e-2 * 1 / u.min]

        def Expression(P, Pdt, rate_expression, rate_leaky):
            P >> P + Pdt[lambda r1: rate_expression if r1.unbound else rate_leaky]

        def Repression(Prom, Rep, rate_binding, rate_unbinding):
            Prom.unbound + 2 * Rep >> Prom.bound[rate_binding]
            Prom.bound >> Prom.unbound + 2 * Rep[rate_unbinding]

        (
            Expression(Pcro, Cro, 5 / u.min, 0.05 / u.min),
            Expression(PcI, CI, 4.25 / u.min, 0.05 / u.min),
        )
        Repression(Pcro, CI, 1, 50**2), Repression(PcI, Cro, 1, 40**2)

        model = set_counts(
            {
                Pcro: 1,
                PcI: 1,
                Cro: 0,
                CI: 10,
                Reporter: 0,
                TtrR: 1,
                TtrS: 1,
                TTR: 0,
                Dummy: 0,
            }
        )
        S = Simulation(model)
        S.duration = 160 * u.hour

        with S.event_time(25 * u.hour):
            TTR(1)
        with S.event_time(60 * u.hour):
            TTR(0)
        with S.event_time(90 * u.hour):
            Dummy(1)
        with S.event_time(130 * u.hour):
            Dummy(0)

        S.plot_data = False
        S.save_data = False
        S.unit_x, S.unit_y = u.hour, 1 / u.ml
        S.run()
        assert S.results is not None


class TestKappaComparison:
    def test_compile(self):
        Link1, Link2, Phos_2 = BaseSpecies()
        Link1.nl_1, Link1.l_1, Link2.nl_2, Link2.l_2
        Phos_2.not_phos, Phos_2.phos_1, Phos_2.phos_2

        A = Link1 * Link2
        B = New(Link1)
        C = Phos_2 * Link2

        Rev[A.nl_1 + B.nl_1 >> A.l_1 + B.l_1][1e-4, lambda r: 0.1 * r]
        A.nl_2.l_1 + C.nl_2.not_phos >> A.l_2.l_1 + C.l_2.not_phos[1e-4]
        C.l_2.not_phos + A.l_2.l_1 >> C.nl_2.phos_1 + A.nl_2.l_1[lambda r: 1 * r]
        A.l_2.nl_1 + C.nl_2.phos_1 >> A.l_2.nl_1 + C.l_2.phos_1[1e-4]
        A.l_2.nl_1 + C.l_2.phos_1 >> A.nl_2.nl_1 + C.nl_2.phos_2[lambda r: r]

        A(1000), B(1000), C(10000)
        S = Simulation(A | B | C)
        result = S.compile()
        assert result is not None


class TestMutualAnnihilation:
    def test_run(self):
        A, B = BaseSpecies()
        A >> 2 * A[1.05 / u.h]
        B >> 2 * B[1 / u.h]
        A(1 / u.ml), B(1 / u.ml)

        S1 = Simulation(A | B)
        S1.duration = 3 * u.h
        S1.volume = 1 * u.ml

        A + B >> Zero[0.1 / u.h]

        S2 = Simulation(A | B)
        S2.duration = (A <= 0) | (B <= 0)
        S2.method = "stochastic"
        S2.volume = 1 * u.ml

        S = S1 + S2
        S.unit_x, S.unit_y = u.h, 1 / u.ml
        S.save_data = False
        S.plot_data = False
        S.repetitions = 2
        S.run()
        assert S.results is not None


class TestPySBComparison:
    def test_compile(self):
        L, R = BaseSpecies()
        L_0, R_0, kf, kr = ModelParameters(100, 200, 1e-3, 1e-3)
        L.sl_0 + R.sr_0 >> L.sl_1 + R.sr_1[kf, lambda r: kr * r]
        L(L_0), R(R_0)
        S = Simulation(L | R)
        result = S.compile()
        assert result is not None


class TestSynchronizedBacterialLysis:
    def test_run(self):
        (
            n_0,
            lysis_0,
            k,
            b,
            mu_g,
            mu,
            c_l,
            alpha_0,
            alpha_h,
            AHL_0,
            gamma_l,
            c_i,
            gamma_i,
            gamma_c,
        ) = (
            10,
            2,
            10 / u.h,
            25 / u.h,
            0.2 / u.h,
            12 / u.h,
            0.5,
            0.5 / u.h,
            35 / u.h,
            5,
            2 / u.h,
            1,
            2 / u.h,
            12 / u.h,
        )

        Cell, Lysis, AHL, LuxI = BaseSpecies()

        Cell >> 2 * Cell[lambda cell: mu_g * cell * (n_0 - cell)]
        (
            Lysis + Cell
            >> Zero[lambda lysis, cell: k * cell / (1 + (lysis_0 / lysis) ** 2)]
        )

        Cell + LuxI >> AHL + Cell + LuxI[b]
        AHL + Cell >> Cell[lambda ahl, cell: mu * ahl / (1 + cell / n_0)]

        (
            AHL
            >> AHL
            + Lysis[
                lambda ahl: (
                    c_l
                    * (
                        alpha_0
                        + alpha_h * (ahl / AHL_0) ** 4 / (1 + (ahl / AHL_0) ** 4)
                    )
                )
            ]
        )
        Lysis >> Zero[gamma_l + mu_g]

        (
            AHL
            >> AHL
            + LuxI[
                lambda ahl: (
                    c_i
                    * (
                        alpha_0
                        + alpha_h * (ahl / AHL_0) ** 4 / (1 + (ahl / AHL_0) ** 4)
                    )
                )
            ]
        )
        LuxI >> Zero[gamma_i + mu_g + gamma_c]

        Cell(5 / u.l), Lysis(0 / u.l), AHL(1e-5 / u.l), LuxI(1e-5 / u.l)

        S = Simulation(Cell | Lysis | AHL | LuxI)
        S.save_data = False
        S.plot_data = False
        S.duration = 10 * u.hour
        S.unit_x, S.unit_y = u.hour, 1 / u.ml
        S.run()
        assert S.results is not None


class TestPhageTransmissionSystem:
    def test_compile(self):
        volume = 1 * u.microliter
        c_donor = 10000 * volume / u.microliter
        c_rec = 1000 * volume / u.microliter
        c_r1 = 3 * 9e6 * volume / u.microliter
        c_r2 = 3e8 * volume / u.microliter
        Mortal, Resource, Age, Dead = BaseSpecies()
        Phage = New(Mortal)
        Donor = New(Mortal * Age)
        Receiver = New(Age)
        R1, R2 = New(Resource)

        Age.young >> Age.old[1 / 4 * (1 / u.min)]

        rm = lambda r: 1 / c_r1 if r.is_a(R1) else 1 / c_r2
        cm = lambda r: 1 / c_donor if r.is_a(Donor) else 1 / c_rec
        grw_r = lambda r1, r2: (
            1 / 20 * cm(r1) * rm(r2) * (u.l / u.s)
            if r2.is_a(R1)
            else 0.08 / 20 * cm(r1) * rm(r2) * (u.l / u.s)
        )
        inf_r = lambda r1, r2: (
            10 * 3e-11 * cm(r1) * rm(r2) * (u.l / u.s)
            if r1.old
            else 0.004 * 10 * 3e-11 * cm(r1) * rm(r2) * (u.l / u.s)
        )
        Receiver.not_infected + Phage >> Receiver.early_infection[inf_r]
        Receiver.early_infection >> Receiver.late_infection[1 / (3 * u.min)]

        Receiver.old + Resource >> Receiver.young + Receiver.not_infected.young[grw_r]
        Receiver >> Dead[1e-4 / u.s]

        phage_rate = lambda r1, r2: (
            cm(r1) * rm(r2) * 850 * (u.l / u.s) if r2.is_a(R1) else 0
        )
        Donor.old + Resource >> 2 * Donor.young[lambda r1, r2: grw_r(r1, r2) / 2]
        Donor + Resource >> Donor + Phage + Resource[phage_rate]

        (
            Mortal
            >> Zero[
                lambda r: (
                    1e-4 * cm(r) * 1 / u.s
                    if r.is_a(Donor)
                    else 0.074 * cm(r) / 24 * (1 / u.min)
                )
            ]
        )
        (
            Dead + Phage
            >> Dead[lambda r1, r2: 0.074 / 24 * cm(r1) * rm(r2) * (u.l / u.min)]
        )

        model = set_counts(
            {
                Receiver.not_infected: c_rec,
                Donor: c_donor,
                Phage: 0,
                R1: c_r1,
                R2: c_r2,
                Dead: 0,
            }
        )
        S = Simulation(model)
        S.duration = 2000 * u.s
        S.volume = volume
        result = S.compile()
        assert result is not None


class TestSBMLExport:
    def test_generate_sbml(self):
        A, B = BaseSpecies()
        A(50)
        A >> B[0.3]
        S = Simulation(A | B)
        S.duration = 10
        sbml = S.generate_sbml()
        assert isinstance(sbml, list)
        assert len(sbml) > 0
        assert "sbml" in sbml[0].lower()
