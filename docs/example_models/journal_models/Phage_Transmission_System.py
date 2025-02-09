from mobspy import *
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    volume = 1 * u.microliter
    c_donor, c_rec = 10000 * volume / u.microliter, 1000 * volume / u.microliter
    c_r1, c_r2 = 3 * 9e6 * volume / u.microliter, 3e8 * volume / u.microliter
    Mortal, Resource, Age, Dead = BaseSpecies()
    Phage = New(Mortal)
    Donor = New(Mortal * Age)
    Receiver = New(Age)
    R1, R2 = New(Resource)

    Age.young >> Age.old[1 / 4 * (1 / u.min)]

    # Receivers
    rm = lambda r: 1 / c_r1 if r.is_a(R1) else 1 / c_r2
    cm = lambda r: 1 / c_donor if r.is_a(Donor) else 1 / c_rec
    grw_r = lambda r1, r2: 1 / 20 * cm(r1) * rm(r2) * (u.l / u.s) \
        if r2.is_a(R1) else 0.08 / 20 * cm(r1) * rm(r2) * (u.l / u.s)
    inf_r = lambda r1, r2: 10 * 3e-11 * cm(r1) * rm(r2) * (u.l / u.s) \
        if r1.old else 0.004 * 10 * 3e-11 * cm(r1) * rm(r2) * (u.l / u.s)
    Receiver.not_infected + Phage >> Receiver.early_infection[inf_r]
    Receiver.early_infection >> Receiver.late_infection[1 / (3 * u.min)]

    # Receiver replication and death
    Receiver.old + Resource >> Receiver.young + Receiver.not_infected.young[grw_r]
    Receiver >> Dead[1e-4 / u.s]

    # Donors Reactions
    phage_rate = lambda r1, r2: cm(r1) * rm(r2) * 850 * (u.l / u.s) if r2.is_a(R1) else 0
    Donor.old + Resource >> 2 * Donor.young[lambda r1, r2: grw_r(r1, r2) / 2]
    Donor + Resource >> Donor + Phage + Resource[phage_rate]

    # Death Reactions
    Mortal >> Zero[lambda r: 1e-4 * cm(r) * 1 / u.s if r.is_a(Donor) else 0.074 * cm(r) / 24 * (1 / u.min)]
    Dead + Phage >> Dead[lambda r1, r2: 0.074 / 24 * cm(r1) * rm(r2) * (u.l / u.min)]

    model = set_counts({Receiver.not_infected: c_rec, Donor: c_donor, Phage: 0, R1: c_r1,
                        R2: c_r2, Dead: 0})
    S1 = Simulation(model)
    S1.duration = 2000 * u.s

    S1.volume = volume
    S1.plot_data = False

    Antibiotics = BaseSpecies()
    Antibiotics + Receiver.not_infected >> Zero[0.015 * 1 / 1e9 * (u.l / u.s)]

    # Change to concentration - Is currently in counts
    Antibiotics(1e7 / u.microliter)

    S2 = Simulation(model | Antibiotics)
    S2.duration = 3 * u.hours
    S2.volume = volume
    Sim = S1 + S2
    Sim.unit_x = u.h
    Sim.unit_y = 1/u.ml
    Sim.run()
    Sim.plot_raw('plot_config_donor.json')

