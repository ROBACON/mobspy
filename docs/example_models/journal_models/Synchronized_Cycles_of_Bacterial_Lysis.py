from mobspy import *
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# n_0 = 10 - dimensionless
# lysis_0 = 2
# n = 2 - dimensionless
# k = 10 - 1/h
# b = 25 - volume/h
# mu_g = 0.2 - 1/h
# mu = 12 - 1/h
# c_l = 0.5 - dimensionless
# alpha_0 = 0.5 - volume/h
# alpha_h = 35 - volume/h
# AHL_0 = 5 - dimensionless
# gamma_l = 2 - 1/h
# c_i = 1 - dimensionless
# gamma_i = 2 - 1/h
# gamma_c = 12 - 1/h

n_0, lysis_0, k, b, mu_g, mu, c_l, alpha_0, alpha_h, AHL_0, gamma_l, c_i, \
    gamma_i, gamma_c = (10, 2, 10/u.h, 25/u.h, 0.2/u.h, 12/u.h, 0.5, 0.5/u.h, 35/u.h, 5,
                        2/u.h, 1, 2/u.h, 12/u.h)

# N, L, H, I - Are the respective conterparts of the meta-species below
# in the original paper's model
Cell, Lysis, AHL, LuxI  = BaseSpecies()

# Cell related reactions
Cell >> 2*Cell [lambda cell: mu_g*cell*(n_0 - cell)]
# Lysis enzyme encounters the cell membrane through the inside of the cell and kills it
Lysis + Cell >> Zero [lambda lysis, cell: k*cell/(1 + (lysis_0/lysis)**2)]

# AHL related reactions
Cell + LuxI >> AHL + Cell + LuxI [b]
AHL + Cell >> Cell [lambda ahl, cell: mu*ahl/(1 + cell/n_0)]

# Lysis related reactions
AHL >> AHL + Lysis [lambda ahl:  c_l*(alpha_0 + alpha_h*(ahl/AHL_0)**4/(1 + (ahl/AHL_0)**4))]
Lysis >> Zero [gamma_l + mu_g]

# LuxI related reactions
AHL >> AHL + LuxI [lambda ahl:  c_i*(alpha_0 + alpha_h*(ahl/AHL_0)**4/(1 + (ahl/AHL_0)**4))]
LuxI >> Zero [gamma_i + mu_g + gamma_c]

Cell(5 / u.l), Lysis(0 / u.l), AHL(1e-5 / u.l), LuxI(1e-5 / u.l)

MySim = Simulation(Cell | Lysis | AHL | LuxI)
MySim.save_data = False
MySim.duration = 10*u.hour
MySim.unit_x, MySim.unit_y = u.hour, 1 / u.ml
MySim.plot_config.title, MySim.plot_config.title_fontsize = 'Synchronized Bacterial Lysis', 18
MySim.plot_config.xlabel_fontsize, MySim.plot_config.ylabel_fontsize = 14, 14
MySim.plot_config.ylabel = r"Conc. (mL$^{-1}$)"
MySim.plot_config.save_to = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + \
                            '/images/Lysis_Clock/Lysis_Clock.pdf'
MySim.run()


