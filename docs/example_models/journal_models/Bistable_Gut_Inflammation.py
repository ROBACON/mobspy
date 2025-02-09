from mobspy import *
import os

# Conversion taking to much time - why?

Dummy, Dilutable, Promoter, TTR, Phosphorable = BaseSpecies()
TF, Reporter = New(Dilutable)
TtrR, TtrS = New(Phosphorable)
Cro, CI = New(TF)
Pcro, PcI = New(Promoter)
Promoter.bound, Promoter.unbound, Phosphorable.dephospho, Phosphorable.phospho

# TtrS triggers the state transition, TTR indicates the presence of gut inflammation
TtrS.dephospho + TTR >> TtrS.phospho + TTR [1/u.min]
TtrS.phospho + TtrR.dephospho >> TtrR.phospho + TtrS.dephospho [1/u.min]
TtrR.phospho >> TtrR.dephospho [5*1e-3*1/u.second]
TtrR.phospho >> Cro + TtrR.phospho [50*0.02*1/u.min]

# Dummy represents the introduction of a constant flow of CI in the system
Dummy >> CI + Dummy [30*0.02*1/u.min]

# Dilution and degradation
Dilutable >> Zero [0.02/u.min]
TF >> Zero [lambda r1: 0 if r1.is_a(CI) else 1.6e-2*1/u.min]

def Expression(P, Pdt, rate_expression, rate_leaky):
    P >> P + Pdt [lambda r1: rate_expression if r1.unbound else rate_leaky]

def Repression (Prom, Rep, rate_binding, rate_unbinding):
    Prom.unbound + 2*Rep >> Prom.bound [rate_binding]
    Prom.bound >> Prom.unbound + 2*Rep [rate_unbinding]

# Pcro produces Cro, and PcI produces CI
Expression(Pcro, Cro, 5/u.min, 0.05/u.min), Expression(PcI, CI, 4.25/u.min, 0.05/u.min)
# Cro represses Pcro, and cI represses PcI
Repression(Pcro, CI, 1, 50**2), Repression(PcI, Cro, 1, 40**2)

model = set_counts({Pcro: 1, PcI: 1, Cro: 0, CI: 10,
                    Reporter: 0, TtrR: 1, TtrS: 1, TTR: 0, Dummy: 0})
S = Simulation(model)
S.duration = 160*u.hour

# Event implementation
with S.event_time(25*u.hour):
    TTR(1)
with S.event_time(60*u.hour):
    TTR(0)
with S.event_time(90*u.hour):
    Dummy(1)
with S.event_time(130*u.hour):
    Dummy(0)

S.plot_data = False
S.unit_x, S.unit_y = u.hour, 1/u.ml
S.step_size, S.a_tol = 0.01*u.hour, 1e-16
S.plot_config.title, S.plot_config.title_fontsize = 'Toggle Switch', 16
S.plot_config.figsize = (6.5, 4)
S.plot_config.vertical_lines = [25, 60, 90, 130]
S.plot_config.ylabel = r"Conc. (mL$^{-1}$)"
S.plot_config.save_to = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + \
                        '/images/Toggle_Switch/Toggle_Switch.pdf'
S.run()
S.plot(CI, Cro)
