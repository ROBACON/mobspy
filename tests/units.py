# from openpyxl import Workbook
import matplotlib.pyplot as plt

from mobspy import BaseSpecies, Simulation, u

# model parameters
duration = 15 * u.min
Ks = 10 * 1 / u.mL
mu_max = 1000 * 1 / u.min

# Species
X, S = BaseSpecies()

# Initial conditions
X(10 * 1 / u.mL)
S(1e8 * 1 / u.mL)
X + S >> 2 * X[lambda x, s: mu_max * (2 * s / (Ks + 2 * s)) * x]

# Simulation
model = Simulation(S | X)
model.volume = 1 * u.mL
model.method = "deterministic"
model.duration = duration
model.unit_y = 1 / u.mL
model.plot_data = False
model.run(step_size=0.01 * u.min)

# plt.subplot(2, 1, 2)
plt.plot(model.results["Time"][0], model.results["S"][0], color="blue", label="S conc")
plt.plot(model.results["Time"][0], model.results["X"][0], color="black", label="X conc")
plt.gca().set_yscale("log")
plt.ylabel("Substrate (S)")
plt.grid()
plt.legend()
plt.show()


# Question: Why X does not asymtpotize at a value when S gets to 0?

# Data generation
# time_value = model.results["Time"][0]
# S_concentration = model.results["S"][0]
# X_concentration = model.results["X"][0]

# print(time_value, S_concentration, X_concentration)

# #Save data
# wb = Workbook()
# ws = wb.active
# ws.append(["Time", "S Concentration", "X Concentration"])
# ws.append([time_value, S_concentration, X_concentration])
# wb.save("model_results.xlsx")
