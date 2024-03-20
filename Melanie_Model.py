import numpy as np
from scipy.optimize import minimize, fsolve, curve_fit, leastsq, differential_evolution
from scipy.integrate import odeint
from scipy.signal import find_peaks, argrelextrema, savgol_filter
import pandas as pd
import matplotlib.pyplot as plt
import tqdm
from mobspy import *


# Use Python env
# New version of mobspy
# python3 -m venv bcrn  to create the environnement
# source bcrn/bin/activate to activate the environnement

# Old version of mobspy
# python3 -m venv modinfalgue  to create the environnement
# source modinfalgue/bin/activate to activate the environnement

### STEP 0 : IMPORT DATA

# Load .csv file and remove 'nan' values
def load_and_remove_nan(file, x_column, y_column, rows_to_skip):
    df = pd.read_csv(file, delimiter=';', skiprows=rows_to_skip)
    x = df[x_column].tolist()
    y = df[y_column].tolist()

    # Remove 'nan' values and corresponding entries
    cleaned_indices = [i for i, val in enumerate(y) if not np.isnan(val)]
    x = [x[i] for i in cleaned_indices]
    y = [y[i] for i in cleaned_indices]

    return x, y


# file = 'growth_excell_format.xlsx'
# df_modif = pd.read_csv(file, sep=';', index_col=0)

df_modif = pd.read_excel('growth_excell_format.xlsx', index_col=0)

def BCRN():
    kappa_1, kappa_2, kappa_3, kappa_4, kappa_5 = ModelParameters(0.5 / u.second, 1 / u.second, 1 / u.second, 1 / u.second, 1 / u.second)

    print('kappa: ', kappa_1, kappa_2, kappa_3, kappa_4, kappa_5)

    Microalga, Nutrient, Light, Complex = BaseSpecies(4)
    Complex.activated, Complex.inactivated

    # Nutrient reaction
    Microalga + Nutrient + Light >> Complex.inactivated + Light[kappa_1]
    Complex.inactivated + Light >> Microalga + Nutrient + Light[kappa_2]

    # Light reaction
    Complex.inactivated + Light >> Complex.activated[kappa_3]
    Complex.activated >> Complex.inactivated + Light[kappa_4]

    # Mitosis
    Complex.activated >> Microalga + Microalga[kappa_5]

    # Simulation
    Microalga(0.125 * 200 * 30.1e6), Nutrient(8.7e12), Complex(0), Light(0)
    S = Simulation(Microalga | Nutrient | Complex | Light)
    S.method = 'deterministic'
    S.duration = 1200 * 3600 * u.second  # à changer selon les données
    S.unit_x = u.second
    # S.plot_data = False
    S.level = 1  # Pour supprimer les messages pendant la simulation
    S.run()

    exit()

    basiCO_parameter_estimation(S, experimental_data=df_modif,
                                parameters_to_estimate=[kappa_1, kappa_2, kappa_3, kappa_4, kappa_5])


BCRN()
exit()

time_exp = df_modif['Time'].tolist()
cell_density_exp = [float(x) for x in df_modif['Microalga'].tolist()]

# plt.plot(time_simulated, cell_density_simulation, color = 'red')
# plt.plot(new_time_simulated, new_cell_density_simulation, color = 'coral')
plt.scatter(time_exp, cell_density_exp, color='blue', alpha=0.7, label='Experimental data')
plt.xlabel('Time (seconds)')
plt.ylabel('Number of cells')
plt.title('Growth of ACCB1808')
plt.legend()
plt.show()

# Stockage des résultats
results = []

# Pour chaque jeu de valeurs d'arguments
for i in range(10):  # Par exemple, pour 10 jeux de valeurs aléatoires
    # Génération de nouvelles valeurs aléatoires pour chaque courbe
    param_values = np.random.uniform(1e-12, 100, size=5)  # Génère des valeurs aléatoires entre 2 et 5

    # Appel de la fonction BCRN avec les valeurs d'arguments actuelles
    time_simulated, cell_density_simulation, new_time_simulated, new_cell_density_simulation = BCRN(param_values)

    # Tracer la courbe correspondante
    plt.plot(new_time_simulated, new_cell_density_simulation, label=f'Curve {i + 1}')

    # Stocker les résultats
    results.append(
        (param_values, time_simulated, cell_density_simulation, new_time_simulated, new_cell_density_simulation))

# Afficher le graphe
plt.xlabel('Time (seconds)')
plt.ylabel('Number of cells')
plt.title('Growth of ACCB1808')
plt.legend()
plt.show()

# Afficher les valeurs des arguments adoptées pour chaque courbe dans la console
for i, (param_values, time_simulated, cell_density_simulation, new_time_simulated, new_cell_density_simulation) \
        in enumerate(results): print(f"Courbe {i + 1} : {param_values}")