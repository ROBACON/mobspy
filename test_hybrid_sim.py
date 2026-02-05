# This script tests the model creation and compilation
# It also test the calculation capabilities
# It uses some simple model and assertions

import mobspy
from mobspy import (
    BaseSpecies,
    Simulation,
    New,
    u,
    ModelParameters,
    simlog,
    Assign,
    Rev,
    All,
    set_counts,
)

# from tellurium import loada as te_load_anti
import numpy as np
from copy import deepcopy
import os
import re

from testutils import compare_model, compare_model_ignore_order


def test_hybrid_sim():
    A, B = BaseSpecies(2)
    A >> 2 * A[1]

    A(1), B(50)
    S1 = Simulation(A)
    S1.save_data = False
    S1.plot_data = False
    S1.duration = 3

    A.reset_reactions()
    A + B >> mobspy.Zero[0.01]

    S2 = Simulation(A | B)
    S2.method = "stochastic"
    S2.duration = (A <= 0) | (B <= 0)
    S2.plot_data = False

    Sim = S1 + S2
    Sim.run(plot_data=False)

    assert compare_model(Sim.compile(), "test_tools/model_8.txt")
    assert Sim.fres[A][-1] == 0 or Sim.fres[B][-1] == 0


if __name__ == "__main__":
    test_hybrid_sim()
