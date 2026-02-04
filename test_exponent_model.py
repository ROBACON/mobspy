"""Model for the exponent tests in MobsPy."""
from mobspy import BaseSpecies, Simulation, u, New
from dataclasses import dataclass
from pint import Quantity

@dataclass
class ModelParameters:
    default_rate: Quantity
    k_tet : Quantity
    coeff : float

@dataclass
class InitialQuantities:
    qty_sender: Quantity
    qty_tet: Quantity

class Model:
    params: ModelParameters

    def __init__(self) -> None:
        default_rate = 3e-11 / u.min
        k_tet = 1208585 * u.counts / u.mL
        coeff = 1.0

        self.params = ModelParameters(
            default_rate=default_rate,
            k_tet=k_tet,
            coeff=coeff,
        )

    def get_simulation_exp(
        self,
        initial: InitialQuantities,
    ) -> Simulation:
        params = ModelParameters(**self.params.__dict__)

        def antibio_uptake_senders(b):
            #breakpoint()
            mu_monod = (params.default_rate)*(params.k_tet**params.coeff/((params.k_tet**params.coeff)+initial.qty_sender**params.coeff))
            rate_monod = mu_monod * b
            return rate_monod

        # Base species definition & reactions under
        Antibiotic = BaseSpecies()
        Tet = New(Antibiotic)
        Tet(initial.qty_tet)

        tet = BaseSpecies()
        tet.tetF, tet.tetT
        Bacterium = tet
        Sender = New(Bacterium)
        Sender(initial.qty_sender)

        Sender.tetF + Tet >> Sender.tetT [antibio_uptake_senders] 

        return Simulation( Sender | Tet )
    


    def get_simulation_default(
        self,
        initial: InitialQuantities,
    ) -> Simulation:
        params = ModelParameters(**self.params.__dict__)

        
        #Without the exponent, to try 
        def antibio_uptake_senders(b):
            #breakpoint()
            mu_monod = (params.default_rate)*(params.k_tet/((params.k_tet)+initial.qty_sender))
            rate_monod = mu_monod * b
            return rate_monod


        # Base species definition & reactions under
        Antibiotic = BaseSpecies()
        Tet = New(Antibiotic)
        Tet(initial.qty_tet)

        tet = BaseSpecies()
        tet.tetF, tet.tetT
        Bacterium = tet
        Sender = New(Bacterium)
        Sender(initial.qty_sender)

        Sender.tetF + Tet >> Sender.tetT [antibio_uptake_senders] 

        return Simulation( Sender | Tet )

