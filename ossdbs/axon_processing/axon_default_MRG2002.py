# Copyright 2024 Konstantin Butenko, Julius Zimmermann
# Copyright 2017 Christian Schmidt
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np


# Disable Ruff warning regarding complexity
# ruff: noqa: C901
def get_axon_parameters_template(diameter: float):
    """Get parameters of axon model for a given diameter.

    Notes
    -----
    Parameters for 2 um and 3 um fiber diameter are taken from:
    Sotiropoulos, S N  and Steinmetz, P N. Assessing the direct effects of deep
    brain stimulation using embedded axon models. J Neural Eng, vol. 4,
    pp. 107-19, Jun 2007.

    Returns
    -------
    dict
        ranvier_nodes: int
            Number of node of Ranvier compartments
        para1_nodes: int
            Number of first paranodal compartments
        para2_nodes: int
            Number of second paranodal compartments
        inter_nodes: int
            Number of internodal compartments
        ranvier_length: float
            Length of the Node of Ranvier
        para1_length: float
            Length of the first paranodal compartments
        para2_length: float
            Length of the second paranodal compartments
        deltax: float
            Length from a node of Ranvier to the next
        axon_diameter: float
            Diameter of the axon
        node_diameter: float
            Diameter of the nodes of Ranvier
        para1_diameter: float
            Diameter of the first paranodal compartments
        para2_diameter: float
            Diameter of the second paranodal compartments
        condg: float
            Axial conductivity
        lamellas: int
            Number of lamellas in the myelin sheath

    """
    # taken from Sotiropoulos and Steinmetz (2007)
    # condg for fiber diameter was determined
    # using the polyfit from the data with 5.7 um and larger.
    # Resulting fitting function: 0.0001x^2+0.0139x+0.5235
    if np.isclose(diameter, 2.0):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 60,
            "deltax": 200.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 10.0,
            "axon_diameter": 1.6,
            "node_diameter": 1.4,
            "para1_diameter": 1.4,
            "para2_diameter": 1.6,
            "condg": 0.552,
            "lamellas": 30,
        }

    # taken from Sotiropoulos and Steinmetz (2007)
    # condg for fiber diameter was determined
    # using the polyfit from the data with 5.7 um and larger.
    # Resulting fitting function: 0.0001x^2+0.0139x+0.5235
    elif np.isclose(diameter, 3.0):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 60,
            "deltax": 278.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 17.0,
            "axon_diameter": 2.1,
            "node_diameter": 1.52,
            "para1_diameter": 1.52,
            "para2_diameter": 2.1,
            "condg": 0.567,
            "lamellas": 43,
        }

    # (fiberD==5.7)
    # {g=0.605 axonD=3.4 nodeD=1.9 paraD1=1.9 paraD2=3.4
    # deltax=500 paralength2=35 nl=80}
    elif np.isclose(diameter, 5.7):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 500.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 35.0,
            "axon_diameter": 3.4,
            "node_diameter": 1.9,
            "para1_diameter": 1.9,
            "para2_diameter": 3.4,
            "condg": 0.605,
            "lamellas": 80,
        }
    # (fiberD==7.3)
    # {g=0.630 axonD=4.6 nodeD=2.4 paraD1=2.4 paraD2=4.6
    # deltax=750 paralength2=38 nl=100}
    elif np.isclose(diameter, 7.3):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 750.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 38.0,
            "axon_diameter": 4.6,
            "node_diameter": 2.4,
            "para1_diameter": 2.4,
            "para2_diameter": 4.6,
            "condg": 0.630,
            "lamellas": 100,
        }
    # (fiberD==8.7)
    # {g=0.661 axonD=5.8 nodeD=2.8 paraD1=2.8 paraD2=5.8
    # deltax=1000 paralength2=40 nl=110}
    elif np.isclose(diameter, 8.7):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1000.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 40.0,
            "axon_diameter": 5.8,
            "node_diameter": 2.8,
            "para1_diameter": 2.8,
            "para2_diameter": 5.8,
            "condg": 0.661,
            "lamellas": 110,
        }
    # (fiberD==10.0)
    # {g=0.690 axonD=6.9 nodeD=3.3 paraD1=3.3 paraD2=6.9
    # deltax=1150 paralength2=46 nl=120}
    elif np.isclose(diameter, 10.0):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1150.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 46.0,
            "axon_diameter": 6.9,
            "node_diameter": 3.3,
            "para1_diameter": 3.3,
            "para2_diameter": 6.9,
            "condg": 0.690,
            "lamellas": 120,
        }
    # (fiberD==11.5)
    # {g=0.700 axonD=8.1 nodeD=3.7 paraD1=3.7 paraD2=8.1
    # deltax=1250 paralength2=50 nl=130}
    elif np.isclose(diameter, 11.5):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1250.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 50.0,
            "axon_diameter": 8.1,
            "node_diameter": 3.7,
            "para1_diameter": 3.7,
            "para2_diameter": 8.1,
            "condg": 0.700,
            "lamellas": 130,
        }
    # (fiberD==12.8)
    # {g=0.719 axonD=9.2 nodeD=4.2 paraD1=4.2 paraD2=9.2
    # deltax=1350 paralength2=54 nl=135}
    elif np.isclose(diameter, 12.8):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1350.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 54.0,
            "axon_diameter": 9.2,
            "node_diameter": 4.2,
            "para1_diameter": 4.2,
            "para2_diameter": 9.2,
            "condg": 0.719,
            "lamellas": 135,
        }
    # (fiberD==14.0)
    # {g=0.739 axonD=10.4 nodeD=4.7 paraD1=4.7 paraD2=10.4
    # deltax=1400 paralength2=56 nl=140}
    elif np.isclose(diameter, 14.0):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1400.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 56.0,
            "axon_diameter": 10.4,
            "node_diameter": 4.7,
            "para1_diameter": 4.7,
            "para2_diameter": 10.4,
            "condg": 0.739,
            "lamellas": 140,
        }
    # (fiberD==15.0)
    # {g=0.767 axonD=11.5 nodeD=5.0 paraD1=5.0 paraD2=11.5
    # deltax=1450 paralength2=58 nl=145}
    elif np.isclose(diameter, 15.0):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1450.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 58.0,
            "axon_diameter": 11.5,
            "node_diameter": 5.0,
            "para1_diameter": 5.0,
            "para2_diameter": 11.5,
            "condg": 0.767,
            "lamellas": 145,
        }
    # (fiberD==16.0)
    # {g=0.791 axonD=12.7 nodeD=5.5 paraD1=5.5 paraD2=12.7
    # deltax=1500 paralength2=60 nl=150}
    elif np.isclose(diameter, 16.0):
        return {
            "ranvier_nodes": 21,
            "para1_nodes": 40,
            "para2_nodes": 40,
            "inter_nodes": 120,
            "deltax": 1500.0,
            "ranvier_length": 1.0,
            "para1_length": 3.0,
            "para2_length": 60.0,
            "axon_diameter": 12.7,
            "node_diameter": 5.5,
            "para1_diameter": 5.5,
            "para2_diameter": 12.7,
            "condg": 0.791,
            "lamellas": 150,
        }
    else:
        raise NotImplementedError(
            f"Axon parameters for diameter {diameter} are not yet implemented."
        )
