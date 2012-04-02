"""The WaveBlocks Project

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

__version__ = 0.1

from Grid import Grid
from DenseGrid import DenseGrid
from TensorProductGrid import TensorProductGrid
from GridFactory import GridFactory

from WaveFunction import WaveFunction

from MatrixPotential import MatrixPotential
from MatrixPotential1S import MatrixPotential1S
from MatrixPotential2S import MatrixPotential2S
from MatrixPotentialMS import MatrixPotentialMS

from PotentialFactory import PotentialFactory



from KineticOperator import KineticOperator


from FourierPropagator import FourierPropagator



from IOManager import IOManager


from BasisShape import BasisShape
from HyperCubicShape import HyperCubicShape

#from Utils import Utils


# TODO: Recheck these files:
from BasisTransformation import BasisTransformation
from BasisTransformationWF import BasisTransformationWF

from ParameterLoader import ParameterLoader
from ParameterProvider import ParameterProvider

from TimeManager import TimeManager

from SimulationLoop import SimulationLoop
from SimulationLoopFourier import SimulationLoopFourier




# Enable dynamic plugin loading for IOManager
import sys
import os

plugin_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(plugin_dir)
