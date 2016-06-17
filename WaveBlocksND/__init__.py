"""The WaveBlocks Project

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 R. Bourquin
@license: Modified BSD License
"""

from __future__ import absolute_import

__version__ = 0.5

# Math
from WaveBlocksND.ComplexMath import ContinuousSqrt

# Grids
from WaveBlocksND.AbstractGrid import AbstractGrid
from WaveBlocksND.DenseGrid import DenseGrid
from WaveBlocksND.TensorProductGrid import TensorProductGrid
from WaveBlocksND.GridWrapper import GridWrapper

from WaveBlocksND.PhaseSpaceLattice import PhaseSpaceLattice

# Wavefunctions
from WaveBlocksND.WaveFunction import WaveFunction

# Potentials
from WaveBlocksND.MatrixPotential import MatrixPotential
from WaveBlocksND.MatrixPotential1S import MatrixPotential1S
from WaveBlocksND.MatrixPotential2S import MatrixPotential2S
from WaveBlocksND.MatrixPotentialMS import MatrixPotentialMS

from WaveBlocksND.KineticOperator import KineticOperator

# Time Propagators
from WaveBlocksND.Propagator import Propagator
from WaveBlocksND.FourierPropagator import FourierPropagator
from WaveBlocksND.HagedornPropagator import HagedornPropagator
from WaveBlocksND.MagnusPropagator import MagnusPropagator
from WaveBlocksND.SemiclassicalPropagator import SemiclassicalPropagator
from WaveBlocksND.McL42scPropagator import McL42scPropagator
from WaveBlocksND.McL84scPropagator import McL84scPropagator
from WaveBlocksND.Pre764scPropagator import Pre764scPropagator
from WaveBlocksND.HagedornPropagatorInhomogeneous import HagedornPropagatorInhomogeneous
from WaveBlocksND.SplittingParameters import SplittingParameters
from WaveBlocksND.PerturbedSplittingParameters import PerturbedSplittingParameters
from WaveBlocksND.ProcessingSplittingParameters import ProcessingSplittingParameters

from WaveBlocksND.IOManager import IOManager

# Basis shapes
from WaveBlocksND.BasisShape import BasisShape
from WaveBlocksND.HyperCubicShape import HyperCubicShape
from WaveBlocksND.SimplexShape import SimplexShape
from WaveBlocksND.HyperbolicCutShape import HyperbolicCutShape
from WaveBlocksND.LimitedHyperbolicCutShape import LimitedHyperbolicCutShape

# Wavepackets
from WaveBlocksND.Wavepacket import Wavepacket
from WaveBlocksND.HagedornBasisEvaluationCommon import HagedornBasisEvaluationCommon
from WaveBlocksND.HagedornBasisEvaluationPhi import HagedornBasisEvaluationPhi
from WaveBlocksND.HagedornBasisEvaluationPsi import HagedornBasisEvaluationPsi
from WaveBlocksND.HagedornWavepacketBase import HagedornWavepacketBase
from WaveBlocksND.HagedornWavepacket import HagedornWavepacket
from WaveBlocksND.HagedornWavepacketNew import HagedornWavepacketNew
from WaveBlocksND.HagedornWavepacketInhomogeneous import HagedornWavepacketInhomogeneous

from WaveBlocksND.Gradient import Gradient
from WaveBlocksND.WavepacketGradient import WavepacketGradient
from WaveBlocksND.GradientHAWP import GradientHAWP
from WaveBlocksND.GradientHAWPnew import GradientHAWPnew
from WaveBlocksND.PositionHAWP import PositionHAWP

# Linear combinations
from WaveBlocksND.LinearCombinationOfWavepackets import LinearCombinationOfWavepackets
from WaveBlocksND.LinearCombinationOfWPs import LinearCombinationOfWPs
from WaveBlocksND.LinearCombinationOfHAWPs import LinearCombinationOfHAWPs

from WaveBlocksND.GradientLCWP import GradientLCWP
from WaveBlocksND.GradientLinearCombinationHAWP import GradientLinearCombinationHAWP

# Quadrature rules
from WaveBlocksND.QuadratureRule import QuadratureRule
from WaveBlocksND.TrapezoidalQR import TrapezoidalQR
from WaveBlocksND.GaussHermiteOriginalQR import GaussHermiteOriginalQR
from WaveBlocksND.GaussHermiteQR import GaussHermiteQR
from WaveBlocksND.GaussLaguerreQR import GaussLaguerreQR
from WaveBlocksND.GenzKeisterOriginalQR import GenzKeisterOriginalQR
from WaveBlocksND.GenzKeisterQR import GenzKeisterQR
from WaveBlocksND.TensorProductQR import TensorProductQR
from WaveBlocksND.SmolyakQR import SmolyakQR

# Inner products
from WaveBlocksND.InnerProduct import InnerProduct
from WaveBlocksND.InnerProduct import InnerProductException
from WaveBlocksND.HomogeneousInnerProduct import HomogeneousInnerProduct
from WaveBlocksND.InhomogeneousInnerProduct import InhomogeneousInnerProduct

from WaveBlocksND.HomogeneousInnerProductLCWP import HomogeneousInnerProductLCWP
from WaveBlocksND.InhomogeneousInnerProductLCWP import InhomogeneousInnerProductLCWP

from WaveBlocksND.Quadrature import Quadrature
from WaveBlocksND.DirectQuadrature import DirectQuadrature
from WaveBlocksND.DirectHomogeneousQuadrature import DirectHomogeneousQuadrature
from WaveBlocksND.DirectInhomogeneousQuadrature import DirectInhomogeneousQuadrature

from WaveBlocksND.SymbolicIntegral import SymbolicIntegral
from WaveBlocksND.SymbolicIntegral0 import SymbolicIntegral0
from WaveBlocksND.GaussianIntegral import GaussianIntegral
from WaveBlocksND.NSDInhomogeneous import NSDInhomogeneous

from WaveBlocksND.SparsityOracle import SparsityOracle
from WaveBlocksND.SparsityOracleTrue import SparsityOracleTrue
from WaveBlocksND.SparsityOraclePSHAWP import SparsityOraclePSHAWP
from WaveBlocksND.SparsityOracleGIHAWP import SparsityOracleGIHAWP

# Basis transformations
from WaveBlocksND.BasisTransformation import BasisTransformation
from WaveBlocksND.BasisTransformationWF import BasisTransformationWF
from WaveBlocksND.BasisTransformationHAWP import BasisTransformationHAWP

# TODO: Recheck these files:
from WaveBlocksND.BlockFactory import BlockFactory

from WaveBlocksND.ParameterLoader import ParameterLoader
from WaveBlocksND.ParameterProvider import ParameterProvider

from WaveBlocksND.TimeManager import TimeManager

from WaveBlocksND.SimulationLoop import SimulationLoop
from WaveBlocksND.SimulationLoopFourier import SimulationLoopFourier
from WaveBlocksND.SimulationLoopHagedorn import SimulationLoopHagedorn
from WaveBlocksND.SimulationLoopHagedornInhomogeneous import SimulationLoopHagedornInhomogeneous

from WaveBlocksND.Observables import Observables
from WaveBlocksND.ObservablesHAWP import ObservablesHAWP
from WaveBlocksND.ObservablesMixedHAWP import ObservablesMixedHAWP
from WaveBlocksND.ObservablesLCWP import ObservablesLCWP

# Enable dynamic plugin loading for IOManager
import sys as _sys
import os as _os

_plugin_dir = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.append(_plugin_dir)
