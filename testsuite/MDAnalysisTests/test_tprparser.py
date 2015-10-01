# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from MDAnalysis.tests.datafiles import TPR, \
    TPR400, TPR402, TPR403, TPR404, TPR405, TPR406, TPR407, \
    TPR450, TPR451, TPR452, TPR453, TPR454, TPR455, TPR455Double, \
    TPR460, TPR461, TPR502, TPR504, TPR505, TPR510, TPR510_bonded

from numpy.testing import TestCase, dec
from test_topology import _TestTopology
import MDAnalysis.topology.TPRParser


@dec.slow
class RefTPR(object):
    """
    this test the data/adk_oplsaa.tpr which is of tpx version 58
    """
    topology = TPR
    parser = MDAnalysis.topology.TPRParser.TPRParser
    ref_n_atoms = 47681
    ref_numresidues = 11302
    ref_proteinatoms = 3341


class TestTPR(_TestTopology, RefTPR):
    """Testing TPR version 58"""


# The follow test the same system grompped by different version of gromacs
# FORMAT: TPRABC, where numbers ABC indicates the version of gromacs that
# generates the corresponding tpr file

class TPRBase(object):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    ref_n_atoms = 2263
    ref_numresidues = 230
    ref_proteinatoms = 1962


@dec.slow
class TPR400(TPRBase):
    topology = TPR400


class TestTPR400(_TestTopology, TPR400):
    """Testing TPR version 58"""


@dec.slow
class TPR402(TPRBase):
    topology = TPR402


class TestTPR402(_TestTopology, TPR402):
    """Testing TPR version 58"""


@dec.slow
class TPR403(TPRBase):
    topology = TPR403


class TestTPR403(_TestTopology, TPR403):
    """Testing TPR version 58"""


@dec.slow
class TPR404(TPRBase):
    topology = TPR404


class TestTPR404(_TestTopology, TPR404):
    """Testing TPR version 58"""


@dec.slow
class TPR405(TPRBase):
    topology = TPR405


class TestTPR405(_TestTopology, TPR405):
    """Testing TPR version 58"""


@dec.slow
class TPR406(TPRBase):
    topology = TPR406


class TestTPR406(_TestTopology, TPR406):
    """Testing TPR version 58"""


@dec.slow
class TPR407(TPRBase):
    topology = TPR407


class TestTPR407(_TestTopology, TPR407):
    """Testing TPR version 58"""


@dec.slow
class TPR450(TPRBase):
    topology = TPR450


class TestTPR450(_TestTopology, TPR450):
    """Testing TPR version 73"""


@dec.slow
class TPR451(TPRBase):
    topology = TPR451


class TestTPR451(_TestTopology, TPR451):
    """Testing TPR version 73"""


@dec.slow
class TPR452(TPRBase):
    topology = TPR452


class TestTPR452(_TestTopology, TPR452):
    """Testing TPR version 73"""


@dec.slow
class TPR453(TPRBase):
    topology = TPR453


class TestTPR453(_TestTopology, TPR453):
    """Testing TPR version 73"""


@dec.slow
class TPR454(TPRBase):
    topology = TPR454


class TestTPR454(_TestTopology, TPR454):
    """Testing TPR version 73"""


@dec.slow
class TPR455(TPRBase):
    topology = TPR455


class TestTPR455(_TestTopology, TPR455):
    """Testing TPR version 73"""


@dec.slow
class TPR455Double(object):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    ref_n_atoms = 21692
    ref_numresidues = 4352
    ref_proteinatoms = 0  # no protein, but DOPC, DPPC, CHOL, SOL
    topology = TPR455Double


class TestTPR455Double(_TestTopology, TPR455Double):
    """Testing TPR version 73, double precision"""


class TPR46xBase(object):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    ref_n_atoms = 44052
    ref_numresidues = 10712
    ref_proteinatoms = 1885


@dec.slow
class TPR460(TPR46xBase):
    topology = TPR460


class TestTPR460(_TestTopology, TPR460):
    """Testing TPR version 83"""


@dec.slow
class TPR461(TPR46xBase):
    topology = TPR461


class TestTPR461(_TestTopology, TPR461):
    """Testing TPR version 83"""

@dec.slow
class TPR502(TPRBase):
    topology = TPR502


class TestTPR502(_TestTopology, TPR502):
    """Testing TPR version 100"""


@dec.slow
class TPR504(TPRBase):
    topology = TPR504


class TestTPR504(_TestTopology, TPR504):
    """Testing TPR version 100"""

@dec.slow
class TPR505(TPRBase):
    topology = TPR505


class TestTPR505(_TestTopology, TPR505):
    """Testing TPR version 100"""

@dec.slow
class TPR510(TPRBase):
    topology = TPR510


class TestTPR510(_TestTopology, TPR510):
    """Testing TPR version 103"""

class TestAllBondedTPR510(TestCase):
    """Test that all bonded interactions are read from a TPR version 103
    
    For a list of the interactions, see table 5.5 of the Gromacs PDF reference.
    """
    topology = TPR510_bonded
    
    def setUp(self):
        """The topology should be parsed previous to each test"""
        self.universe = MDAnalysis.Universe(self.topology)

    def test_regular_bond(self):
        """Test for regular bounds (BONDS)"""
        self.assertIn((0, 1), self.universe._topology['bonds'])

    def test_G96_bonds(self):
        """Test for G96 bonds (G96BONDS)"""
        self.assertIn((1, 2), self.universe._topology['bonds'])

    def test_morse(self):
        """Test for Morse bonds (MORSE)"""
        self.assertIn((2, 3), self.universe._topology['bonds'])

    def test_cubic_bonds(self):
        """Test for cubic bonds (CUBICBONDS)"""
        self.assertIn((3, 4), self.universe._topology['bonds'])

    def test_connection(self):
        """Test for connections (CONNBONDS)"""
        self.assertIn((4, 5), self.universe._topology['bonds'])

    def test_harmonic_potentials(self):
        """Test for harmonic potentials (HARMONIC)"""
        self.assertIn((5, 6), self.universe._topology['bonds'])

    def test_fene_bonds(self):
        """Test for FENE bonds (FENEBONDS)"""
        self.assertIn((6, 7), self.universe._topology['bonds'])

    def test_restraint_potentials(self):
        """Test for restraint potentials (RESTRAINTPOT)"""
        self.assertIn((7, 8), self.universe._topology['bonds'])
 
    def test_tabulated_bond(self):
        """Test for tabulated bonds (TABBONDS)"""
        self.assertIn((8, 9), self.universe._topology['bonds'])

    def test_tabulated_bond_nc(self):
        """Test for tabulated bonds without exclusion (TABBONDSNC)"""
        self.assertIn((9, 10), self.universe._topology['bonds'])

    def test_constraints(self):
        """Test for constraints (CONSTR)"""
        self.assertIn((10, 11), self.universe._topology['bonds'])
 
    def test_constraints_nc(self):
        """Test for constraints without connection (CONSTRNC)"""
        self.assertIn((11, 12), self.universe._topology['bonds'])

    def test_regular_angles(self):
        """Test for regular angles (ANGLES)"""
        self.assertIn((0, 1, 2), self.universe._topology['angles'])

    def test_G96_angles(self):
        """Test for G96 angles (G96ANGLES)"""
        self.assertIn((1, 2, 3), self.universe._topology['angles'])

    def test_cross_bond_bond(self):
        """Test for cross bond-bond angles (CROSS_BOND_BOND)"""
        self.assertIn((2, 3, 4), self.universe._topology['angles'])

    def test_cross_bond_angle(self):
        """Test for cross bond-angle angles (CROSS_BOND_ANGLE)"""
        self.assertIn((3, 4, 5), self.universe._topology['angles'])
    
    def test_urey_bradley(self):
        """Test for Urey-Bradley (UREY_BRADLEY)"""
        self.assertIn((4, 5, 6), self.universe._topology['angles'])

    def test_quartic_angles(self):
        """Test for quartic angles (QANGLES)"""
        self.assertIn((5, 6, 7), self.universe._topology['angles'])

    def test_restricted_bending_potential(self):
        """Test for restricted bending potentials (RESTRANGLES)"""
        self.assertIn((6, 7, 8), self.universe._topology['angles'])

    def test_tabulated_angle(self):
        """Test for tabulated angle (TABANGLES)"""
        self.assertIn((7, 8, 9), self.universe._topology['angles'])
