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
from numpy.testing import TestCase, dec, assert_
import functools

from MDAnalysis.tests.datafiles import (
    TPR,
    TPR400, TPR402, TPR403, TPR404, TPR405, TPR406, TPR407,
    TPR450, TPR451, TPR452, TPR453, TPR454, TPR455, TPR455Double,
    TPR460, TPR461, TPR502, TPR504, TPR505, TPR510, TPR2016,
    TPR510_bonded, TPR2016_bonded,
)

from MDAnalysisTests.topology.base import ParserBase
import MDAnalysis.topology.TPRParser


class TestTPR(ParserBase):
    """
    this test the data/adk_oplsaa.tpr which is of tpx version 58
    """
    filename = TPR
    parser = MDAnalysis.topology.TPRParser.TPRParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames']
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 3


# The follow test the same system grompped by different version of gromacs
# FORMAT: TPRABC, where numbers ABC indicates the version of gromacs that
# generates the corresponding tpr file

class TPRBase(ParserBase):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames']
    expected_n_atoms = 2263
    expected_n_residues = 230
    expected_n_segments = 2


# All these classes should be generated in a loop. Yet, nose test generation
# seems to work only with functions, and not with classes.
class TestTPR400(TPRBase):
    filename = TPR400

class TestTPR402(TPRBase):
    filename = TPR402

class TestTPR403(TPRBase):
    filename = TPR403

class TestTPR404(TPRBase):
    filename = TPR404

class TestTPR405(TPRBase):
    filename = TPR405

class TestTPR406(TPRBase):
    filename = TPR406

class TestTPR450(TPRBase):
    filename = TPR450

class TestTPR451(TPRBase):
    filename = TPR451

class TestTPR452(TPRBase):
    filename = TPR452

class TestTPR453(TPRBase):
    filename = TPR453

class TestTPR454(TPRBase):
    filename = TPR454

class TestTPR455(TPRBase):
    filename = TPR455

class TPRDouble(ParserBase):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames']
    expected_n_atoms = 21692
    expected_n_residues = 4352
    expected_n_segments = 7

class TestTPR455Double(TPRDouble):
    filename = TPR455Double


class TPR46xBase(ParserBase):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames']
    expected_n_atoms = 44052
    expected_n_residues = 10712
    expected_n_segments = 8


class TestTPR460(TPR46xBase):
    filename = TPR460

class TestTPR461(TPR46xBase):
    filename = TPR461


class TestTPR502(TPRBase):
    filename = TPR502

class TestTPR504(TPRBase):
    filename = TPR504

class TestTPR505(TPRBase):
    filename = TPR505

class TestTPR510(TPRBase):
    filename = TPR510

class TestTPR2016(TPRBase):
    filename = TPR2016


def _test_is_in_topology(name, elements, topology_path, topology_section):
    """
    Test if an interaction appears as expected in the topology
    """
    universe = MDAnalysis.Universe(topology_path)
    parser = MDAnalysis.topology.TPRParser.TPRParser(topology_path)
    top = parser.parse()
    for element in elements:
        assert_(element in getattr(top, topology_section).values,
                'Interaction type "{}" not found'.format(name))


def test_all_bonds():
    topologies = (TPR510_bonded, TPR2016_bonded)
    bonds = {'BONDS':[(0, 1)], 'G96BONDS':[(1, 2)], 'MORSE':[(2, 3)],
             'CUBICBONDS':[(3, 4)], 'CONNBONDS':[(4, 5)], 'HARMONIC':[(5, 6)],
             'FENEBONDS':[(6, 7)], 'RESTRAINTPOT':[(7, 8)],
             'TABBONDS':[(8, 9)], 'TABBONDSNC':[(9, 10)],
             'CONSTR':[(10, 11)], 'CONSTRNC':[(11, 12)],}
    bond_type_in_topology = functools.partial(_test_is_in_topology,
                                              topology_section='bonds')
    for topology in topologies:
        for bond_type, elements in bonds.items():
            yield (bond_type_in_topology, bond_type, elements, topology)


def test_all_angles():
    topologies = (TPR510_bonded, TPR2016_bonded)
    angles = {'ANGLES':[(0, 1, 2)], 'G96ANGLES':[(1, 2, 3)],
              'CROSS_BOND_BOND':[(2, 3, 4)], 'CROSS_BOND_ANGLE':[(3, 4, 5)],
              'UREY_BRADLEY':[(4, 5, 6)], 'QANGLES':[(5, 6, 7)],
              'RESTRANGLES':[(6, 7, 8)], 'TABANGLES':[(7, 8, 9)],}
    angle_type_in_topology = functools.partial(_test_is_in_topology,
                                               topology_section='angles')
    for topology in topologies:
        for angle_type, elements in angles.items():
            yield (angle_type_in_topology, angle_type, elements, topology)


def test_all_dihedrals():
    topologies = (TPR510_bonded, TPR2016_bonded)
    dihs = {'PDIHS':[(0, 1, 2, 3), (1, 2, 3, 4), (7, 8, 9, 10)],
            'RBDIHS':[(4, 5, 6, 7)], 'RESTRDIHS':[(8, 9, 10, 11)],
            'CBTDIHS':[(9, 10, 11, 12)], 'FOURDIHS':[(6, 7, 8, 9)],
            'TABDIHS':[(10, 11, 12, 13)],}
    dih_type_in_topology = functools.partial(_test_is_in_topology,
                                             topology_section='dihedrals')
    for topology in topologies:
        for dih_type, elements in dihs.items():
            yield (dih_type_in_topology, dih_type, elements, topology)


def test_all_impropers():
    topologies = (TPR510_bonded, TPR2016_bonded)
    imprs = {'IDIHS':[(2, 3, 4, 5), (3, 4, 5, 6)], 'PIDIHS':[(5, 6, 7, 8)]}
    impr_type_in_topology = functools.partial(_test_is_in_topology,
                                              topology_section='impropers')
    for topology in topologies:
        for impr_type, elements in imprs.items():
            yield (impr_type_in_topology, impr_type, elements, topology)

