# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


from . import _AUXREADERS
from ..lib import util

def get_auxreader_for(auxdata, format=None):
    if format is None:
        if type(auxdata) == str:
            ## assume it's a filename?
            format = util.guess_format(auxdata)
        else:
            ## arrays etc
            pass
    format = format.upper()
    try:
        return _AUXREADERS[format]
    except KeyError:
        raise ValueError("Unknown auxiliary data format for {0}".format(auxdata))
