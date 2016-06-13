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
