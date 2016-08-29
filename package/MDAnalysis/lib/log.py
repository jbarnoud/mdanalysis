# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Setting up logging --- :mod:`MDAnalysis.lib.log`
====================================================

Configure logging for MDAnalysis. Import this module if logging is
desired in application code.

Logging to a file and the console is set up by default as described
under `logging to multiple destinations`_.

The top level logger of the library is named *MDAnalysis* by
convention; a simple logger that writes to the console and logfile can
be created with the :func:`create` function. This only has to be done
*once*. For convenience, the default MDAnalysis logger can be created
with :func:`MDAnalysis.start_logging`::

 import MDAnalysis
 MDAnalysis.start_logging()

Once this has been done, MDAnalysis will write messages to the logfile
(named `MDAnalysis.log` by default but this can be changed with the
optional argument to :func:`~MDAnalysis.start_logging`).

Any code can log to the MDAnalysis logger by using ::

 import logging
 logger = logging.getLogger('MDAnalysis.MODULENAME')

 # use the logger, for example at info level:
 logger.info("Starting task ...")

The important point is that the name of the logger begins with
"MDAnalysis.".

.. _logging to multiple destinations:
   http://docs.python.org/library/logging.html?#logging-to-multiple-destinations

.. SeeAlso:: The :mod:`logging` module in the standard library contains
             in depth documentation about using logging.


Convenience functions
---------------------

Two convenience functions at the top level make it easy to start and
stop the default *MDAnalysis* logger.

.. autofunction:: MDAnalysis.start_logging
.. autofunction:: MDAnalysis.stop_logging


Other functions and classes for logging purposes
------------------------------------------------

.. autogenerated, see Online Docs
"""
from __future__ import print_function, division

import datetime
import sys
import logging
import re
import warnings

from .. import version


def start_logging(logfile="MDAnalysis.log", version=version.__version__):
    """Start logging of messages to file and console.

    The default logfile is named `MDAnalysis.log` and messages are
    logged with the tag *MDAnalysis*.
    """
    create("MDAnalysis", logfile=logfile)
    logging.getLogger("MDAnalysis").info(
        "MDAnalysis %s STARTED logging to %r", version, logfile)


def stop_logging():
    """Stop logging to logfile and console."""
    logger = logging.getLogger("MDAnalysis")
    logger.info("MDAnalysis STOPPED logging")
    clear_handlers(logger)  # this _should_ do the job...


def create(logger_name="MDAnalysis", logfile="MDAnalysis.log"):
    """Create a top level logger.

    - The file logger logs everything (including DEBUG).
    - The console logger only logs INFO and above.

    Logging to a file and the console as described under `logging to
    multiple destinations`_.

    The top level logger of MDAnalysis is named *MDAnalysis*.  Note
    that we are configuring this logger with console output. If a root
    logger also does this then we will get two output lines to the
    console.

    .. _logging to multiple destinations:
       http://docs.python.org/library/logging.html?#logging-to-multiple-destinations
    """

    logger = logging.getLogger(logger_name)

    logger.setLevel(logging.DEBUG)

    # handler that writes to logfile
    logfile_handler = logging.FileHandler(logfile)
    logfile_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logfile_handler.setFormatter(logfile_formatter)
    logger.addHandler(logfile_handler)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def clear_handlers(logger):
    """clean out handlers in the library top level logger

    (only important for reload/debug cycles...)

    """
    for h in logger.handlers:
        logger.removeHandler(h)


class NullHandler(logging.Handler):
    """Silent Handler.

    Useful as a default::

      h = NullHandler()
      logging.getLogger("MDAnalysis").addHandler(h)
      del h

    see the advice on logging and libraries in
    http://docs.python.org/library/logging.html?#configuring-logging-for-a-library
    """
    def emit(self, record):
        pass


def echo(s='', replace=False, newline=True):
    """Simple string output that immediately prints to the console.

    Parameters
    ==========

    s: str
        The string to output.
    replace: bool
        If ``True``, the string will be printed on top of the current line. In
        practice, ``\r`` is added at the beginning og the string.
    newline: bool
        If ``True``, a new line id added at the end of the string.
    """
    if replace:
        s = '\r' + s
    if newline:
        end='\n'
    else:
        end=''
    print(s, file=sys.stderr, end=end)
    sys.stderr.flush()


def _new_format(template, variables):
    """Format a string that follows the {}-based syntax.
    """
    return template.format(**variables)


def _legacy_format(template, variables):
    """Format a string that follows the %-based syntax.
    """
    return template % variables


def _guess_string_format(template):
    """Guess if the template follow {}-based or %-based syntax.
    """
    match = re.search(r'%\((step|numsteps|percentage)\)', template)
    if match is None:
        return _new_format
    else:
        return _legacy_format


def _set_verbose(verbose, quiet, default=True,
                 was='quiet', now='verbose'):
    """Return the expected value of verbosity

    This function aims at handling the deprecation of the *quiet* keyword in
    versin 0.16.

    This function issues a deprecation warning if *quiet* was set (is not
    None), and raises a ValueError if *verbose* and *quiet* are set to
    contradicting values.

    If *verbose* is set, then the function returns the set value of *verbose*.
    If it is not set, but *quiet* is set, then the function returns
    `not quiet`. Finally, if none of *verbose* nor *quiet* is set, then
    *default* is returned.

    During the deprecation phase of the *quiet* keyword, this function is
    expected to be used as follow:

    .. code-block:: python

       def method(verbose=None, quiet=None):
           # *verbose* and *quiet* are set to None to distinguish explicitly
           # set values.
           self.verbose = _set_verbose(verbose, quiet, default=True)

    At the end of the deprecation period, the code above should be replaced by:

    .. code-block:: python

       def method(verbose=True):
           # The *quiet* keyword disapeard and the default value for *verbose*
           # is set to the actual default value.
           self.verbose = verbose

    In `MDAnalysis.analysis.hbonds.hbonds_analysis`, the deprecation scheme is
    more complex: *quiet* becomes *verbose*, and *verbose* becomes *debug*.
    Hence, this function allows to use diffrent argument names to display in
    error messages and deprecation warnings.
    """
    if quiet is not None:
        warnings.warn("Keyword *{}* is deprecated (from version 0.16); "
                      "use *{}* instead.".format(was, now), DeprecationWarning)
        if verbose is not None and verbose == quiet:
            raise ValueError("Keywords *{}* and *{}* are contradicting each other."
                             .format(now, was))
        return not quiet
    elif verbose is None:
        return default
    else:
        return verbose


class ProgressMeter(object):
    r"""Simple progress meter

    Usage::

       u = Universe(PSF, DCD)
       pm = ProgressMeter(u.trajectory.n_frames, interval=100)
       for ts in u.trajectory:
           pm.echo(ts.frame)
           ...

    For a trajectory with 10000 frames this will produce output such
    as ::

       Step   100/10000 [  1.0%]
       Step   200/10000 [  2.0%]
       ...

    The default *format* string is::

        Step {step:5d}/{numsteps} [{percentage:5.1f}%]

    By default, each line of the progress meter is displayed on top of the
    previous one. To prevent this behaviour, set the *dynamic* keyword to
    ``False``.

    It is possible to embed (almost) arbitrary additional data in the
    format string, for example a current RMSD value::

        format_line = "RMSD {rmsd:5.2f} at {step:5d}/{numsteps} [{percentage:5.1f}%]"
        pm = ProgressMeter(u.trajectory.n_frames,
                           interval=100, format=format_line)
        for ts in u.trajectory:
           pm.echo(ts.frame, rmsd=current_rmsd)
           ...

    This will print something like::

       RMSD   1.02 at  100/10000 [  1.0%]
       RMSD   1.89 at  200/10000 [  2.0%]
       ...

    """

    def __init__(self, numsteps, format=None, interval=10, offset=1,
                 verbose=None, dynamic=True,
                 format_handling='auto', quiet=None):
        r"""
        Parameters
        ==========

        numsteps: int
            total number of steps
        interval: int
            only calculate progress every *interval* steps [10]
        format: str
            a format string with Python variable interpolation. Allowed
            values:

            * *step*: current step
            * *numsteps*: total number of steps as supplied in *numsteps*
            * *percentage*: percentage of total
            * *elapsed*: time elapsed since the first update
            * *time_per_iteration*: average time between updates
            * *estimate_time_to_completion* or *etc*: estimate time to
              completion

            Time tracking and estimate time to completion are not available
            with {}-based format syntax. *estimate_time_to_completion* is a
            :class:`datetime.datetime` instance, *elapsed* and
            *time_per_iteration* are :class:`datetime.timedelta` instances;
            they can be formated as such.

            The last call to :meth:`ProgressMeter.echo` will automatically
            issue a newline ``\n``.

            If *format* is ``None`` then the default is used::

                Step {step:5d}/{numsteps} [{percentage:5.1f}%]

        offset: int
            number to add to *step*; e.g. if *step* is 0-based (as in MDAnalysis)
            then one should set *offset* = 1; for 1-based steps, choose 0. [1]
        verbose: bool
            If ``False``, disable all output, ``True`` print all messages as
            specified, [``True``]
        dynamic: bool
            If ``True``, each line will be printed on top of the previous one.
            This is done by prepedind the format with ``\r``. [``True``]
        format_handling: str
            how to handle the format string. Allowed values are:

            * *new*: the format string uses {}-based formating
            * *legacy*: the format string uses %-basedd formating
            * *auto*: default, try to guess how the format string should be
              handled

            When using the *legacy*, %-based, formated syntax, the time
            tracking and estimate time to completion are not available.

        .. versionchanged:: 0.8
           Keyword argument *quiet* was added.

        .. deprecated:: 0.16
           Keyword argument *quiet* is deprecated in favor of *verbose*.

        .. versionchanged:: 0.16
           Keyword argument *dynamic* replaces ``\r`` in the format.
           Keyword argument *dynamic* replaces ``\r`` in the format. The
           {}-based formating syntax is recommended for the *format* argument,
           and some time tracking (*elapsed*, *time_per_iteration*, and
           *estimate_time_to_completion* template variable and attributes) are
           added.
        """
        self.numsteps = numsteps
        self.interval = int(interval)
        self.offset = int(offset)
        self.dynamic = dynamic
        self.numouts = -1

        # The *quiet* keyword argument is deprecated.
        self.verbose = _set_verbose(verbose, quiet, default=True)

        if format is None:
            format = "Step {step:5d}/{numsteps} [{percentage:5.1f}%]"
            self.format_handler = _new_format
        else:
            if format_handling == 'auto':
                self.format_handler = _guess_string_format(format)
            else:
                self.format_handler = {'new': _new_format,
                                       'legacy': _legacy_format}[format_handling]
        self.format = format
        self.step = 0
        self.percentage = 0.0
        assert numsteps > 0, "numsteps step must be >0"
        assert interval > 0, "interval step must be >0"

        # Time tracking stuff
        self._start_time = None
        self._last_time = None
        self.elapsed = None
        self.time_per_iteration = None
        self.estimate_time_to_completion = None
        self.etc = None

    def update(self, step, **kwargs):
        """Update the state of the ProgressMeter.

        *kwargs* are additional attributes that can be references in
        the format string.
        """
        self.step = step + self.offset
        self.percentage = 100. * self.step / self.numsteps
        for k, v in kwargs.items():
            setattr(self, k, v)

        self.numouts += 1

        # Update time tracking and estimate time to completion
        self._last_time = datetime.datetime.now()
        if self._start_time is None:
            self._start_time = self._last_time
        self.elapsed = self._last_time - self._start_time
        self.time_per_iteration = datetime.timedelta(
            seconds=self.elapsed.total_seconds() / self.step)
        remaining = self.numsteps - self.step
        self.estimate_time_to_completion = (self._last_time
                                            + remaining
                                            * self.time_per_iteration)
        self.etc = self.estimate_time_to_completion
       
    def echo(self, step, **kwargs):
        """Print the state to stderr, but only every *interval* steps.

        1) calls :meth:`~ProgressMeter.update`
        2) writes step and percentage to stderr with :func:`echo`,
           using the format string (in :attr:`ProgressMeter.format`)

        The last step is always shown, even if not on an *interval*, and a
        carriage return is replaced with a new line for a cleaner display.

        *kwargs* are additional attributes that can be references in
        the format string.

        .. Note:: If *verbose* = ``False`` has been set in the
                  constructor or if :attr:`ProgressMeter.verbose` has
                  been set to ``False``, then no messages will be
                  printed.
        """
        if not self.verbose:
            return
        self.update(step, **kwargs)
        format = self.format
        newline = not self.dynamic
        if self.step == self.numsteps:
            newline = True
        elif self.numouts % self.interval == 0:
            pass
        else:
            return
        echo(self.format_handler(format, vars(self)),
             replace=self.dynamic, newline=newline)
