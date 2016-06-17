#!/usr/bin/env python

#
# Copyright (C) 2010, 2016  Smithsonian Astrophysical Observatory
#
#
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""

"""

TOOLNAME = "correct_periscope_drift"
VERSION = "0.1"

# import standard python modules as required
import os
import sys
import numpy as np
from numpy import sin, cos, tan, arctan2, radians, degrees, sqrt



# Import the CIAO contributed modules.
# 
# Here I explicitly name what is imported for documentation
# purposes. You can use the 'from foo import *' syntax but
# in general it is better to either be explicit or load
# into a separate namespace, so that you can find out
# where a routine or symbol is defined
#


from ciao_contrib.logger_wrapper import initialize_logger, make_verbose_level, set_verbosity, handle_ciao_errors
from ciao_contrib.param_wrapper import open_param_file
from ciao_contrib.runtool import dmstat
import pycrates
from pychips import (add_curve, print_window, set_plot_xlabel, set_plot_ylabel, clear_plot,
                     add_window, set_plot_title)
from sherpa import ui


# Set up the logging/verbose code
initialize_logger(TOOLNAME)

# Use v<n> to display messages at the given verbose level.
# You can pick other names than v0,v1, ... v5 if desired.
#
v1 = make_verbose_level(TOOLNAME, 1)
v2 = make_verbose_level(TOOLNAME, 2)
v5 = make_verbose_level(TOOLNAME, 5)


def process_command_line(argv):

    import paramio as pio
    """Handle the parameter input for this script."""

    if argv is None or argv == []:
        raise ValueError("argv argument is None or empty")

    # open_param_file, from ciao_contrib.param_wrapper,
    # opens the parameter file and sorts out the command
    # line. It returns a dictionary with a few fields;
    # fp is the handle to use in paramio calls, as shown
    # below; and parnams is the name to use if you need
    # to re-open the parameter file to set things
    # (which is unlikely, speak to Doug if you think you
    # might need this).
    #
    pinfo = open_param_file(argv, toolname=TOOLNAME)
    fp = pinfo["fp"]

    # Use the parameter library to get the arguments
    # and perform any desired validity checks (e.g.
    # the outfile check/clobber check shown below).
    #
    # Parameters should be queried in the same order
    # as the parameter file
    #
    mypars = {'progname': pinfo['progname'],
              'parname': pinfo['parname']}

    for stringpar in ['evtfile', 'input_asolfile',
                      'corr_asolfile', 'corr_plot_prefix']:
        mypars[stringpar] = pio.pgetstr(fp, stringpar)
        if mypars[stringpar].strip() == "":
            raise ValueError("{} parameter is empty".format(stringpar))

    mypars['corr_poly_degree'] = pio.pgeti(fp, "corr_poly_degree")

    for floatpar in ['x', 'y', 'radius', 'src_min_counts']:
        mypars[floatpar] = pio.pgetd(fp, floatpar)

    clobber = pio.pgetb(fp, "clobber")
    verbose = pio.pgeti(fp, "verbose")

    # We close the parameter file here; if you need to write
    # values to the file you can either leave it open or -
    # possibly better - is to close it and re-open it later.
    #
    pio.paramclose(fp)

    # Set tool and module verbosity
    set_verbosity(verbose)

    # check outfile, exiting if it exists and clobber is
    # not set. Note: the error message should be updated to better
    # match that used by CIAO tools but I have not done this.
    #
    #
    if not clobber and os.path.exists(mypars['corr_asolfile']):
        raise IOError("clobber is no and outfile ({0}) exists".format(
                mypars['corr_asolfile']))

    mypars['clobber'] = clobber
    mypars['verbose'] = verbose

    # Return a dictionary with useful info (you could return a tuple or
    # something else, depending on the parameter file)
    #
    return mypars


# Display parameter info to the user.
# The format and choice of information to display is up to you;
# I've chosen to include some ancilalry information but it
# depends on the tool what should be used.
#
def display_start_info(opts):
    v1("Running: {0}".format(opts["progname"]))
    v2("  version = {0}".format(VERSION))
    v2("with parameters:")
    v2("  evtfile={0}".format(opts["evtfile"]))
    v2("  input_asolfile={0}".format(opts["input_asolfile"]))
    # probably other values here too
    v2("  verbose={0}".format(opts["verbose"]))
    v2("  and CALDB is set to  {0}".format(os.environ["CALDB"]))
    v2("  and ASCDS_INSTALL is {0}".format(os.environ["ASCDS_INSTALL"]))
    v2("-" * 60)


def equatorial2transform(ra, dec, roll):
    """Construct the transform/rotation matrix from RA,Dec,Roll (given in degrees)
    :returns: transform matrix
    :rtype: Nx3x3 numpy array

    """
    ra = np.radians(ra)
    dec = np.radians(dec)
    roll = np.radians(roll)
    ca = np.cos(ra)
    sa = np.sin(ra)
    cd = np.cos(dec)
    sd = np.sin(dec)
    cr = np.cos(roll)
    sr = np.sin(roll)
    # This is the transpose of the transformation matrix (related to
    # translation of original perl code
    rmat = np.array(
        [[ca * cd,                  sa * cd,                sd     ],
         [-ca * sd * sr - sa * cr, -sa * sd * sr + ca * cr, cd * sr],
         [-ca * sd * cr + sa * sr, -sa * sd * cr - ca * sr, cd * cr]])

    return rmat.transpose()


def radec2eci(ra, dec):
    """
    Convert from RA,Dec to ECI.  The input ``ra`` and ``dec`` values can be 1-d
    arrays of length N in which case the output ``ECI`` will be an array with
    shape (3,N).

    Borrowed from Ska.quatutil

    :param ra: Right Ascension (degrees)
    :param dec: Declination (degrees)
    :returns: numpy array ECI (3-vector or 3xN array)
    """
    r = np.radians(ra)
    d = np.radians(dec)
    return np.array([np.cos(r) * np.cos(d), np.sin(r) * np.cos(d), np.sin(d)])


#def roll_range(asp_roll):
#    max_roll = np.max(asp_roll)
#    min_roll = np.min(asp_roll)
#    if (max_roll - min_roll) > 180:
#        roll_180 = asp_roll.copy()
#        roll_180[roll_180 > 180] -= 360
#        max_roll = np.max(roll_180)
#        min_roll = np.min(roll_180)
#        if max_roll - min_roll > 180:
#            raise ValueError("> 180 deg roll range")
#    return max_roll - min_roll
#

def extract_events(event_file, src_x, src_y, src_radius):
    """
    Get events from specified source circle

    :param event_file: Chandra event 1 or 2 file
    :param src_x: Sky X coordinate of source region to extract
    :param src_y: Sky Y coordinate of source region to extract
    :param src_radius: Source region/circle radius in pixels
    :returns: CRATE of events
    """
    import pycrates
    regstring = "circle({},{},{})".format(src_x, src_y, src_radius)
    events = pycrates.read_file("{}[sky={}]".format(
            event_file, regstring))
    return events


def get_event_yag_zag(evt_ra, evt_dec, ra_nom, dec_nom, roll_nom):
    """
    Convert RA and Dec positions into Y and Z angle

    This takes the RA and Dec positions of the events and takes a "nominal" pointing
    (which maybe taken from the nominal values provided in the event file header
    and uses that reference nominal pointing to convert to Y and Z angle
    that should relate to Aspect Camera Y and Z angle for the purposes of
    fitting Y and Z angle periscope drift.

    :param evt_ra: event RA
    :param evt_dec: event Dec
    :param ra_nom: A single "nominal" reference value for RA
    :param dec_nom: A single "nominal" reference value for Dec
    :param roll_nom: A single "nominal" reference value for Roll
    :returns: yag, zag in arcsecs
    """
    if len(evt_ra) != len(evt_dec):
        raise ValueError("len(evt_ra) != len(evt_dec), {} != {}".format(
                len(evt_ra), len(evt_dec)))

    # Transform to Earth Centered Inertial
    eci = radec2eci(evt_ra, evt_dec)
    # transform from 3 x N to N x 3
    eci = eci.transpose()

    att_stack = np.repeat(np.array([ra_nom, dec_nom, roll_nom]),
                          len(evt_ra)).reshape(3, len(evt_ra))
    # Transforms
    Ts = equatorial2transform(att_stack[0], att_stack[1], att_stack[2])

    # The position of the events rotated into the frame of the
    d_aca = np.sum(Ts.transpose(2, 0, 1) * eci, axis=-1).transpose()

    R2A = 206264.81
    yag = np.arctan2(d_aca[:, 1], d_aca[:, 0]) * R2A
    zag = np.arctan2(d_aca[:, 2], d_aca[:, 0]) * R2A
    return yag, zag


def time_bins(times, x, nbins=20):
    h, bins = np.histogram(times, bins=nbins)
    bin_centers = (bins[:-1] + bins[1:]) / 2.0
    inds = np.digitize(times, bins) - 1
    bin_x = []
    bin_std = []
    for idx in range(0, nbins):
        data = x[inds == idx]

        bin_x.append(np.mean(data))
        bin_std.append(np.std(data)/np.sqrt(len(data)))
    return bin_centers, np.array(bin_x), np.array(bin_std)


def fit(fit_data, evt_times, data_id, opt):

    init_error = 5

    ui.clean()
    ui.load_arrays(data_id, evt_times - evt_times[0], fit_data,
                   np.zeros_like(fit_data) + init_error)
    v2("Fitting a line to the data to get reduced stat errors")
    # First just fit a line to get reduced errors on this set
    ui.polynom1d.line
    ui.set_model(data_id, 'line')
    ui.thaw('line.c1')
    ui.fit(data_id)
    fit = ui.get_fit_results()
    calc_error = init_error * np.sqrt(fit.rstat)
    ui.set_staterror(data_id, calc_error)
    # Then fit the specified model
    v2("Fitting a polynomial of degree {} to the data".format(opt['corr_poly_degree']))
    ui.polynom1d.fitpoly
    ui.freeze('fitpoly')
    # Thaw the coefficients requested by the degree of the desired polynomial
    for deg in range(1, 1 + opt['corr_poly_degree']):
        ui.thaw("fitpoly.c{}".format(deg))
    if opt['corr_poly_degree'] > 1:
        ui.thaw('fitpoly.offset')
    ui.set_model(data_id, 'fitpoly')
    ui.fit(data_id)
    mp = ui.get_model_plot(data_id)
    return mp


# The '@handle_ciao_errors' decorator will catch any error thrown and
# display it in a format like that of a CIAO tool, then exit with a
# non-zero status. So, if at any point in your code you need to exit just
# raise an error and let handle_ciao_errors bother with the display.
#
# This also means that routines can be used from ChIPS/Sherpa/other scripts
# without having to worry about them calling sys.exit.


@handle_ciao_errors(TOOLNAME, VERSION)
def main(opt):
    events = extract_events(opt['evtfile'],
                            opt['x'], opt['y'], opt['radius'])

    evt_ra_nom = events.get_key('RA_NOM').value
    evt_dec_nom = events.get_key('DEC_NOM').value
    evt_roll_nom = events.get_key('ROLL_NOM').value

    asol = pycrates.read_file(opt['input_asolfile'])
    asol_times = asol.get_column('time').values

    # Sanity check the two input files
    asol_obsid = asol.get_key('OBS_ID').value
    evt_obsid = events.get_key('OBS_ID').value
    if asol_obsid != evt_obsid:
        v1("Error Aspect solution obsid {} != event file obsid {}".format(asol_obsid, evt_obsid))

    evt_ra = events.get_column('RA').values
    evt_dec = events.get_column('Dec').values
    evt_times = events.get_column('Time').values

    # Limit to only using events contained within the range of the aspect solution
    ok_times = (evt_times > asol_times[0]) & (evt_times < asol_times[-1])
    if not np.any(ok_times):
        raise ValueError("No events in region are contained within time range of aspect solution.")
    # Limit this *in place*
    evt_ra = evt_ra[ok_times]
    evt_dec = evt_dec[ok_times]
    evt_times = evt_times[ok_times]

    if len(evt_times) < opt['src_min_counts']:
        v1("Warning only {} counts in src region.  {} minimum suggested 'src_min_counts'".format(
                len(evt_times), opt['src_min_counts']))

    ax_data = {}
    ax_map = {'yag': 'dy',
              'zag': 'dz'}

    ax_data['yag'], ax_data['zag'] = get_event_yag_zag(evt_ra, evt_dec,
                                                       evt_ra_nom, evt_dec_nom, evt_roll_nom)

    # Store comments to print in block after all of the sherpa fit output
    fit_comments = []
    plot_list = []

    for data_id, ax in enumerate(['yag', 'zag']):
        fit_data = ax_data[ax] - np.mean(ax_data[ax])
        mp = fit(fit_data, evt_times, data_id, opt)

        bin_centers, bin_mean, bin_std = time_bins(evt_times, fit_data)

        add_window()
        add_curve((bin_centers - evt_times[0]) / 1000., bin_mean, [bin_std, +bin_std],
                  ["line.style", "none", "symbol.style", "none", "err.style", "cap"])
        add_curve(mp.x / 1000., mp.y, ["symbol.style", "none"])
        set_plot_xlabel("Observation elapsed/delta time (ks)")
        set_plot_ylabel("Position offset from mean, {} (arcsec)".format(ax))
        set_plot_title("Fit of {} data (with time-binned event offsets)".format(ax))
        fit_plot = "{}_fit_{}.png".format(opt['corr_plot_prefix'], ax)
        if os.path.exists(fit_plot) and opt['clobber']:
            os.unlink(fit_plot)
        plot_list.append(fit_plot)
        print_window(fit_plot)

        add_window()
        data_plot = "{}_data_{}.png".format(opt['corr_plot_prefix'], ax)
        ui.plot_fit(data_id)
        if os.path.exists(data_plot) and opt['clobber']:
            os.unlink(data_plot)
        set_plot_xlabel("Observation elapsed/delta time (ks)")
        set_plot_ylabel("Position offset from mean, {} (arcsec)".format(ax))
        set_plot_title("Raw data and fit in {}".format(ax))
        plot_list.append(data_plot)
        print_window(data_plot)

        asol_corr = np.interp(asol_times, mp.x + evt_times[0], mp.y)
        asol_col_to_fix = asol.get_column(ax_map[ax])
        fit_comments.append("Events show drift range of {:.2f} arcsec in {} axis".format(
                np.max(asol_corr) - np.min(asol_corr), ax))
        fit_comments.append("Max absolute correction of {:.2f} arcsec for {} axis".format(
                np.max(np.abs(asol_corr)), ax))

        asol_col_to_fix.values += (asol_corr / 20)

    v1("-" * 60)
    v1("Fit results")
    for c in fit_comments:
        v1("\t{}".format(c))
    v1("-" * 60)
    v2("Writing out corrected aspect solution file to {}".format(opt['corr_asolfile']))
    v2("\tTo review fit see correction plots in:")
    for p in plot_list:
        v2("\t\t{}".format(p))

    # Actually write out the new aspect solution file
    asol.write(opt['corr_asolfile'], clobber=opt['clobber'])


if __name__ == "__main__":
    opt = process_command_line(sys.argv)
    display_start_info(opt)
    main(opt)





