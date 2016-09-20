# Tests to be run from ska env not ciao
import os
from glob import glob
from astropy.table import Table
from Ska.astro import sph_dist
import Ska.Shell
from Ska.Shell import bash
import Ska.File
import numpy as np
import pytest
import tempfile


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


def get_event_yag_zag(evt_ra, evt_dec, ra_pnt, dec_pnt, roll_pnt):
    """
    Convert RA and Dec positions into Y and Z angle

    This takes the RA and Dec positions of the events and takes a reference pointing
    (which maybe taken from the NOM or PNT values provided in the event file header)
    and uses that reference pointing to convert to Y and Z angle
    that should relate to Aspect Camera Y and Z angle for the purposes of
    fitting Y and Z angle periscope drift.

    :param evt_ra: event RA
    :param evt_dec: event Dec
    :param ra_pnt: A single reference value for RA
    :param dec_pnt: A single reference value for Dec
    :param roll_pnt: A single reference value for Roll
    :returns: yag, zag in arcsecs
    """
    if len(evt_ra) != len(evt_dec):
        raise ValueError("len(evt_ra) != len(evt_dec), {} != {}".format(
                len(evt_ra), len(evt_dec)))

    # Transform to Earth Centered Inertial
    eci = radec2eci(evt_ra, evt_dec)
    # transform from 3 x N to N x 3
    eci = eci.transpose()

    att_stack = np.repeat(np.array([ra_pnt, dec_pnt, roll_pnt]),
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
    """
    Bin 'x' by the times in 'times'.

    :param times: times used for time bins
    :param x: dataset binned in equal time chunks
    :param n_bins: number of time bins to use
    :returns: bin time centers, bin data mean, bin data std
    """

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


def p2p(times, ras, decs, ra_pnt, dec_pnt, roll_pnt):
    """
    Get the Y and Z peak-to-peak values for the binned mean of the data.  *_pnt values
    used for reference point for yag/zag transformation.
    """
    yag, zag = get_event_yag_zag(ras, decs, ra_pnt, dec_pnt, roll_pnt)
    bins_yag = time_bins(times, yag, nbins=40)
    bins_zag = time_bins(times, zag, nbins=40)
    return np.max(bins_yag[1]) - np.min(bins_yag[1]), np.max(bins_zag[1]) - np.min(bins_zag[1])


def max_offset(times, ras, decs, ra_pnt, dec_pnt, roll_pnt):
    """
    Get the Y and Z max offset values for the binned mean of the data.  *_pnt values
    used for reference point for yag/zag transformation.
    """
    yag, zag = get_event_yag_zag(ras, decs, ra_pnt, dec_pnt, roll_pnt)
    bins_yag = time_bins(times, yag)
    bins_zag = time_bins(times, zag)
    return np.max(np.abs(bins_yag[1] - np.mean(yag))), np.max(np.abs(bins_zag[1] - np.mean(zag)))


def cpd(asol, evtfile, x, y, radius, out):
    """
    Run the correct periscope drift tool.
    """
    bash("""./correct_periscope_drift infile= {asol} \
evtfile= {evtfile} \
x={x} y={y} radius={radius} \
outfile= {out}_asol1.fits \
corr_plot_root= {out} \
clobber+""".format(testdir=TESTDIR, evtfile=evtfile, asol=asol, y=y, x=x, radius=radius, out=out),
         env=ciaoenv)


ciaoenv = Ska.Shell.getenv('source /soft/ciao/bin/ciao.sh')

TESTDIR = os.path.abspath(".")

obsids = [{'obsid': 17128, 'x': 4064.88, 'y': 4071.87, 'radius': 6},
          {'obsid': 9926, 'x': 4103.6, 'y': 4062.3, 'radius': 6},
          {'obsid': 10228, 'x': 32888.17, 'y': 32710.78, 'radius': 20}]


@pytest.mark.parametrize("src_info", obsids)
def test_process(src_info):
    """
    Run the correct periscope tool on data and reprocess events.  Test that the outputs are somewhat
    reasonable.
    """
    tempdir = tempfile.mkdtemp()
    with Ska.File.chdir(tempdir):
        print os.path.abspath(os.curdir)
        # Work with a local copy of the tools.  Also, the parameters for correct periscope drift
        # don't work unless it is 'named' that, so use a symbolic link as needed.
        obsid = src_info['obsid']
        x = src_info['x']
        y = src_info['y']
        radius = src_info['radius']

        bash("cp {}/correct_periscope_drift.par .".format(TESTDIR))
        bash("ln -s {}/../correct_periscope_drift.py correct_periscope_drift".format(TESTDIR))

        bash("download_chandra_obsid {}".format(obsid), env=ciaoenv)

        asol = glob("{}/primary/*asol1*fits*".format(obsid))[0]
        evt2file = glob("{}/primary/*evt2*fits*".format(obsid))[0]
        evt1file = glob("{}/secondary/*evt1*fits*".format(obsid))[0]

        cpd(asol, evt2file, x, y, radius, 'cpd_from_evt2')
        cpd(asol, evt1file, x, y, radius, 'cpd_from_evt1')

        # reduce the evt1 file to just the source region in place
        bash("""dmcopy {evt1file}"[(x,y)=circle({x}, {y}, {radius})]" {obsid}_orig_src_evt1.fits clobber+""".format(evt1file=evt1file, x=x, y=y, radius=radius, obsid=obsid),
             env=ciaoenv)
        bash("gzip --fast -f {}_orig_src_evt1.fits".format(obsid))
        bash("cp {}_orig_src_evt1.fits.gz {}".format(obsid, evt1file))

        bash("mv {} {}/".format(asol, obsid))
        # reprocess with the aspect solution from cpd via evt2 fit
        bash("cp cpd_from_evt2_asol1.fits {}/primary/".format(obsid))
        bash("chandra_repro indir={}/ outdir={}/repro_evt2 verbose=1".format(obsid, obsid),
             env=ciaoenv)

        # reprocess with the aspect solution from cpd via evt1 fit
        bash("rm {}/primary/cpd_from_evt2_asol1.fits".format(obsid))
        bash("cp cpd_from_evt1_asol1.fits {}/primary/".format(obsid))
        bash("chandra_repro indir={}/ outdir={}/repro_evt1 verbose=1".format(obsid, obsid),
             env=ciaoenv)

        # Extract RA,Dec,Time from the event files
        bash('dmcopy {evt2file}"[(x,y)=circle({x},{y},{radius})][cols ra,dec,time]" {obsid}_orig_evt2_radectime.fits clobber+'.format(
                evt2file=evt2file, x=x, y=y, radius=radius, obsid=obsid), env=ciaoenv)
        evtfile_from_evt2 = glob("{}/repro_evt2/*evt2.fits*".format(obsid))[0]
        bash('dmcopy {evt2file}"[(x,y)=circle({x},{y},{radius})][cols ra,dec,time]" {obsid}_fixed_from_evt2_radectime.fits clobber+'.format(
                evt2file=evtfile_from_evt2, x=x, y=y, radius=radius, obsid=obsid), env=ciaoenv)
        evtfile_from_evt1 = glob("{}/repro_evt1/*evt2.fits*".format(obsid))[0]
        bash('dmcopy {evt2file}"[(x,y)=circle({x},{y},{radius})][cols ra,dec,time]" {obsid}_fixed_from_evt1_radectime.fits clobber+'.format(
                evt2file=evtfile_from_evt1, x=x, y=y, radius=radius, obsid=obsid), env=ciaoenv)

        # Read those files
        uncorr_data = Table.read("{}_orig_evt2_radectime.fits".format(obsid))
        evt1_corr_data = Table.read("{}_fixed_from_evt1_radectime.fits".format(obsid))
        evt2_corr_data = Table.read("{}_fixed_from_evt2_radectime.fits".format(obsid))
        # ignore any values of 0 that seem to be from dmcopy with ra,dec
        evt1_corr_data = evt1_corr_data[evt1_corr_data['RA'] != 0]
        evt2_corr_data = evt2_corr_data[evt2_corr_data['RA'] != 0]
        orig_aspect = Table.read(glob("{}/*asol1*".format(obsid))[0])
        ra_pnt = np.mean(orig_aspect['ra'])
        dec_pnt = np.mean(orig_aspect['dec'])
        roll_pnt = np.mean(orig_aspect['roll'])

        # Peak to peak of source data
        uncorr_p2p = p2p(uncorr_data['time'], uncorr_data['RA'], uncorr_data['DEC'],
                         ra_pnt, dec_pnt, roll_pnt)
        evt1_corr_p2p = p2p(evt1_corr_data['time'], evt1_corr_data['RA'], evt1_corr_data['DEC'],
                            ra_pnt, dec_pnt, roll_pnt)
        evt2_corr_p2p = p2p(evt2_corr_data['time'], evt2_corr_data['RA'], evt2_corr_data['DEC'],
                            ra_pnt, dec_pnt, roll_pnt)

        # yags are better or close in peak to peak
        assert evt1_corr_p2p[0] < uncorr_p2p[0] + 0.03
        assert evt2_corr_p2p[0] < uncorr_p2p[0] + 0.03
        # zags are better or close in peak to peak
        assert evt1_corr_p2p[1] < uncorr_p2p[1] + 0.03
        assert evt2_corr_p2p[1] < uncorr_p2p[1] + 0.03

        # Offsets from mean
        uncorr_offset = max_offset(uncorr_data['time'], uncorr_data['RA'], uncorr_data['DEC'],
                         ra_pnt, dec_pnt, roll_pnt)
        evt1_corr_offset = max_offset(evt1_corr_data['time'], evt1_corr_data['RA'], evt1_corr_data['DEC'],
                            ra_pnt, dec_pnt, roll_pnt)
        evt2_corr_offset = max_offset(evt2_corr_data['time'], evt2_corr_data['RA'], evt2_corr_data['DEC'],
                            ra_pnt, dec_pnt, roll_pnt)

        # yags are better or close
        assert evt1_corr_offset[0] < uncorr_offset[0] + 0.03
        assert evt2_corr_offset[0] < uncorr_offset[0] + 0.03
        # zags are better or close
        assert evt1_corr_offset[1] < uncorr_offset[1] + 0.03
        assert evt2_corr_offset[1] < uncorr_offset[1] + 0.03


        # means are unchanged
        evt1_off = 3600 * sph_dist(np.mean(uncorr_data['RA']), np.mean(uncorr_data['DEC']),
                                   np.mean(evt1_corr_data['RA']), np.mean(evt1_corr_data['DEC']))
        assert evt1_off < 0.1
        evt2_off = 3600 * sph_dist(np.mean(uncorr_data['RA']), np.mean(uncorr_data['DEC']),
                                   np.mean(evt2_corr_data['RA']), np.mean(evt2_corr_data['DEC']))
        assert evt2_off < 0.1


def test_fix_introduced_offset():
    """
    Introduce an arbitrary drift in a data set and show that the correction does not change the
    mean astrometry (from the drifted set to the corrected set)
    """

    ## Check that the tool can fix really drifted data and not change mean astrometry
    # Do this by constructing drifted data (just in sky x, y) and undoing with the tool
    obsid = 16659
    x = 4133.76
    y = 4078.74
    radius = 6
    bash("download_chandra_obsid {}".format(obsid), env=ciaoenv)

    asol = glob("{}/primary/*asol1.fits*".format(obsid))[0]
    evt2file = glob("{}/primary/*evt2.fits*".format(obsid))[0]
    evt1file = glob("{}/secondary/*evt1.fits*".format(obsid))[0]

    # reduce the evt1 file to just the source region in place
    bash("""dmcopy {evt1file}"[(x,y)=circle({x}, {y}, {radius})]" {obsid}_src_evt1.fits clobber+""".format(evt1file=evt1file, x=x, y=y, radius=radius, obsid=obsid),
         env=ciaoenv)

    # Reprocess it so it is just "like" other reprocessing we'll do later.
    bpix = glob("{}/primary/*_bpix1.fits*".format(obsid))[0]
    mtl = glob("{}/secondary/*_mtl1.fits*".format(obsid))[0]
    bash("""acis_process_events acaofffile= {asol} \
    infile="{obsid}_src_evt1.fits[cols -status]" \
    outfile= {obsid}_src_repro_evt1.fits \
    badpixfile= {bpix} \
    mtlfile= {mtl} \
    eventdef=")stdlev1" pix_adj=RANDOMIZE rand_pha=yes apply_tgain=yes apply_cti=yes clobber+""".format(
            asol=asol, obsid=obsid, bpix=bpix, mtl=mtl),
         env=ciaoenv)

    evt1_src_file = "{}_src_repro_evt1.fits".format(obsid)

    # Construct an aspect solution with artificially bad offsets
    # first, link to the version of the tool that just introduces bad offsets
    bash("rm correct_periscope_drift")
    bash("ln -s {}/introduce_bad_offset.py correct_periscope_drift".format(TESTDIR))
    cpd(asol, evt1_src_file, x, y, radius, 'bad_offsets')
    bad_asol = "bad_offsets_asol1.fits"

    # Apply those offsets to the source file to make evt1 with artificial "jet"
    bash("""acis_process_events acaofffile= {asol} \
    infile="{obsid}_src_repro_evt1.fits[cols -status]" \
    outfile= {obsid}_bad_offsets_evt1.fits \
    badpixfile= {bpix} \
    mtlfile= {mtl} \
    eventdef=")stdlev1" pix_adj=RANDOMIZE rand_pha=yes apply_tgain=yes apply_cti=yes clobber+""".format(
            asol=bad_asol, obsid=obsid, bpix=bpix, mtl=mtl),
         env=ciaoenv)

    # Rerun the correction tool on this bad file with the bad aspect solution file.
    # First link to the right/good version of the tool
    bash("rm correct_periscope_drift")
    bash("ln -s {}/../correct_periscope_drift.py correct_periscope_drift".format(TESTDIR))
    # Note a larger radius is needed to get all the events for this drifted data
    cpd(asol=bad_asol, evtfile="{}_bad_offsets_evt1.fits".format(obsid),
        x=x, y=y, radius=40, out='fixed_offsets')


    # Apply those offsets to the source file to correct artificial "jet"
    bash("""acis_process_events acaofffile= {asol} \
    infile="{obsid}_bad_offsets_evt1.fits[cols -status]" \
    outfile= {obsid}_fixed_evt1.fits \
    badpixfile= {bpix} \
    mtlfile= {mtl} \
    eventdef=")stdlev1" pix_adj=RANDOMIZE rand_pha=yes apply_tgain=yes apply_cti=yes clobber+""".format(
            asol="fixed_offsets_asol1.fits", obsid=obsid, bpix=bpix, mtl=mtl),
         env=ciaoenv)

    bash('dmcopy {evt2file}"[(x,y)=circle({x},{y},{radius})][cols ra,dec,time]" {obsid}_orig_evt2_radectime.fits clobber+'.format(
            evt2file=evt2file, x=x, y=y, radius=radius, obsid=obsid), env=ciaoenv)
    bash("dmcopy {}_fixed_evt1.fits'[cols ra,dec,time]' {}_fixed_evt1_radectime.fits clobber+".format(obsid, obsid), env=ciaoenv)
    bash("dmcopy {}_bad_offsets_evt1.fits'[cols ra,dec,time]' {}_bad_offsets_evt1_radectime.fits clobber+".format(obsid, obsid), env=ciaoenv)


    orig_aspect = Table.read(glob("{}/primary/*asol1*".format(obsid))[0])
    fixed_aspect = Table.read("fixed_offsets_asol1.fits".format(obsid))

    # Is the fixed solution a lot like the original?  This is specific to this obsid/case.
    # dy in mm not arcsec
    assert np.max(orig_aspect['dy'] - fixed_aspect['dy']) < .15

    uncorr_data = Table.read("{}_bad_offsets_evt1_radectime.fits".format(obsid))
    corr_data = Table.read("{}_fixed_evt1_radectime.fits".format(obsid))

    # means are unchanged
    corr = 3600 * sph_dist(np.mean(uncorr_data['RA']), np.mean(uncorr_data['DEC']),
                           np.mean(corr_data['RA']), np.mean(corr_data['DEC']))
    assert corr < 0.1


