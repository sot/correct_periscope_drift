# Correct Events for Aspect Periscope Drift

## Introduction and History

Thermal cycling on the spacecraft can result in an apparent temporal drift of the sky
position of an X-ray source during an observation.  This appears as a drift of up to about
0.5 arcsec in X-ray event sky X, Y coordinates over time.  Because of the thermal
variation time scales, this effect is usually most prominent in long observations (more
than about 50 ksec).

As of DS 8.4, a drift correction is applied to the aspect solution using the periscope
gradients telemetry.  However, temporal drifts have continued to increase with thermal
variation of the spacecraft.  Therefore the Aspect team suggests that to accomplish
science related to sub-arcsec source structure, users should follow this thread to correct
residual drift induced by the periscope.  This requires a relatively bright, on-axis source (within
a few arcmin off-axis angle) to perform a "self-calibration" of the aspect solution.

## Overview of determining a correction

To aid in determining the drift during an observation, we provide a script which takes as input:

 * Chandra evt1 or evt2 file
 * Chandra aspect solution
 * Coordinates and radius of a bright, point-like X-ray source
 * Degree of the desired fit polynomial

This script returns:

 * Plots of the fits of the drift in Aspect Camera Y and Z angle
 * New aspect solution file with offsets updated to include drift corrections

Users of the script may use the plots to evaluate the goodness of fit to their data and to
determine if applying the correction will have value.

## Applying a correction

This new aspect solution may be applied via the chandra_repro (via acis_process_events or
hrc_process_events) procedure to correct the sky positions in a new evt2 file.  See

http://cxc.harvard.edu/ciao/threads/createL2/

and

http://cxc.harvard.edu/ciao/ahelp/chandra_repro.html

## Note on coordinate systems

The correction is performed by directly updating the dy and dz values in the aspect
solution.  The dy/dz values are aligned with the Aspect Camera Y and Z axes, and
therefore, to fit the X-ray events, this tool transforms the event coordinates into
positions in the Aspect Camera frame, and then independently fits offsets in those Y and Z
axes.  Note that the positions in this frame are abbreviated 'yag' and 'zag' in the fit and
data plots output by the tool.

## Walkthrough/Example

### Get data

Fetch the data for an observation.

For gratings data, obviously only use events around the zeroth order position.

    download_chandra_obsid 17128

### Select a source

View and select a region to use for the autocorrection. Celldetect is one rough method to view and select a source.

    punlearn celldetect
    celldetect infile= 17128/primary/acisf17128N002_evt2.fits.gz outfile=mysrc.fits
    punlearn dmsort
    dmsort mysrc.fits'[SRCLIST]' key=-snr outfile=mysrc_snrsort.fits
    dmlist mysrc_snrsort.fits'[cols net_counts, x, y, snr, r][net_counts > 500]' data | head

    --------------------------------------------------------------------------------
    Data for Table Block SRCLIST
    --------------------------------------------------------------------------------
    ROW    NET_COUNTS           POS(X,Y)                         SNR   R[2]

     1             4730.250 (     4064.8883702839,     4071.8739789965)
     63.5722618103 [        1.9038963318         1.7280120850]
     2      1916.0833740234 (     5050.8106140536,     6133.5909474573)
     21.5115623474 [       98.9364013672        96.6939392090]
     3      2097.0068359375 (     4743.0030283711,     6043.9606311763)
     20.2611484528 [       78.0483093262        70.0926284790]


We are looking for a high count, high SNR source, relatively close to the optical axis (which for ACIS is
defined at x,y = 4096.5, 4096.5).  The first source will do.  We select a radius larger
than the larger reported value of R of the celldetect shape, and view in ds9.

    ds9 17128/primary/acisf17128N002_evt2.fits.gz \
    -bin about 4064.8883702839 4071.8739789965 \
    -scale log \
    -regions command "circle 4064.8883702839  4071.8739789965 # color=red" \
    -zoom 8

![ds9 screenshot of src](ds9_src.png)

The source looks point-like and is basically contained within the specified region.

### Setup the correction tool using the source for reference events for the correction

    punlearn correct_periscope_drift

Set the source for the tool.

    pset correct_periscope_drift x=4064.8883702839 y=4071.8739789965 radius=6

Set the other input and output files as desired

    pset correct_periscope_drift infile= 17128/primary/pcadf557756838N002_asol1.fits.gz
    pset correct_periscope_drift evtfile= 17128/primary/acisf17128N002_evt2.fits.gz
    pset correct_periscope_drift outfile=driftcorr_asol1.fits
    pset correct_periscope_drift corr_plot_root=demo_corr

### Run the tool

    correct_periscope_drift

The tool will write out an updated aspect solution to 'driftcorr_asol1.fits' as requested
by the outfile parameter and will save the plots of the fits into the working directory.

    Running: correct_periscope_drift
       version = 0.1
    with parameters:
      infile=17128/primary/pcadf557756838N002_asol1.fits.gz
      evtfile=17128/primary/acisf17128N002_evt2.fits.gz
      outfile=driftcorr_asol1.fits
      verbose=2
      and ASCDS_INSTALL is /soft/ciao-4.8
    ------------------------------------------------------------
    Fitting a line to the data to get reduced stat errors
    Fitting a polynomial of degree 2 to the data
    Fitting a line to the data to get reduced stat errors
    Fitting a polynomial of degree 2 to the data
    ------------------------------------------------------------
    Fit results
        Events show drift range of 0.36 arcsec in yag axis
        Max absolute correction of 0.21 arcsec for yag axis
        Events show drift range of 0.02 arcsec in zag axis
        Max absolute correction of 0.01 arcsec for zag axis
    ------------------------------------------------------------
    Writing out corrected aspect solution file to driftcorr_asol1.fits
        To review fit see correction plots in:
           demo_corr_fit_yag.png
           demo_corr_data_yag.png
           demo_corr_fit_zag.png
           demo_corr_data_zag.png

![y-angle fit and binned data](demo_corr_fit_yag.png)
![y-angle fit and raw data](demo_corr_data_yag.png)
![z-angle fit and binned data](demo_corr_fit_zag.png)
![z-angle fit and raw data](demo_corr_data_zag.png)

### Run chandra_repro with the new aspect solution.

    # move the original aspect solution out of the primary directory and rename
    mv 17128/primary/pcadf557756838N002_asol1.fits.gz 17128/pcadf557756838N002_asol1.fits.gz.ORIG
    cp driftcorr_asol1.fits 17128/primary/pcadf557756838_driftcorr_asol1.fits
    cd 17128
    chandra_repro
    cd ..

Note that the duration of the event processing (acis_process_events, hrc_process_events)
task in chandra_repro scales with the number of events and, for the long observations with
bright sources for which this tool is useful, may take more than 30 minutes to run on a
modern CPU.

### Verify the outputs

    pset correct_periscope_drift evtfile= 17128/repro/acisf17128_repro_evt2.fits
    pset correct_periscope_drift infile= 17128/primary/pcadf557756838_driftcorr_asol1.fits
    pset correct_periscope_drift outfile=already_fixed.fits
    pset correct_periscope_drift corr_plot_root="already_fixed"
    correct_periscope_drift

    Running: correct_periscope_drift
      version = 0.1
    with parameters:
      infile=17128/primary/pcadf557756838_driftcorr_asol1.fits
      evtfile=17128/repro/acisf17128_repro_evt2.fits
      outfile=already_fixed.fits
      verbose=2
      and ASCDS_INSTALL is /soft/ciao-4.8
    ------------------------------------------------------------
    Fitting a line to the data to get reduced stat errors
    Fitting a polynomial of degree 2 to the data
    Fitting a line to the data to get reduced stat errors
    Fitting a polynomial of degree 2 to the data
    ------------------------------------------------------------
    Fit results
        Events show drift range of 0.00 arcsec in yag axis
        Max absolute correction of 0.00 arcsec for yag axis
        Events show drift range of 0.00 arcsec in zag axis
        Max absolute correction of 0.00 arcsec for zag axis
    ------------------------------------------------------------
    Writing out corrected aspect solution file to already_fixed.fits
            To review fit see correction plots in:
                   already_fixed_fit_yag.png
                   already_fixed_data_yag.png
                   already_fixed_fit_zag.png
                   already_fixed_data_zag.png


![y-angle fit and binned data](already_fixed_fit_yag.png)
![y-angle fit and raw data](already_fixed_data_yag.png)
![z-angle fit and binned data](already_fixed_fit_zag.png)
![z-angle fit and raw data](already_fixed_data_zag.png)


## Details


The correction and fitting routine extracts the X-ray events from the provided source
region (a circle centered at the provided coordinates with the supplied radius) and converts
the RA, Dec of the X-ray events into approximately Aspect Camera/PCAD frame Y and Z
relative to the RA_PNT, DEC_PNT, ROLL_PNT supplied in the event list.  
We then fit a two independent curves to the mean-subtracted Y and Z angle data using a sherpa fit model.
A polynomial is used as the fit model; users may specify the degree of the desired
polynomial as an option to the tool. Advanced users may directly edit the Python fitting script to use a custom
model.

The Sherpa fit is then applied to the aspect solution 'dy', 'dz' columns and a new aspect
solution file with those updated columns is written out.

## API

The script uses the standard parameter interface with these allowed parameters:

 * infile - input aspect solution file
 * evtfile - event file
 * outfile - corrected/output aspect solution file
 * corr_plot_root - prefix for correction evaluation plots
 * x - src sky x
 * y - src sky y
 * radius - src circle radius in pixels
 * src_min_counts - minimum required src counts
 * corr_poly_degree - Degree of sherpa fit polynomial
 * clobber - Overwrite the output files if they exist?
 * verbose - Debug level (0=no debug information)
 * mode


## References

http://cxc.harvard.edu/mta/ASPECT/ECR_perifidcorr/ECR_pipe_perifidcorr.html
