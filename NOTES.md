Correct Events for Aspect Periscope Drift
-----------------------------------------


Introduction and History
========================

Temporal drift in X-ray source positions can occur due to changes in the periscope
alignment which induce drift in fid light positions as observed by the Aspect camera.

As of DS ????, a dynamic correction is applied to the aspect solution using the periscope
gradients.  This correction reduces ???? .  Based on observed data, there may be
uncorrected drift of up to ??? per ks.  These observed drifts are increasing ....
At this time, the Aspect team suggests that to accomplish
science requiring fine position or structure that users should use their fixed-position X-ray data to
autocorrect residual drift induced by the periscope.

Overview of determining a correction
====================================

To aid in determining the drift during an observation, we provide a script ???? which takes as input:

 * Chandra evt1 or evt2 file
 * Chandra aspect solution
 * the coordinates and radius of a bright, point-like X-ray source
 * the degree of the desired fit polynomial

This script returns:

 * Plots of the fits of the drift in Aspect Camera Y and Z angle
 * A new aspect solution file with offsets updated to include drift corrections

Users of the script may use the plots to evaluate the goodness of fit to their data and to
determine if applying the correction will have value.


Applying a correction
=====================

This new aspect solution may be applied via the _process_events procedure to correct the
sky positions in a new evt2 file.  See

????


Walkthrough/Example
===================

Fetch the data for an observation.

For gratings data, obviously only use events around the zeroth order position.

export obsid=16659
download_chandra_obsid $obsid

View and select a region to use for the autocorrection:

export evt2=$obsid/primary/acisf16659N001_evt2.fits.gz

Celldetect is one rough method to view and select a source

punlearn celldetect
celldetect infile=$evt2 outfile=mysrc.fits
punlearn dmsort
dmsort mysrc.fits'[SRCLIST]' key=-snr outfile=mysrc_snrsort.fits
bash-4.1$ dmlist mysrc_snrsort.fits'[cols net_counts, x, y, snr, r][net_counts > 500]' data |
head

--------------------------------------------------------------------------------
Data for Table Block SRCLIST
--------------------------------------------------------------------------------
 
ROW    NET_COUNTS           POS(X,Y)                                 SNR
R[2]
 
     1      6187.6665039062 (     4557.7041829672,     4545.8621727355)
     72.2447509766 [       10.6039648056         6.1736688614]
     2            3222.8750 (     4133.7619994320,     4078.7477989208)
     52.1157989502 [        1.8717634678         1.7822588682]
     3      1616.1666259766 (     3451.3598281418,     4158.0214822771)
     34.1368141174 [       11.2211189270         8.3461980820]


Looking for a high count, high SNR source, relatively close to the optical axis (which for ACIS is
defined at x,y = 4096.5, 4096.5).  The second source will do.  We select a radius larger
than the larger reported value of R of the celldetect shape, and view in ds9.

ds9 $evt2 \
-bin about 4133.7619994320 4078.7477989208 \
-scale log \
-regions command "circle 4133.7619994320 4078.7477989208 6 # color=red" \
-zoom 8

![ds9 screenshot of src](ds9_src.png)

The source looks point-like and is contained within the specified region.

punlearn correct_periscope_drift

Set the source for the tool.

pset correct_periscope_drift src_x=4133.7619994320 src_y=4078.7477989208 src_radius=6

Set the other input and output files as desired

pset correct_periscope_drift evt2file=$evt2
pset correct_periscope_drift input_asolfile=16659/primary/pcadf537654279N001_asol1.fits.gz

correct_periscope_drift

The tool will write out an updated aspect solution to 'corr.fits' by default and plots of
the fits and the event data.

![y-angle fit and binned data](corr_fit_yag.png)
![y-angle fit and raw data](corr_data_yag.png)
![z-angle fit and binned data](corr_fit_zag.png)
![z-angle fit and raw data](corr_data_zag.png)


Details
=======

The correction and fitting routine extracts the X-ray events from the provided source
region (a circle centered at the provided coordinates with the supplied radius) and converts
the RA, Dec of the X-ray events into approximately Aspect Camera/PCAD frame Y and Z
relative to the RA_NOM, DEC_NOM, ROLL_NOM supplied in the event list.  
We then fit a two independent curves to the mean-subtracted Y and Z angle data using a sherpa fit model.
A polynomial is used as the fit model; users may specify the degree of the desired
polynomial as an option to the tool. Advanced users may directly edit the Python fitting script to use a custom
model.

The Sherpa fit is then applied to the aspect solution 'dy', 'dz' columns and a new aspect
solution file with those updated columns is written out.

API
===

The script ???? uses the standard parameter interface with these allowed parameters:





References
==========

http://cxc.harvard.edu/mta/ASPECT/ECR_perifidcorr/ECR_pipe_perifidcorr.html
