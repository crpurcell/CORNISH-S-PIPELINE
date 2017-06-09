#-----------------------------------------------------------------------------#
# Description of the ATCA CABB continuum data reduction pipeline              #
#                                                                             #
# Cormac Purcell                                                              #
# 19-May-2017                                                                 #
#                                                                             #
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# 1_indxObsDay.py
#-----------------------------------------------------------------------------#

# USAGE: 
$> ./1_indxObsDay.py <YYYY-MM-DD.freq>
$> ./1_indxObsDay.py 2010-12-22.5500

# REQUIRES: 
One days worth of observing that has been loaded into MIRIAD format using the
ATLOD command, e.g.:
$> atlod in=2010-12-22/*.C1977 out=2010-12-22.uv options=birdie,noauto,rfiflag,xycor edge=5
and then split into files containing a single frequency window, e.g.:
$> uvsplit in=2010-12-22.uv options=nosource,clobber

# DESCRIPTION:
- Parse the table definition file "Imports/DBSchema.sql".
- parse the config file, including any gross flag statements.
- Backup the original tables (flags, header, history, vartable).
- Apply any global UVFLAG commands with the tag ".onindex".
- Connect to existing database or create a new one based on the schema file.
- Run UVINDEX and parse to find observation parameters.
- Fill the tables 'observations', 'field_coords' and 'pointings'.
- Create the entries in the 'mosaic_obs' table and insert into the DB.
  Note: this table is necessary as the pointing names are not unique.

# Script creates the following files:
DATA/DB_<freq>.sqlite (sqlite3 database, tables described below).
DATA/uvdata/<YYYY-MM-DD>.<freq>.uvindex (index file for the day).

# Database tables:
# Column name       Example entry

# observations:
timeStamp (KEY)     10DEC22:19:34:34.9
dayName             2010-12-22
pntName             0823-500
calCode             c 
nAnts               6
fieldName           082526.87-501038.49 (named for X and Y position)
RA_deg              126.36195833333333
Dec_deg             -50.17735833333333
splitFlag           0
useFlag             0

# field_coords:
fieldName (KEY)     082526.87-501038.49
RA_deg              126.36195833333333
Dec_deg             -50.17735833333333
calCode             c 

# mosaic_obs:
dayPntID            2010-12-22.0823-500
dayName             2010-12-22
pntName             0823-500
calCode             c
RA_deg              126.36195833333333
Dec_deg             -50.17735833333333
splitFlag           0
useFlag             0

# pointings:
pntName             0823-500
calCode             c 
RA_deg              126.36195833333333
Dec_deg             -50.17735833333333
fieldName           082526.87-501038.49


#-----------------------------------------------------------------------------#
# 2_calObsDay.py
#-----------------------------------------------------------------------------#

# USAGE:
$> ./2_calObsDay.py [-p] <YYYY-MM-DD.config>
$> ./2_calObsDay.py -p 2010-12-22.5500

# REQUIRES:
The outputs of 1_indxObsDay.py and a config file.

# DESCRIPTION:
1) Read and parse the config file. Flagging commands are given special
treatment and are run at different stages, depending on the extension
of their 'keyword' in the config file.
2) Backup the tables in the split file (flags, header, history, vartable).
3) Apply flagging commands with the extension "precal".
4) Create some diagnostic plots in ./PLOTS: 
   * uvCov.ps = uv-coverage plot
   * el.ps = elevation vs time
   * ampCalRaw.ps = uncalibrated XX amplitudes vs time
   * TsysX.ps and TsysY.ps = antenna Tsys vs time for X and Y feeds
   * Smon.ps = rms seeing monitor measurement versus time
   * XYphase.ps = XY-phase versus time for each antenna
   * axisRMS.ps = RMS tracking error versus time for each baseline
   * axisMax.ps = maximum tracking error versus time for each baseline
5) Split off a copy of the bandpass calibrator (bpcal), unflagged data
only.
6) Run flagging commands with extension "bpcalN" on the bpcal.
For PGFLAG, flag in order: XX, YY, XY, YX & apply to all.
7) Calibrate the bandpass (and gains) of the bpcal. Plot the solutions.
   MFCAL in=bpcal options=interpolate
   GPCAL options=xyvary nfbin=<nfbinBP>
   Plots:
         * SolnBPamp.ps = amplitude solutions
         * SolnBPpha.ps = phase solutions
         * BPcalSpec.ps = average spectrum of bpcal after gains applied
         * BPcalReIm.ps = Re vs Im after gains have been applied
8) Split off copies of the phase calibrators (pcal<n>), copy over bandpass.
9) Run flagging commands with extension "phcalN" on each phase calibrator.
For PGFLAG, flag in order: XX, YY, XY, YX & apply to all.
10) Loop through each phase-cal in turn and run gpcal:
   GPCAL vis=phcal options=xyvary,qusolve nfbin=<nfbin>
11) Merge solutions into first pcal: gpcopy mode=merge options=nopass
12) Plot the solutions:
    * SolnGainAmpXX/YY/XY = XX/YY/XY amplitude solutions versus time
    * SolnGainPhaXX/YY/XY = XX/YY/XY phase solutions versus time
    * pcal<n>.GAIN.<freq>_Amp.ps = amplitude vs time of pcal after calibration
    * pcal<n>.GAIN.<freq>_Pha.ps = phase vs time of pcal after calibration
    * PcalSpec.ps = average spectrum of each calibrator after gains applied
13) Split off a copy of the flux calibrator (fcal), copy over bandpass.
14) Run flagging commands with extension "fluxcalN" on each phase calibrator.
For PGFLAG, flag in order: XX, YY, XY, YX & apply to all.
15) Calibrate and bootstrap the flux calibrator:
   GPCAL options=xyvary nfbin=<nfbin>
   GPBOOT vis=pcal cal=fcal
   MFBOOT vis=fcal,pcal select=source(fluxCalName)
16) Optional: apply calibration to a secondary "test" calibrator and plot the
spectrum.
17) Run 2 passes of steps 5-16. The flagging commands applied during pass 2
operate on data calibrated in the 1st pass and may be defined differently.
18) Copy all calibration tables back to the multi-source file.
19) Loop through each integration/pointing and split into a uv-file
file.
   * Query the start & end timestamps from the DB, create select statement.
   * Split out the pointing using UVCAT and the selection statement
   * Set a flag in the DB indicating a successful split (observations table).
   * Run flagging commands on each integration/pointing with labels "onsplit".


#-----------------------------------------------------------------------------#
# 3_imageField.py
#-----------------------------------------------------------------------------#
# NOTE: This script needs to be updated before using in anger as part of the
pipleine process:
* Convert the method for uv-data concatenation to the link method.
* Update the database with parameters of the images.

# USAGE:
$> ./3_imageField.py [-p] [-f <freqExt>] [-d <dayName>] <pntName>    #

# e.g., to image # 1729-37:
$> ./3_imageField.py -p -f 5500 173315.19-372232.40 # 1729-37

# REQUIRES:
Data from individual integrations split into the directory "DATA/uvdata_split"
by the 2_calObsDay.py script. Book-keeping data is stored in the database for
the relevant frequency configuration, e.g. "DATA/DB_5500.sqlite"

# DESCRIPTION:
Image a single field from the CORNISH-South dataset. Supplying the day name
(using the -d argument) limits the data to a particular day only.
1) Read a parse the default imaging configuration file (set at the top of the
script, e.g., "imaging5500.config" for frequency extension 5500).
2) Query the database for the uv-data in the field.
3) Concatinate the integrations into a single time-order dataset.
4) Check for a custom imaging configuration file e.g.,
"./customImgConfigs/173315.19-372232.40_img5500.config"
if found, replace default entries with custom values (can be a small subset).
5) Run flagging commands with extension "preimg" on the field.
For PGFLAG, flag in order: XX, YY, XY, YX & apply to all.
6) Calculate the default imaging parameters from the dataset.
* Use UVLIST with options="spectral" to get the data parameters.
* Determine the primary beam size from 1.22 * C/(max(freq)*22m)
* Determine the resolution from 1.22 * C / (max(freq) * max(baseline))
* Cell size: default to 5 pixels across the sybthesised beam.
* Image size: set from fov_FWHM in config file
7) For each output channel plot the spectrum of that chunk before averaging.
8) Image the field in Stokes V and measure the RMS.
* Call INVERT using the same parameters as Stokes I imaging
* Call IMSTAT to measure the RMS noise and parse the log to read.
9) Image the field in Stokes I using INVERT
* invert cell=<cellsize_asec> robust=<robust> options=mfs,sdb,mosaic stokes=i line=<channel_selection>
10) Deconvolve the dirty image using MFCLEAN down to RMSv*mult.
* mfclean gain=<gain> niter=<niter> region=perc(33) minpatch=<minpatch> speed=<speed>
11) Output a clean model and residual using RESTOR.


#-----------------------------------------------------------------------------#
# 4_genEquTileGrid.py
#-----------------------------------------------------------------------------#

# USAGE:
$> ./4_genEquTileGrid.py

# REQUIRES:
A DATA directory containing .sqlite database files that have been created by
the script 1_indxObsDay.py

# DESCRIPTION:
 Variables set at the begining of the script define the Galactic
 border of the CORNISH-South survey and the properties of the
 tiles. This script sets down tiles along constant lines of
 declination, overlapping by ~60 arcsec. Once the tile centres have
 been calculated a 'tile_coords' table is written to all '.sqlite'
 database files found in the '../DATA/' directory.

# Note 1: The script only needs to be run once to fill the "tile_coords" table
  but can be run multiple times without harm (e.g., to visualise the tile
  positions).

# Note 2: The script does not try to correct for the curvature of the
  lines of constant Dec. towards the Equatorial poles, but does
  correct for the cos(Dec.) factor. I think we get away with this 
  simplification in CORNISH-South as we do not approach the pole
  closely and the tiles are relatively small.

# Note 3: You should set the Galactic borders of the region large
  enough so that the pattern of tiles overlap the imaged data out to
  the edge of the primary beam FWHM at the lowest frequency.


#-----------------------------------------------------------------------------#
# 5_imageTileLin.py
#-----------------------------------------------------------------------------#

# USAGE:
$> ./5_imageTileLin.py [-p] [-f <freqExt>] [-d <dayName>] <tileID>    #

# e.g., to image Tile 1481
$> ./5_imageTileLin.py -f 5500 1481

# REQUIRES:
Data from individual integrations split into the directory "DATA/uvdata_split"
by the 2_calObsDay.py script. Book-keeping data is stored in the database for
the relevant frequency configuration, e.g. "DATA/DB_5500.sqlite"

# DESCRIPTION:
Image one or more tiles from the CORNISH-South dataset using MIRIAD. Images all
possible fields that overlap the tile footprint, combines them onto a tile and
crops the tile to a 4000 x 4000 pixel image.

1) Parse the imaging configuration file (imaging<ext>.config) and open the
database for the current frequency extension.
2) Query the database for parameters of the requested tile (coordinate, extent,
pixel-scale).
3) Query the database for the fields overlapping the tile and their parameters
(coordinate, extent, pixel-scale)
4) Loop through each pointing querying the uv-data recorded for each. Construct
a table of uv-data files for each field with columns indicating whether the
data has been split out to disk or flagged not to use.
5) Run UVLIST to get the spectral & uv parameters of the data. Calculate the
default image parameters from the output and parameters in the config file
(field of view, frequency range, pimary beam FWHM range, resolution range)
6) Create a KVIS primary beam annotation file showing the pointing pattern the
boundary of the tile.
7) Set the pixel size as 3x the high-frequency resolution and the image size
as the low-frequency field-of-view x 'fov_FWHM' multiplier in the configuration
file.
8) Loop through the pointings:
  - Create a temp working directory containing symbolic links to the uv-data.
  - Format selection string to accomplish MFS into one channel
  - Plot the spectrum of the data
  - Create a small (0.5xFWHM) Stokes V image and measure the RMS noise.
  - Create a full-size Stokes I image.
  - Run MFCLEAN down to a cutoff of 'rms_mult' x RMS.
  - Restore the CC map and create a residual.
  - Delete the temp working directory.
9) Run LINMOS to stitch the field images into a square tile.
