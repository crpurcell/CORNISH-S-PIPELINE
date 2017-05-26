#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     5_imageTileLin.py                                                 #
#                                                                             #
# USAGE:    ./5_imageTileLin.py [-p] [-f <freqExt>] [-d <dayName>] <tileID>   #
#           ./5_imageTileLin.py -f 5500 1481                                  #
#                                                                             #
# PURPOSE:  Image a tile from a mosaicked ATCA dataset by inverting and       #
#           deconvolving each field separately and then mosaicing.            #
#                                                                             #
# NOTE:     If <freqExt> is not supplied then it defaults to '5500'.          #
#                                                                             #
# MODIFIED: 24-May-2017 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Hardcoded paths
dataRootDir = "../DATA"

# Default frequency extension
defIFext = '5500'

# Default imaging configuration file for each frequency
defImgConfigs = {"5500": "imaging5500.config",
                 "9000": "imaging9000.config"}

# Constants
C = 2.998e8

# Radius of the fields to use when determining overlap
radFld_deg = 12.7/60.0 # @4.5 GHz

#-----------------------------------------------------------------------------#

import os
import sys
import shutil
import copy
import re
import random
import string
import optparse
import getopt
import time
import math as m
import numpy as np
from numpy.lib import recfunctions as rec
import sqlite3

from Imports.util_ATCA_PIPE import *
from Imports.util_tables import *
register_sqlite3_numpy_dtypes()   # Map numpy data-types to sqlite3 data-types
print "Setting up the MIRpy environment ...",
sys.stdout.flush()
from Imports.mirpy import miriad as mir
print "done."


#-----------------------------------------------------------------------------#
def main():
    
    # Misc setup
    IFext = defIFext
    doPause = False
    redoImg = False
    dayName = None
    startTime = time.time()
    
    # Check for required files and directories
    startDir = os.getcwd()
    fail_not_exists(dataRootDir, 'directory')
    uvDataSplitDir = dataRootDir + "/uvdata_split"
    fail_not_exists(uvDataSplitDir, 'directory')
    imageFieldDir = dataRootDir + "/images_fields"
    imageTileDir =  dataRootDir + "/images_tiles"
    reqDirLst = ['LOGS', 'PLOTS', imageFieldDir, imageTileDir]
    for reqDir in reqDirLst:
        if not os.path.exists(reqDir):
            os.mkdir(reqDir)
            
    # Get the options from the command line
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'pf:d:r')
    except getopt.GetoptError:
        print_usage()
        sys.exit()
    if len(args)< 1:
        print "\nPlease specify a tile ID on the command line:"
        print_usage()
        sys.exit()
    tileID = args[0]
    for o, a in opts:
        if o == "-f":
            IFext = a
        if o == "-d":
            dayName = a
        if o == "-p":
            doPause = True
        if o == "-r":
            redoImg = True
    log_wr(None, "\n> Using the frequency extension '%s'." % IFext)
    
    # Open a log file
    LF = open('LOGS/Tile_' + tileID + '.' + IFext + '_img.log', 'w')

    # Read and parse the default imaging configuration file
    configFile = defImgConfigs[IFext]
    fail_not_exists(configFile, 'file')
    try:
        pDict, flgTab = config_read_flag(configFile)
        log_wr(LF, "\n> Successfully parsed the default imaging config file.")
    except Exception:
        log_fail(LF, "\n> Err: Failed to parse the imaging config file.")

        
    # Connect to the database file for the current frequency extension
    dbFile = dataRootDir + '/DB_' + str(IFext) + '.sqlite'
    if os.path.exists(dbFile):
        try:
            log_wr(LF, "\n> Connecting to DB file '%s' ..." % dbFile)
            conn = sqlite3.connect(dbFile)
            conn.create_function("REGEXP", 2, regexp)
            cursor = conn.cursor()
        except Exception:
            log_fail(LF, "\n> Err: Failed to connect to DB '%s'" % dbFile)

    # Query the database for the tile parameters
    log_wr(LF, "\n> Querying tile parameters ...")
    tlTab = query_tile_parameters(cursor, tileID)
    if len(tlTab)==0:
        log_fail(LF, "\n> Err: No parameters found for that tile.")

    # Query the database for the pointings in the tile
    # Return {pntName: [RA_deg, Dec_deg], ...}
    log_wr(LF, "\n> Querying the pointings in tile %s ..." % tileID)
    pntParmDict = query_fields_in_tile(cursor, tileID, radFld_deg)
    if len(tlTab)==0:
        log_fail(LF, "\n> Err: No fields found underneath tile footprint.")
    else:
        log_wr(LF, "\n> Found %d fields underneath the tile footprint" % 
               len(pntParmDict))
        
    # Loop through the pointings in the tile, querying uv-data in each
    # Return RecArr:
    # dayPntID, dayName, pntName, calcode, RA_deg, Dec_deg, splitFlag, useFlag
    # Want to construct a dictionary of existing uv-data per pntName and
    # a flat list of uv-data
    obsTabLst = []
    for pntName in pntParmDict.keys():
        log_wr(LF, "\n> Querying uv-data for '%s' ..." % pntName)

        # Query observations in each pointing & return a recArray
        # timeStamp, dayName, pntName, calCode, nAnts, fieldName, RA_deg, 
        # Dec_deg, splitFlag, useFlag 
        obsTab = query_obs_in_field(cursor, pntName, splitFlag=False,
                                       useFlag=False, dayName=dayName)

        # Check that the data exists. If it does append the row to a list
        for row in obsTab:
            uvDataPnt = row['timeStamp'] +  '_' + IFext + '.uv'
            uvDataPntPath =  uvDataSplitDir + '/' + uvDataPnt
            if os.path.exists(uvDataPntPath):
                exists = 'Yes'
                obsTabLst.append(row)
            else:
                exists = 'No'
            log_wr(LF, "     %s exists? %s." %  (uvDataPntPath, exists))

    # Stack the list of recArrays into a new master observation table
    obsTab =  rec.stack_arrays(obsTabLst, autoconvert=True)
    obsTab.sort(order=['dayName','pntName'])
    log_wr(LF, ">> Will image the following uv-data:")
    for row in obsTab:
        log_wr(LF, "%s" % row)

    # Calculate the default image parameters
    uvData0 = uvDataSplitDir + '/' + obsTab['timeStamp'][0] +  '_' + \
              IFext + '.uv'
    outLogTmp = str(tileID) + '_' + IFext + '.uvlist.dat'
    mir.uvlist(vis=uvData0, options='array', log=outLogTmp)
    baselineDict = parse_uvlist_array(outLogTmp)
    os.remove(outLogTmp)
    outLogTmp = str(tileID) + '_' + IFext + '.uvlist.dat'
    mir.uvlist(vis=uvData0, options='spectral', log=outLogTmp)
    specDict = parse_uvlist_spectral(outLogTmp)
    pDict['nChan'] = specDict['nChan']
    os.remove(outLogTmp)
    baselineMax_m = max(baselineDict.keys())
    freqArr_Hz = make_lin_axis(naxis=specDict['nChan'],
                                crpix=specDict['startChan'],
                                crval=specDict['startFreq_GHz'],
                                cdelt=specDict['freqInc_GHz']) * 1e9
    freqMed_Hz = np.median(freqArr_Hz)
    freqMax_Hz = np.max(freqArr_Hz)
    freqMin_Hz = np.min(freqArr_Hz)
    FWHMmin_deg = m.degrees( 1.22 * C / (freqMax_Hz * 22.0) )
    FWHMmin_asec = FWHMmin_deg * 3600.0
    FWHMmax_deg = m.degrees( 1.22 * C / (freqMin_Hz * 22.0) )
    FWHMmax_asec = FWHMmax_deg * 3600.0
    resolutionHigh_deg = m.degrees( 1.22 * C / (freqMax_Hz * baselineMax_m) )
    resolutionHigh_asec = resolutionHigh_deg * 3600.0
    resolutionLow_deg = m.degrees( 1.22 * C / (freqMin_Hz * baselineMax_m) )
    resolutionLow_asec = resolutionLow_deg * 3600.0
    log_wr(LF, "\n> Field size:")
    log_wr(LF, "  %.2f - %.2f arcmin (%.2f - %.2f GHz)" % (FWHMmin_asec/60,
                                                           FWHMmax_asec/60,
                                                           freqMax_Hz/1e9,
                                                           freqMin_Hz/1e9))
    log_wr(LF, "\n> Resolution:")
    log_wr(LF, "  %.2f - %.2f arcsec (%.2f - %.2f GHz)" % (resolutionHigh_asec,
                                                          resolutionLow_asec,
                                                          freqMax_Hz/1e9,
                                                          freqMin_Hz/1e9))
    
    # Calculate the positions of the tile corners
    dDec_deg = abs(tlTab['nPixY'] * tlTab['pixscaleY_asec'] / 3600.0)
    decHigh_deg = tlTab['Dec_deg'] + dDec_deg / 2.0
    decLow_deg = tlTab['Dec_deg'] - dDec_deg / 2.0
    dRA_deg = tlTab['nPixX'] * tlTab['pixScaleX_asec'] /3600.0/2.0
    dRAtop_deg = dRA_deg/(m.cos(m.radians(decHigh_deg)))
    dRAbot_deg = dRA_deg/(m.cos(m.radians(decLow_deg)))
    tlRAtopL_deg = tlTab['RA_deg'] + dRAtop_deg # RA, top-left
    tlRAtopR_deg = tlTab['RA_deg'] - dRAtop_deg # RA, top-right
    tlRAbotL_deg = tlTab['RA_deg'] + dRAbot_deg # RA, bottom-left
    tlRAbotR_deg = tlTab['RA_deg'] - dRAbot_deg # RA, bottom-right
    
    # Create a primary beam annotation file    
    ANN = open(imageTileDir + '/Tile_' + str(tileID) + '_' + str(IFext) + \
               '.ann', 'w')
    for row in obsTab:        
        ANN.write("COLOUR MAGENTA\n")
        ANN.write("CIRCLE W %f %f %f\n" % (row['RA_deg'],
                                           row['Dec_deg'],
                                           FWHMmin_deg/2.0 ))
        ANN.write("COLOUR GREEN\n")
        ANN.write("CIRCLE W %f %f %f\n" % (row['RA_deg'],
                                           row['Dec_deg'],
                                           FWHMmax_deg/2.0 ))
        ANN.write("COLOUR WHITE\n")
        ANN.write("CLINES W %f %f %f %f %f %f %f %f %f %f\n" % \
                  (tlRAtopL_deg, decHigh_deg,
                   tlRAtopR_deg, decHigh_deg,
                   tlRAbotR_deg, decLow_deg,
                   tlRAbotL_deg, decLow_deg,
                   tlRAtopL_deg, decHigh_deg))
    ANN.close()

    # Set the image resolution parameters: 3 pixels across a beam,
    # Image size specified in the driving parameter file in terms of the
    # primary beam at the lowest frequency
    if not 'cellsize_asec' in pDict:
        pDict['cellsize_asec'] = resolutionHigh_asec / 3.0
    FWHMmax_pix = int(m.ceil(FWHMmax_asec/float(pDict['cellsize_asec'])))
    pDict['imsize_pix'] = int(m.ceil(FWHMmax_pix * float(pDict['fov_FWHM'])))
    pDict['imsizeV_pix'] = int(m.ceil(FWHMmax_pix * 0.5 ))
    log_wr(LF, "\n> Imaging parameters:")
    log_wr(LF, "  Pixel size: %.2f arcsec"  % float(pDict['cellsize_asec']))
    log_wr(LF, "  Image size: %d pixels (%.2f arcmin)"  % \
           (pDict['imsize_pix'],
            (pDict['imsize_pix']*float(pDict['cellsize_asec'])/60.0)))

    # Loop through the mosaiced pointings and image
    chImgDict = {}
    imgTabLst = []
    
    pntNameLst = pntParmDict.keys()
    print pntNameLst
    pntNameLst.sort()
    log_wr(LF, "\n> Will image %d pointings ..." % len(pntNameLst))
    for pntName in pntNameLst:

        # Format the input data string
        log_wr(LF, "\n> Processing mosaicked pointing '%s'" % pntName)        
        uvSplitDataLst = obsTab[obsTab['fieldName']==pntName]['timeStamp']
        uvSplitDataLst = uvSplitDataLst.tolist()
        uvSplitDataLst = [e + '_' + str(IFext) + '.uv' for e in uvSplitDataLst]

        # Create a temporary working directory and create symbolic links to
        # the uv-data. This is necessary to circumvent the limited input
        # string length for file names in MIRIAD. We will clean up aftewards.
        tmpDir = ''.join(random.sample(string.letters+string.digits, 8))
        log_wr(LF, "\n> Creating working dir '%s' and linking in data ..." %
               tmpDir)
        if os.path.exists(uvDataSplitDir + '/' + tmpDir):
            shutil.rmtree(uvDataSplitDir + '/' + tmpDir, ignore_errors=True)
        os.mkdir(uvDataSplitDir + '/' + tmpDir)
        os.chdir(uvDataSplitDir + '/' + tmpDir)
        for uvSplitData in uvSplitDataLst:
            os.symlink('../' + uvSplitData, uvSplitData)
        os.chdir(startDir)
        uvData = uvDataSplitDir + '/' + tmpDir + '/*' 

        # Image a pointing in Stokes I
        imgTab = image_pointing(uvData, pntName, imageFieldDir, IFext, pDict,
                                stokes="I", redo=redoImg, LF=LF)
        imgTabLst.append(imgTab)
        
        # Clean up
        if os.path.exists(uvDataSplitDir + '/' + tmpDir):
            shutil.rmtree(uvDataSplitDir + '/' + tmpDir, ignore_errors=True)

        # Insert the field paramaters into the database
        insert_arr_db(cursor, imgTab, 'field_images')
        conn.commit()
            
    # Merge the imaging results into a master image table.
    # Filter for unsuccessful imaging attempts
    imgTab =  rec.stack_arrays(imgTabLst, autoconvert=True)
    imgTab = imgTab[imgTab['success']==1]

    # Linmos the pointing images files into a mosaic file for each channel
    log_wr(LF, "\n> Creating the mosaicked maps ...")

    # Linmos the pointing images files into a mosaic file and crop
    log_wr(LF, "\n> Creating the mosaicked maps ...")
    chanLst = np.unique(imgTab['chan']).tolist()
    for chan in chanLst:
        log_wr(LF, ">> Stitching maps for channel %d, Stokes I ..." % chan)
        imgTabFilt = imgTab[imgTab['chan']==chan]
        imgDataLst = imgTabFilt['imgClnData'].tolist()
        log_wr(LF, "   %s" % imgDataLst)
        tmpDir = ''.join(random.sample(string.letters+string.digits, 8))
        if os.path.exists(imageFieldDir + '/' + tmpDir):
            shutil.rmtree(imageFieldDir + '/' + tmpDir, ignore_errors=True)
        os.mkdir(imageFieldDir + '/' + tmpDir)
        os.chdir(imageFieldDir + '/' + tmpDir)
        for imgData in imgDataLst:
            os.symlink('../' + imgData, imgData[:28])
        os.chdir(startDir)
        imgIMos = imageTileDir + '/Mos' + str(tileID) + '_' + \
                  'CH%d' % int(chan) + '_' + str(IFext) + '_' + "I" + \
                  '.xy'        
        if os.path.exists(imgIMos):
            shutil.rmtree(imgIMos, ignore_errors=True)
        print mir.linmos(_in=imageFieldDir + '/' + tmpDir + '/*',
                         out=imgIMos)
        
        # Find the reference pixel of the mosaic file
        outLogTmp = imgIMos + ".PRTHD.log"
        if os.path.exists(outLogTmp):
            os.remove(outLogTmp)
        mir.prthd(_in=imgIMos, log=outLogTmp)
        cDict = parse_prthd(outLogTmp)
        os.remove(outLogTmp)
        
        # Calculate the crop boundaries
        dx = (tlTab['RA_deg'] - cDict["x"]*15.0)/cDict["dx"]
        dy = (tlTab['Dec_deg'] - cDict["y"])/cDict["dy"]
        regionStr = "relpixel,boxes(-2000,-2000,1999,1999)"
        regionStr = "relpixel,boxes(%d,%d,%d,%d)" % (-2000+dx, -2000+dy,
                                                     1999+dx, 1999+dy)

        
        ANN = open(imageTileDir + '/Tile_' + str(tileID) + '_' + str(IFext) + \
                   'REFPIX.ann', 'w')
        ANN.write("COLOUR GREEN\n")
        ANN.write("CROSS W %f %f %f %f\n" % (cDict["x"]*15.0,
                                             cDict["y"],
                                             FWHMmax_deg/20.0,
                                             FWHMmax_deg/20.0 ))
        ANN.write("COLOUR WHITE\n")
        ANN.write("CROSS W %f %f %f %f\n" % (tlTab['RA_deg'],
                                             tlTab['Dec_deg'],
                                             FWHMmax_deg/20.0,
                                             FWHMmax_deg/20.0 ))
        ANN.close()
        # REF:    9493 9150
        # CENTRE: 5071 5161
        # DIFF:   4422 3989
        
        # Crop the tile to a square boundary
#        imgITile = imageTileDir + '/Tile' + str(tileID) + '_' + \
#                   'CH%d' % int(chan) + '_' + str(IFext) + '_' + "I" + \
#                   '.xy'        
#        if os.path.exists(imgITile):
#            shutil.rmtree(imgITile, ignore_errors=True)
#        print mir.imsub(_in=imgIMos, out=imgITile, region=regionStr)

        # Clean up
        if os.path.exists(imageFieldDir + '/' + tmpDir):
            shutil.rmtree(imageFieldDir + '/' + tmpDir, ignore_errors=True)

            
            
    #region=relcenter,boxes(-2000,-2000,1999,1999)    
    endTime = time.time()
    print "Duration = %.2f min" % ((endTime-startTime)/60.0)


#-----------------------------------------------------------------------------#
def image_pointing(uvData, pntName, imageFieldDir, IFext, pDict, stokes="I",
                   redo=True, LF=None):

    # Create a recarray to hold the results
    dtype = [('fieldName', 'a50'),
             ('chan', 'i8'),
             ('lineStr', 'a20'),
             ('stokes', 'a2'),
             ('rms', 'f4'),
             ('rmsV', 'f4'),
             ('cutoff', 'f4'),
             ('imsize_pix','i8'),
             ('cellsize_asec','f4'),
             ('robust','f4'),
             ('niter','i8'),
             ('minpatch','i8'),
             ('speed','f4'),
             ('imgClnData', 'a50'),
             ('beamMajInt_asec','f4'),
             ('beamMinInt_asec','f4'),
             ('beamPAint_asec','f4'),
             ('beamMaj_asec','f4'),
             ('beamMin_asec','f4'),
             ('beamPA_asec','f4'),
             ('success','i8')]
    outTab = np.zeros((int(pDict['nchanImg'])), dtype=dtype)
    outTab['fieldName'] = pntName
    outTab['stokes'] = stokes
    outTab['imsize_pix'] = int(pDict['imsize_pix'])
    outTab['cellsize_asec'] = float(pDict['cellsize_asec'])
    outTab['robust'] = float(pDict['robust'])
    outTab['niter'] = int(pDict['niter'])
    outTab['minpatch'] = int(pDict['minpatch'])
    outTab['speed'] = float(pDict['speed'])
    
    # Parse the restoring beam size or set to blank
    beamFWHM = ""
    beamPA = ""
    if "beam" in pDict:
        try:
            beamFWHM = [float(x) for x in pDict["beam"][:2]]
            beamPA = float(pDict["beam"][2])
            log_wr(LF, ">> Forcing restoring beam to %s " % pDict["beam"])
        except Exception:
            log_wr(LF, ">> Missformed beam value imaging config file ...")
            log_wr(LF, ">> Setting to default ...")
    
    # Loop through the chosen number of image bins and use the line
    # selection parameter to choose and average channels to image
    specPltLst = []
    nchanBin = int(m.floor(pDict['nChan']/int(pDict['nchanImg'])))
    for i in range(int(pDict['nchanImg'])):
        
        # Line selection
        lineStr = "channel,%d,%d,%d" % (nchanBin,(1+i*nchanBin),1)
        outTab[i]['lineStr'] = lineStr
        outTab[i]['chan'] = int(i+1)
        
        # Plot the spectrum of each chunk
        specPlt = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                  '_' + str(IFext) + '_spec.ps'
        if os.path.exists(specPlt):
            os.remove(specPlt)
        log_wr(LF, "> PLOTTING: '%s'" % specPlt)
        try:
            print mir.uvspec(vis=uvData, line=lineStr, stokes='xx,yy',
                             interval=1e3, options='nobase,avall', nxy='1,1',
                             device=specPlt + '/cps')
        except Exception:
            log_wr(LF, "Plotting failed!")            
        if os.path.exists(specPlt):
            specPltLst.append(specPlt)

        # Call invert to create a small Stokes V dirty image
        log_wr(LF, "\n> Imaging a Stokes V patch.")
        imgVData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                   '_' + str(IFext) + '_V.xy'
        if redo and os.path.exists(imgVData):
            log_wr(LF, ">> Deleting previous version first ...")
            shutil.rmtree(imgVData, ignore_errors=True)
        if os.path.exists(imgVData):
            log_wr(LF, ">> Skipping: previous version exists ...")
        else:
            try:
                logInvert = mir.invert(vis=uvData,
                                       map=imgVData,
                                       imsize=pDict['imsizeV_pix'],
                                       cell=pDict['cellsize_asec'],
                                       robust=pDict['robust'],
                                       options='mfs,sdb,mosaic',
                                       stokes='v',
                                       line=lineStr)
            except Exception:
                log_wr(LF, ">> Err: failed to image '%s_CH%s' in Stokes V" %
                       (pntName, (i+1)))
                continue            
        print

        # Measure the RMS noise in the centre of the V image
        log_wr(LF, "\n> Measuring the RMS noise in Stokes V")
        outLogTmp = pntName + '_' + 'CH%d' % (i+1) + '_' + str(IFext) + \
                    '.imstat.dat'
        print mir.imstat(_in=imgVData, log=outLogTmp)
        statDict = parse_imstat(outLogTmp)
        os.remove(outLogTmp)
        cutoff = float(pDict['rms_mult']) * float(statDict['rms'])
        log_wr(LF, "\n> RMS=%f, MULT=%.2f, CUTOFF=%f" % \
               (float(statDict['rms']), float(pDict['rms_mult']), cutoff))
        outTab[i]['rmsV'] = float(statDict['rms'])
        outTab[i]['cutoff'] = cutoff
        
        # Call invert to create a dirty map
        log_wr(LF, "\n> Creating a dirty image ...")
        imgData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                  '_' + str(IFext) + '_' + stokes + '.xy'
        imgBeam = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                  '_' + str(IFext) + '_' + stokes + '.beam'
        if redo and (os.path.exists(imgData) or os.path.exists(imgBeam)):
            log_wr(LF, ">> Deleting previous version first ...")
            shutil.rmtree(imgData, ignore_errors=True)
            shutil.rmtree(imgBeam, ignore_errors=True)
        if os.path.exists(imgData):
            log_wr(LF, ">> Skipping: previous version exists ...")
        else:
            try:
                print mir.invert(vis=uvData,
                                 map=imgData,
                                 beam=imgBeam,
                                 imsize=pDict['imsize_pix'],
                                 cell=pDict['cellsize_asec'],
                                 robust=pDict['robust'],
                                 options='mfs,sdb,mosaic',
                                 stokes=stokes,
                                 line=lineStr)
            except Exception:
                log_wr(LF, ">> Err: failed to image '%s_CH%s' in Stokes %s" %
                       (pntName, (i+1), stokes))
                continue
            
        # Call MFCLEAN to deconvolve the field and create a clean model
        log_wr(LF, "\n> Deconvolving using MFCLEAN ...")
        imgModData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) \
                      + '_' + str(IFext) + '_' + stokes + '.model'
        if redo and os.path.exists(imgModData):
            log_wr(LF, ">> Deleting previous version first ...")
            shutil.rmtree(imgModData, ignore_errors=True)
        if os.path.exists(imgModData):
            log_wr(LF, ">> Skipping: previous version exists ...")
        else:
            try:
                print mir.mfclean(map=imgData,
                                  beam=imgBeam,
                                  out=imgModData,
                                  gain=pDict['gain'],
                                  cutoff=cutoff,
                                  niters=pDict['niter'],
                                  region='perc(33)',
                                  minpatch=pDict['minpatch'],
                                  speed=pDict['speed'])
            except Exception:
                log_wr(LF, ">> Err: failed to clean '%s_CH%s' in Stokes %s" %
                       (pntName, (i+1), stokes))
                continue
        
        # Convolve the clean components and restore the residual
        log_wr(LF, "\n> Saving a CLEANed map ...")
        imgClnData = pntName + '_' + 'CH%d' % (i+1) + '_' + str(IFext) + \
                      '_' + stokes + '.clean'
        imgClnPath = imageFieldDir + '/' + imgClnData
        outTab[i]['imgClnData'] = imgClnData
        if redo and os.path.exists(imgClnPath):
            log_wr(LF, ">> Deleting previous version first ...")
            shutil.rmtree(imgClnPath, ignore_errors=True)
        if os.path.exists(imgClnPath):
            log_wr(LF, ">> Skipping: previous version exists ...")
        else:
            try:
                restorStr = mir.restor(map=imgData,
                                       beam=imgBeam,
                                       out=imgClnPath,
                                       model=imgModData,
                                       mode='clean',
                                       fwhm=beamFWHM,
                                       pa=beamPA,
                                       options='mfs')
                print restorStr
                beamMaj, beamMin, beamPA = parse_restorStr(restorStr)
                outTab[i]['beamMaj_asec'] = beamMaj
                outTab[i]['beamMin_asec'] = beamMin
                outTab[i]['beamPA_asec'] = beamPA
            except Exception:
                log_wr(LF, ">> Err: failed to restore '%s_CH%s' in Stokes %s" %
                       (pntName, (i+1), stokes))
                continue

        # Measure the intrinsic beam
        if not beamPA == "":
            log_wr(LF, "\n> Measuring the un-forced beam ...")
            tmpClnData = pntName + '_' + 'CH%d' % (i+1) + '_' + str(IFext) + \
                         '_' + stokes + '.temp'
            tmpClnPath = imageFieldDir + '/' + tmpClnData
            if os.path.exists(tmpClnPath):
                shutil.rmtree(tmpClnPath, ignore_errors=True)
            try:
                restorStr = mir.restor(map=imgData,
                                       beam=imgBeam,
                                       out=tmpClnPath,
                                       model=imgModData,
                                       mode='clean',
                                       options='mfs')
                beamMajInt, beamMinInt, beamPAint = parse_restorStr(restorStr)
                outTab[i]['beamMajInt_asec'] = beamMajInt
                outTab[i]['beamMinInt_asec'] = beamMinInt
                outTab[i]['beamPAint_asec'] = beamPAint
                log_wr(LF, ">> Intrinsic beam size is [%s, %s, %s]" % \
                       (beamMajInt, beamMinInt, beamPAint))
            except Exception:
                log_wr(LF, ">> Err: failed to measure the un-forced beam size.")

        # Create a residual map
        log_wr(LF, "\n> Saving a residual map ...")
        imgResData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) \
                      + '_' + str(IFext) + '_' + stokes + '.resid'
        if redo and os.path.exists(imgResData):
            log_wr(LF, ">> Deleting previous version first ...")
            shutil.rmtree(imgResData, ignore_errors=True)
        if os.path.exists(imgResData):
            log_wr(LF, ">> Skipping: previous version exists ...")
        else:
            try:
                print mir.restor(map=imgData,
                                 beam=imgBeam,
                                 out=imgResData,
                                 model=imgModData,
                                 mode='residual')
            except Exception:
                log_wr(LF, ">> Err: failed to residual '%s_CH%s' in Stokes %s" %
                       (pntName, (i+1), stokes))
                continue
    
        # Measure the RMS noise in the centre of the I residual
        log_wr(LF, "\n> Measuring the RMS noise in Stokes I residual")
        outLogTmp = pntName + '_' + 'CH%d' % (i+1) + '_' + str(IFext) + \
                    '.imstatI.dat'
        print mir.imstat(_in=imgResData, log=outLogTmp)
        statDict = parse_imstat(outLogTmp)
        os.remove(outLogTmp)
        log_wr(LF, "\n> RMS=%f" % float(statDict['rms']))
        outTab[i]['rms'] = float(statDict['rms'])
        
        # Stitch the spectrum plots together into a single file
        try:
            specPlt = imageFieldDir + '/' + pntName + '_' + str(IFext) + \
                      '_spec.ps'
            log_wr(LF, "\n> Saving spectra plot to %s" % specPlt)
            cmd = "gs -dNOPAUSE -sDEVICE=ps2write -sOUTPUTFILE=%s " % specPlt
            cmd += "-dBATCH " + " ".join(specPltLst)
            os.system(cmd)
            for e in specPltLst:
                os.remove(e)
        except Exception:
            pass

        # Update the success flag
        outTab[i]['success'] = 1
        
    return outTab
    
    
#-----------------------------------------------------------------------------#
def run_stored_flags(flgDict, flgGrp, vis, LF):
    """
    Run the flagging commands specified in the driving file.
    """
    
    if flgGrp in flgDict:
        log_wr(LF, "\n> Applying FLAG commands in '%s' group ..." % flgGrp)
        for argStr in flgDict[flgGrp]:
            log_wr(LF, "> %s vis=%s" % (argStr, vis))
            argDict = argstr2dict(argStr)
            cmd = argDict.pop('cmd')
            argDict['vis'] = vis
            if cmd=='uvflag':
                if not 'flagval' in argDict:
                    argDict['flagval'] = 'flag'
                print mir.uvflag(**argDict)
            if cmd=='pgflag':
                print mir.pgflag(**argDict)
            if cmd=='mirflag':
                print mir.mirflag(**argDict)


#-----------------------------------------------------------------------------#
def print_usage():
    print "\n\tUSAGE: ./6_imageTileLin.py [-p] [-f <freqExt>] [-d <dayName>]"
    print "\t                                                       <tileID>"
    print "\t       ./6_imageTileLin.py -f 5500 5\n"
    sys.exit()


#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    main()
