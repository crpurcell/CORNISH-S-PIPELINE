#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     3_imageField.py                                                   #
#                                                                             #
# USAGE:    ./3_imageField.py [-p] [-f <freqExt>] [-d <dayName>] <pntName>    #
#           ./3_imageField.py -p -f 5500 173315.19-372232.40 # 1729-37        #
#           ./3_imageField.py -p -f 5500 171141.33-384224.77                  #
#           ./3_imageField.py -p -f 5500 -d 2010-12-22 171141.33-384224.77    #
#                                                                             #
# PURPOSE:  Image a field from an ATCA dataset                                #
#                                                                             #
# NOTE:     If <freqExt> is not supplied then it defaults to '5500'.          #
#                                                                             #
# MODIFIED: 04-Jan-2017 by C. Purcell                                         #
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

#-----------------------------------------------------------------------------#

import os
import sys
import shutil
import copy
import re
import optparse
import getopt
import time
import math as m
import numpy as np
import sqlite3

from Imports.util_ATCA_PIPE import *
from Imports.util_tables import regexp
print "Setting up the MIRpy environment ...",
sys.stdout.flush()
from Imports.mirpy import miriad as mir
print "done."


#-----------------------------------------------------------------------------#
def main():
    
    # Misc setup
    IFext = defIFext
    doPause = False
    dayName = None
    startTime = time.time()
    
    # Check for required files and directories
    startDir = os.getcwd()
    fail_not_exists(dataRootDir, 'directory')
    uvDataSplitDir = dataRootDir + "/uvdata_split"
    fail_not_exists(uvDataSplitDir, 'directory')
    imageFieldDir = dataRootDir + "/images_fields"
    reqDirLst = ['LOGS', 'PLOTS', imageFieldDir]
    for reqDir in reqDirLst:
        if not os.path.exists(reqDir):
            os.mkdir(reqDir)

    # Get the options from the command line
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'pf:d:')
    except getopt.GetoptError:
        print_usage()
        sys.exit()
    if len(args)< 1:
        print "\nPlease specify a field name on the command line:"
        print_usage()
        sys.exit()
    pntName = args[0]
    for o, a in opts:
        if o == "-f":
            IFext = a
        if o == "-d":
            dayName = a
        if o == "-p":
            doPause = True
    log_wr(None, "\n> Using the frequency extension '%s'." % IFext)

    # Open a log file
    LF = open('LOGS/' + pntName + '.' + IFext + '_img.log', 'w')

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

    # Query the database for the uv-data in the field
    obsTab = query_obs_in_field(cursor=cursor, fieldName=pntName,
                                splitFlag=True, useFlag=True, dayName=dayName)
    if len(obsTab)==0:
        log_fail(LF, "\n> Err: No entries for this field in DB.")

    # Change the timestamps column to a list of paths to the uvsplit data
    # Format as a string which MIRIAD can read.
    uvSplitDataLst = obsTab['timeStamp'].tolist()
    uvSplitDataLst = [e + '_' + str(IFext) + '.uv' for e in uvSplitDataLst]
    uvSplitStr = ','.join(uvSplitDataLst)
    
    # Feedback:
    message = ">Found the following data for field '%s':\n" % pntName
    for e in uvSplitDataLst:
        message += "%s\n" % e
    log_wr(LF, message)

    # Concatenate the uv-data into a single dataset
    # Need to move to the directory containing the uv-data as MIRIAD will
    # complain if the vis input string is too long. Move back after uvcat.
    os.chdir(uvDataSplitDir)
    uvDataCat = pntName+ '_' + str(IFext) + '.uv'
    if os.path.exists(uvDataCat):
        shutil.rmtree(uvDataCat, ignore_errors=True)
    log_wr(LF, ">Concatenating data into single file '%s'"% uvDataCat) 
    print uvDataSplitDir
    print mir.uvcat(vis=uvSplitStr, out=uvDataCat, options='unflagged')
    uvDataCat = uvDataSplitDir + '/' + uvDataCat
    os.chdir(startDir)

    # Check for custom imaging configuration file
    customCfg =  'customImgConfigs/' + pntName + '_img' + IFext + '.config'
    log_wr(LF, "\n> Checking for custom config:\n  '%s' ..." % customCfg)
    if os.path.exists(customCfg):
        log_wr(LF, "  Found! Replacing default entries:")
        pDictCustom, flgTabCustom = config_read_flag(customCfg)
        for e, v in pDictCustom.iteritems():
            pDict[e] = v
            log_wr(LF, "  %s = %s" % (e, v))
        log_wr(LF, "  Adding flagging commands:")
        for k, v in flgTabCustom.iteritems():
            log_wr(LF, "   Group '%s':" % k)
            for e in v:
                log_wr(LF, "   '%s'" % e)
            if k in flgTab:
                flgTab[k] = flgTab[k] + v
            else:
                flgTab[k] = v

    # Run pre-imaging flags
    run_stored_flags(flgTab, 'preimg', uvDataCat, LF)
                
    # Calculate the default image parameters
    outLogTmp = pntName + '_' + IFext + '.uvlist.dat'
    mir.uvlist(vis=uvDataCat, options='array', log=outLogTmp)
    baselineDict = parse_uvlist_array(outLogTmp)
    os.remove(outLogTmp)
    outLogTmp = pntName + '_' + IFext + '.uvlist.dat'
    mir.uvlist(vis=uvDataCat, options='spectral', log=outLogTmp)
    specDict = parse_uvlist_spectral(outLogTmp)
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
    
    # Create a primary beam annotation file
    ANN = open(imageFieldDir + '/' + pntName + '_' + str(IFext) + \
               '.ann', 'w')
    ANN.write("COLOUR MAGENTA\n")
    ANN.write("CIRCLE W %f %f %f\n" % (obsTab['RA_deg'][0],
                                       obsTab['Dec_deg'][0],
                                       FWHMmin_deg/2.0 ))
    ANN.write("COLOUR GREEN\n")
    ANN.write("CIRCLE W %f %f %f\n" % (obsTab['RA_deg'][0],
                                       obsTab['Dec_deg'][0],
                                       FWHMmax_deg/2.0 ))
    ANN.write("COLOUR YELLOW\n")
    ANN.write("CIRCLE W %f %f %f\n" % (obsTab['RA_deg'][0],
                                       obsTab['Dec_deg'][0],
                                       (FWHMmax_deg 
                                        * float(pDict['fov_FWHM'])
                                        /2.0 / 3.0 )))
    ANN.close()
    
    # Set the image resolution parameters: 3 pixels across a beam,
    # Image size specified in the driving parameter file in terms of the
    # primary beam at the lowest frequency
    if not 'cellsize_asec' in pDict:
        pDict['cellsize_asec'] = resolutionHigh_asec / 3.0
    FWHMmax_pix = int(m.ceil(FWHMmax_asec/float(pDict['cellsize_asec'])))
    imsize_pix = int(m.ceil(FWHMmax_pix * float(pDict['fov_FWHM'])))
    log_wr(LF, "\n> Imaging parameters:")
    log_wr(LF, "  Pixel size: %.2f arcsec"  % float(pDict['cellsize_asec']))
    log_wr(LF, "  Image size: %d pixels (%.2f arcmin)"  % \
           (imsize_pix, (imsize_pix*float(pDict['cellsize_asec'])/60.0)))

    # Loop through the chosen number of image bins and use the line selection
    # parameter to choose and average channels to image
    specPltLst = []
    nchanBin = int(m.floor(specDict['nChan']/int(pDict['nchanImg'])))
    for i in range(int(pDict['nchanImg'])):
        
        # Line selection
        lineStr = "channel,%d,%d,%d" % (nchanBin,(1+i*nchanBin),1)

        # Plot the spectrum of each chunk
        specPlt = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                  '_' + str(IFext) + '_spec.ps'
        if os.path.exists(specPlt):
            os.remove(specPlt)
        log_wr(LF, "> PLOTTING: '%s'" % specPlt)
        print mir.uvspec(vis=uvDataCat, line=lineStr, stokes='xx,yy',
                         interval=1e3, options='nobase,avall', nxy='1,1',
                         device=specPlt + '/cps')
        if os.path.exists(specPlt):
            specPltLst.append(specPlt)
            
        # Call invert to create a Stokes V dirty image
        log_wr(LF, "\n> Imaging a Stokes V patch.")
        imgVData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                   '_' + str(IFext) + '_V.xy'
        if os.path.exists(imgVData):
            shutil.rmtree(imgVData, ignore_errors=True)
        print mir.invert(vis=uvDataCat,
                         map=imgVData,
                         imsize=imsize_pix,
                         cell=pDict['cellsize_asec'],
                         robust=pDict['robust'],
                         options='mfs,sdb,mosaic',
                         stokes='v',
                         line=lineStr)

        # Measure the RMS noise in the centre of the V image
        log_wr(LF, "\n> Measuring the RMS noise in Stokes V")
        outLogTmp = pntName + '_' + IFext + '.imstat.dat'
        print mir.imstat(_in=imgVData, log=outLogTmp)
        statDict = parse_imstat(outLogTmp)
        os.remove(outLogTmp)
        cutoff = float(pDict['rms_mult']) * float(statDict['rms'])
        log_wr(LF, "\n> RMS=%f, MULT=%.2f, CUTOFF=%f" % \
               (float(statDict['rms']), float(pDict['rms_mult']), cutoff))
    
        # Call invert to create a dirty Stokes I spectral beam
        log_wr(LF, "\n> Creating a dirty image ...")
        imgIData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                   '_' + str(IFext) + '_I.xy'
        if os.path.exists(imgIData):
            shutil.rmtree(imgIData, ignore_errors=True)
        imgIBeam = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                   '_' + str(IFext)  + '_I.beam'
        if os.path.exists(imgIBeam):
            shutil.rmtree(imgIBeam, ignore_errors=True)
        print mir.invert(vis=uvDataCat,
                         map=imgIData,
                         beam=imgIBeam,
                         imsize=imsize_pix,
                         cell=pDict['cellsize_asec'],
                         robust=pDict['robust'],
                         options='mfs,sdb',
                         stokes='i',
                         line=lineStr)
    
        # Call MFCLEAN to deconvolve the field and create a clean model
        log_wr(LF, "\n> Deconvolving using MFCLEAN ...")
        imgIModData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                      '_' + str(IFext) + '_I.model'
        if os.path.exists(imgIModData):
            shutil.rmtree(imgIModData, ignore_errors=True)
        print mir.mfclean(map=imgIData,
                          beam=imgIBeam,
                          out=imgIModData,
                          gain=pDict['gain'],
                          cutoff=cutoff,
                          niters=pDict['niter'],
                          region='perc(33)',
                          minpatch=pDict['minpatch'],
                          speed=pDict['speed'])
    
        # Convolve the clean components and restore the residual
        log_wr(LF, "\n> Saving a CLEANed map ...")
        imgIClnData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                      '_' + str(IFext) + '_I.clean'
        if os.path.exists(imgIClnData):
            shutil.rmtree(imgIClnData, ignore_errors=True)
        print mir.restor(map=imgIData,
                         beam=imgIBeam,
                         out=imgIClnData,
                         model=imgIModData,
                         mode='clean',
                         options='mfs')

        # Create a residual map
        log_wr(LF, "\n> Saving a residual map ...")
        imgIResData = imageFieldDir + '/' + pntName + '_' + 'CH%d' % (i+1) + \
                      '_' + str(IFext) + '_I.resid'
        if os.path.exists(imgIResData):
            shutil.rmtree(imgIResData, ignore_errors=True)
        print mir.restor(map=imgIData,
                         beam=imgIBeam,
                         out=imgIResData,
                         model=imgIModData,
                         mode='residual')

    # Stitch the spec plots together into a single file
    specPlt = imageFieldDir + '/' + pntName + '_' + str(IFext) + '_spec.ps'
    log_wr(LF, "\n> Saving spectra plot to %s" % specPlt)
    cmd = "gs -dNOPAUSE -sDEVICE=ps2write -sOUTPUTFILE=%s " % specPlt
    cmd += "-dBATCH " + " ".join(specPltLst)
    os.system(cmd)
    for e in specPltLst:
        os.remove(e)
    
    
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
    print "\n\tUSAGE: ./4_imageField.py [-f <IFext>] [-p] <pntName>"
    print "\t       ./4_imageField.py -f 5500 -p 171052.72-381135.35\n"


#-----------------------------------------------------------------------------#
def pause():
    print "\nPress <RETURN> to continue ...",
    raw_input()


#-----------------------------------------------------------------------------#
if __name__ == '__main__':
    main()
