#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     2_calObsDay.py                                                    #
#                                                                             #
# USAGE:    ./2_calObsDay.py [-p] <YYYY-MM-DD.freq>                           #
#           ./2_calObsDay.py -p 2010-12-22.5500                               #
#                                                                             #
# PURPOSE:  Calibrate, flag and split out mosaic data from a single day of    #
#           ATCA mosaic observations.                                         #
#                                                                             #
# MODIFIED: 04-Jan-2017 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Hardcoded paths
dataRootDir = '../DATA'

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

from Imports.util_tables import *
from Imports.util_ATCA_PIPE import *
from Imports.util_plotting import *
print "Setting up the MIRpy environment ...",
sys.stdout.flush()
from Imports.mirpy import miriad as mir
print "done."


#-----------------------------------------------------------------------------#
def main():

    # Misc setup
    doPause = False
    startTime = time.time()

    # Check for required files and directories
    startDir = os.getcwd()
    fail_not_exists(dataRootDir, 'directory')
    atlodDir = dataRootDir + "/uvdata_atlod"
    fail_not_exists(atlodDir, 'directory')
    intDataDir =  dataRootDir + "/int_data"
    uvDataSplitDir = dataRootDir + "/uvdata_split"
    reqDirLst = ['LOGS', 'PLOTS', intDataDir, uvDataSplitDir]
    for reqDir in reqDirLst:
        if not os.path.exists(reqDir):
            os.mkdir(reqDir)

    # Get the options from the command line
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'p')
    except getopt.GetoptError:
        print_usage()
    if len(args)< 1:
        print "\nPlease specify a day name on the command line:"
        print_usage()
    dayName, IFext = args[0].rsplit('.', 1)
    for o, a in opts:
        if o == "-p":
            doPause = True

    # Open a log file
    LF = open('LOGS/' + dayName + '.' + IFext + '_cal.log', 'w')
            
    # Read and parse the pipeline input file
    configFile = dayName + '.' + IFext + '.config'
    fail_not_exists(configFile, 'file')
    try:
        pDict, flgDict = config_read_flag(configFile)
        log_wr(LF, "\n> Successfully parsed the input parameter file.")
    except Exception:
        log_fail(LF, "\n> Err: Failed to parse the input parameter file.")
        
    # Check that the multi-source, single-frequency uv-data exists
    uvDataMS = atlodDir + '/' + dayName + '.' + IFext
    fail_not_exists(uvDataMS, 'file', LF)
    
    # Backup the original flag-table, header, var-table and history
    backupDir = uvDataMS + '/tabsOriginal'
    if os.path.exists(backupDir):
        log_wr(LF, "\n> Restoring  original tables ...")
        remLst = ['bandpass', 'flags', 'gains', 'gainsf', 'header', 'history',
                  'leakage', ' leakagef', 'vartable']
        for e in remLst:
            if os.path.exists(uvDataMS + '/' + e):
                os.remove(uvDataMS + '/' + e)
        cpLst = ['flags', 'header', 'history', 'vartable']
        for e in cpLst:
            shutil.copy2(uvDataMS + '/tabsOriginal/' + e, uvDataMS + '/' + e)
    else:
        log_wr(LF, "\n> Backing up original tables ...")
        os.mkdir(backupDir)
        shutil.copy(uvDataMS + '/flags', backupDir + '/flags')
        shutil.copy(uvDataMS + '/header', backupDir + '/header')
        shutil.copy(uvDataMS + '/history', backupDir + '/history')
        shutil.copy(uvDataMS + '/vartable', backupDir + '/vartable')
    
    # Run the stored FLAG commands in the config file [with label: precal]
    run_stored_flags(flgDict, 'precal', uvDataMS, LF)
    
    # Create a directory for the plots from the current day
    outPlotDir = 'PLOTS/' + dayName + '.' + IFext 
    if not os.path.exists(outPlotDir):
        os.mkdir(outPlotDir)

    # Make some initial diagnostic plots on the MS file
    outPS = outPlotDir + '/' + dayName + '.' + IFext + '_uvCov.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the uv-coverage for XX ...')
        mir.uvplt(vis=uvDataMS, select='pol(xx)', axis='uc,vc',
                  options='nobase,equal', device=outPS + '/cps')
    outPS =  outPlotDir + '/' + dayName + '.' + IFext + '_el.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the elevation vs time for XX ...')
        mir.uvplt(vis=uvDataMS, select='pol(xx)', axis='time,el',
                  options='nobase', device=outPS + '/cps')
    outPS =  outPlotDir + '/' + dayName + '.' + IFext + '_ampCalRaw.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the amplitude vs time for XX & YY ...')
        mir.uvplt(vis=uvDataMS, select='pol(xx),pol(yy)', axis='time,amp',
                  nxy='1,1', device=outPS + '/cps')
    outPS = outPlotDir + '/' + dayName + '.' + IFext + '_TsysX.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the Tsys vs time for X ...')
        mir.varplt(vis=uvDataMS, xaxis='time', yaxis='xtsys',
                   device=outPS + '/cps')
    outPS = outPlotDir + '/' + dayName + '.' + IFext + '_TsysY.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the Tsys vs time for Y ...')
        mir.varplt(vis=uvDataMS, xaxis='time', yaxis='ytsys',
                   device=outPS + '/cps')
    outPS = outPlotDir + '/' + dayName + '.' + IFext + '_SMon.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the seeing monitor vs time ...')
        mir.varplt(vis=uvDataMS, xaxis='time', yaxis='smonrms',
                   device=outPS + '/cps')
    outPS = outPlotDir + '/' + dayName + '.' + IFext + '_XYphase.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the XY-phase vs time ...')
        mir.varplt(vis=uvDataMS, xaxis='time', yaxis='xyphase',
                   device=outPS + '/cps')
    outPS = outPlotDir + '/' + dayName + '.' + IFext + '_axisRMS.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the RMS tracking error vs time ...')
        mir.varplt(vis=uvDataMS, xaxis='time', yaxis='axisrms',
                   device=outPS + '/cps')
    outPS = outPlotDir + '/'  + dayName + '.' + IFext + '_axisMax.ps'
    if not os.path.exists(outPS):
        log_wr(LF, '\n> Plotting the max tracking error vs time ...')
        mir.varplt(vis=uvDataMS, xaxis='time', yaxis='axismax',
                   device=outPS + '/cps')

    # Pause
    if doPause:
        pause()

    # Run the calibration routine in two passes -----------------------------#
    # The flagging commands applied before pass 2 are performed on data loaded
    # and calibrated by the first pass. These are specified in the ASCII
    # driving file with a 2 appended.

    # Create a sub-directory for intermediate calibration files
    dataCalDir = intDataDir + '/' + dayName  + '.' + IFext + '_CalDir'
    if not os.path.exists(dataCalDir):
        os.mkdir(dataCalDir)
    
    print
    print '-'*80
    print '>>> CALIBRATION PASS 1'
    print '-'*80
    print

    # Pass 1 of the calibration routine works on uncalibrated data
    dataCalDict = calibrate_atca(uvDataMS, pDict, workingDir=dataCalDir,
                                 outPlotDir=outPlotDir, flgDict=flgDict,
                                 passNum=1, doPlots=True, LF=LF)

    # Pause
    if doPause:
        pause()
    
#    print
#    print '-'*80
#    print '>>> PASS 2'
#    print '-'*80
#    print
    
    # Pass 2 of the calibration routine works on calibrated data
#    dataCalDict = calibrate_atca(uvDataMS, pDict, workingDir=dataCalDir,
#                                 outPlotDir=outPlotDir, flgDict=flgDict,
#                                 passNum=2, doPlots=True, LF=LF)
    
    # Copy from gaincal[0] back to the multisource file
    log_wr(LF, "\n> Copying all calibration tables from '%s' to '%s'" % \
           (dataCalDict['phaseCals'][0], uvDataMS))
    print mir.gpcopy(vis=dataCalDict['phaseCals'][0], out=uvDataMS,
                     mode='copy')

    # Plot the calibrated amps
#    outPS =  outPlotDir + '/' + dayName + '.' + IFext + '_ampCalFlg.ps'
#    if not os.path.exists(outPS):
#        log_wr(LF, '\n> Plotting the amplitude vs time for I ...')
#        mir.uvplt(vis=uvDataMS, select='pol(i)', axis='time,amp',
#                  nxy='1,1', device=outPS + '/cps')
        
    # Copy the flags to the multisource file?
    
    # Split out individual sources / mosaic pointings -----------------------#
    
    print
    print "-"*80
    print
    log_wr(LF, '\n> Ready to split out individual cuts ...')

    # Pause
    if doPause:
        pause()

    # Connect to the database for the current IF/frequency extension
    dbFile = dataRootDir + '/DB_' + str(IFext) + '.sqlite'
    if os.path.exists(dbFile):
        try:
            log_wr(LF, "\n> Connecting to DB file '%s' ..." % dbFile)
            conn = sqlite3.connect(dbFile)
            cursor = conn.cursor()
        except Exception:
            log_fail(LF, "\n> Err: Failed to connect to DB '%s'" % dbFile)
            
    # Fetch a recarray containing the source names and point names
    sql ="SELECT timeStamp, pntName FROM observations where dayName = ?"    
    obsTab = select_into_arr(cursor, sql, args=([dayName]))

    # Loop through and split out the individual times
    for i in range(len(obsTab)):
        startTimeStamp = obsTab['timeStamp'][i]
        try:
            endTimeStamp = obsTab['timeStamp'][i+1]
        except Exception:
            dtObj1 = timestamp_to_datetime(obsTab['timeStamp'][i])
            dtObj2 = dtObj1 + datetime.timedelta(minutes=10)
            endTimeStamp = datetime_to_timestamp(dtObj2)

        # Construct a selection statement
        selStr = 'source(%s)' % (obsTab['pntName'][i])
        selStr += ',time(%s,%s)' % (startTimeStamp, endTimeStamp)
        uvDataSplit = uvDataSplitDir + '/' + startTimeStamp + \
                      '_' + str(IFext) + '.uv'
        
        log_wr(LF, "\n> Splitting '%s' into '%s' ..." % (obsTab['pntName'][i],
                                                         uvDataSplit))
        try:
            shutil.rmtree(uvDataSplit, ignore_errors=True)
            print mir.uvcat(vis=uvDataMS, out=uvDataSplit, options='unflagged',
                            select=selStr)
        except Exception:
            log_wr(LF, "Failed to split!")
            continue
        
        # Update the database flag indicating a successful split
        log_wr(LF, ">>> Updating database with sucessful split ...")
        sql = """
        UPDATE observations
        SET splitFlag=1, useFlag=1
        WHERE timeStamp=?;
        """
        cursor.execute(sql, tuple([startTimeStamp]))
        conn.commit()
        
        # Perform flagging operations on split source
        try:
            log_wr(LF, "\n> Running stored flags ...")
            run_stored_flags(flgDict, 'onsplit', uvDataSplit, LF)
        except Exception:
            log_wr(LF, "Failed to flag!")
            continue

    # Clean up
    cursor.close()
    conn.close()

    # End
    endTime = time.time()    
    print "Duration = %.2f min" % ((endTime-startTime)/60.0)


#-----------------------------------------------------------------------------#
def calibrate_atca(uvDataMS, pDict, workingDir='.', outPlotDir='.',
                   flgDict={}, passNum=1, doPlots=True, LF=None):
    """
    Calibrate an ATCA dataset: bandpass, gains vs time and absolute scale.
    Apply the flagging steps set in the flgDict in between calibration
    operations.
    """

    # Figure out the day name and extension
    dummy, dayIFname = os.path.split(uvDataMS)
    dayName, IFext = os.path.splitext(dayIFname)

    # Flag and calibrate the bandpass ---------------------------------------#
    
    # Split out the bandpass calibrator source
    bpCalSplit = workingDir + '/' +  pDict['bpcal'] + '.BP' + IFext
    if passNum==1:
        log_wr(LF, "\n> Splitting out the bandpass cal '%s' ..." % \
               pDict['bpcal'])
        if os.path.exists(bpCalSplit):
            log_wr(LF, ">>> Removing previous version ...")
            shutil.rmtree(bpCalSplit, ignore_errors=True)
        print mir.uvcat(vis=uvDataMS, out=bpCalSplit,
                        options='nocal,nopol,nopass,unflagged',
                        select='source(%s)' % pDict['bpcal'])
    else:
        log_wr(LF, "\n> Using previously split bandpass cal '%s' ..." % \
               pDict['bpcal'])

    # Make an un-calibrated backup
    if not os.path.exists(bpCalSplit + '.bak'):
        shutil.copytree(bpCalSplit, bpCalSplit + '.bak')

    # Run the stored FLAG commands in the config file [with label: bpcal<N>]
    run_stored_flags(flgDict, 'bpcal' + str(passNum), bpCalSplit, LF)

    # Run MFCAL to calibrate the bandpass
    # Run GPCAL to calibrate the gains so we can copy to the MS file later
    log_wr(LF, "\n> Calibrating the bandpass solutions ...")
    print mir.mfcal(vis=bpCalSplit, interval=pDict['intervalBP'],
              refant=pDict['refant'], options='delay')
    print mir.gpcal(vis=bpCalSplit, interval=pDict['intervalBP'],
                    refant=pDict['refant'], options='xyvary,reset',
                    nfbin=pDict['nfbinBP'])

    if doPlots:
        
        # Plot the bandpass solutions
        outPS = outPlotDir + '/' + dayIFname + "_P%d" % passNum + \
                '_SolnBPamp.ps'
        print mir.gpplt(vis=bpCalSplit, yaxis='amp', options='bandpass',
                        yrange='0,2', nxy='1,1', device=outPS + '/cps')
        outPS = outPlotDir + '/' + dayIFname + "_P%d" % passNum + \
                '_SolnBPpha.ps'
        print mir.gpplt(vis=bpCalSplit, yaxis='pha', options='bandpass',
                        nxy='1,1', yrange='-180,180', device=outPS + '/cps')
        
        # Plot the spectrum of the BP calibrator
        outPS = outPlotDir + '/' + pDict['bpcal'] + '.BP' + IFext + \
                "_P%d" % passNum + '_BPcalSpec.ps'
        if os.path.exists(outPS):
            os.remove(outPS)
        print mir.uvspec(vis=bpCalSplit, stokes='xx,yy', interval=1e3,
                         nxy='1,1', axis='frequency,amplitude',
                         device=outPS + '/cps')

        # Plot the Real versus Imaginary visibilities
        outPS = outPlotDir + '/' + pDict['bpcal'] + '.BP' + IFext + \
                "_P%d" % passNum + '_BPcalReIm.ps'
        if os.path.exists(outPS):
            os.remove(outPS)
        print mir.uvplt(vis=bpCalSplit, stokes='xx,yy', axis='real,imag',
                            options='source,avall,nobase', nxy='1,1',
                            device=outPS + '/cps')

    # Flag and calibrate the Phase Calibrator(s) ----------------------------#

    # Make a list out of single phasecals
    if not isinstance(pDict['phasecals'], list):
        pDict['phasecals'] = [pDict['phasecals']]

    # Split out and flag the phase calibrator(s)
    phCalSplitLst = []
    for phCalName in pDict['phasecals']:
        phCalSplit = workingDir + '/' +  phCalName + '.GAIN' + IFext
        if passNum==1:
            log_wr(LF, "\n> Splitting out the Phase-cal '%s' ..." % phCalName)
            phCalSplitLst.append(phCalSplit)
            if os.path.exists(phCalSplit):
                log_wr(LF, ">>> Removing previous version ...")
                shutil.rmtree(phCalSplit, ignore_errors=True)
            print mir.uvcat(vis=uvDataMS, out=phCalSplit,
                            options='nocal,nopol,nopass,unflagged',
                            select='source(%s)' % phCalName)

            # Make an un-calibrated backup
            if not os.path.exists(phCalSplit + '.bak'):
                shutil.copytree(phCalSplit, phCalSplit + '.bak')
                
            # Copy over the bandpass table
            print mir.gpcopy(vis=bpCalSplit, out=phCalSplit, mode='copy',
                             options='nocal,nopol')

        else:
            log_wr(LF, "\n> Using previously split phase cal '%s' ..." % \
                   phCalName)
            phCalSplitLst.append(phCalSplit)
            
        # Run the stored FLAG commands in the config file [label: phcal]
        run_stored_flags(flgDict, 'phcal' + str(passNum), phCalSplit, LF)

    # Calibrate each phasecal in turn 
    for phCalName, phCalSplit in zip(pDict['phasecals'], phCalSplitLst):
        log_wr(LF, "\n> Calibrating gains-vs-time on '%s' ..." % phCalName)
        print mir.gpcal(vis=phCalSplit, interval=pDict['intervalGC'],
                        refant=pDict['refant'], options='xyvary,qusolve',
                        nfbin=pDict['nfbin'])

    # Merge the solutions from other calibrators into the 1st phasecal
    if len(pDict['phasecals'])>1:
        log_wr(LF, "\n> Merging solutions into the 1st gain-calibrator")
        for i in range(1, len(pDict['phasecals'])):
            mir.gpcopy(vis=phCalSplitLst[i], out=phCalSplitLst[0],
                       mode='merge', options='nopass')
            
    if doPlots:
        
        # Plot the XX & YY solutions
        log_wr(LF, "\n> Plotting the solutions")
        outAmpLog = outPlotDir + '/' + dayIFname + "_P%d" % passNum + \
                    '_SolnGainAmp'
        print mir.gpplt(vis=phCalSplitLst[0], yaxis='amp', options='gains',
                        yrange='0,2', log=outAmpLog + '.dat')
        timeArr_hrs, gainsXXArr, gainsYYArr, gainType, polType = \
                     parse_gainlog(outAmpLog + '.dat')
        figXX = plot_gains_vs_time(timeArr_hrs, gainsXXArr, gainType, "XX")
        figXX.savefig(outAmpLog + 'XX.eps')
        figYY = plot_gains_vs_time(timeArr_hrs, gainsYYArr, gainType, "YY")    
        figYY.savefig(outAmpLog + 'YY.eps')
        os.remove(outAmpLog + '.dat')
        outPhaLog = outPlotDir + '/' + dayIFname + "_P%d" % passNum + \
                    '_SolnGainPha'
        print mir.gpplt(vis=phCalSplitLst[0], yaxis='pha', options='wrap',
                        yrange='-180,180', log=outPhaLog + '.dat')
        timeArr_hrs, gainsXXArr, gainsYYArr, gainType, polType = \
                     parse_gainlog(outPhaLog + '.dat')
        figXX = plot_gains_vs_time(timeArr_hrs, gainsXXArr, gainType, "XX")
        figXX.savefig(outPhaLog + 'XX.eps')
        figYY = plot_gains_vs_time(timeArr_hrs, gainsYYArr, gainType, "YY")
        figYY.savefig(outPhaLog + 'YY.eps')
        os.remove(outPhaLog + '.dat')
        
        # Plot the XY solutions
        outAmpLog = outPlotDir + '/' + dayIFname + "_P%d" % passNum + \
                    '_SolnGainAmpXY'
        print mir.gpplt(vis=phCalSplitLst[0], yaxis='amp', options='xygains',
                        yrange='0,2', log=outAmpLog + '.dat')
        timeArr_hrs, gainsXXArr, gainsYYArr, gainType, polType = \
                     parse_gainlog(outAmpLog + '.dat')
        figXY = plot_gains_vs_time(timeArr_hrs, gainsXXArr, gainType, "XY")
        figXY.savefig(outAmpLog + '.eps')
        os.remove(outAmpLog + '.dat')
        outPhaLog = outPlotDir + '/' + dayIFname + "_P%d" % passNum + \
                    '_SolnGainPhaXY'
        print mir.gpplt(vis=phCalSplitLst[0], yaxis='pha', options='xygains',
                        yrange='-180,180', log=outPhaLog + '.dat')
        timeArr_hrs, gainsXXArr, gainsYYArr, gainType, polType = \
                     parse_gainlog(outPhaLog + '.dat')
        figXY = plot_gains_vs_time(timeArr_hrs, gainsXXArr, gainType, "XY")
        figXY.savefig(outPhaLog + '.eps')
        os.remove(outPhaLog + '.dat')
        
        # Plot the phasecals after applying the solutions
        for phCalName, phCalSplit in zip(pDict['phasecals'], phCalSplitLst):
            log_wr(LF, "\n> Plotting the calibrated '%s' ..." % phCalName)
            outPS = outPlotDir + '/' + phCalName + IFext + "_P%d" % passNum +\
                    '.GAIN_Amp.ps'
            print mir.uvplt(vis=phCalSplit, stokes='xx,yy', axis='time,amp',
                            options='source', nxy='1,3',
                            device=outPS + '/cps')
            outPS = outPlotDir + '/' + phCalName + IFext + "_P%d" % passNum +\
                    '.GAIN_Pha.ps'
            print mir.uvplt(vis=phCalSplit, stokes='xx,yy', axis='time,pha',
                            options='source', yrange='-80,80', nxy='1,3',
                            device=outPS + '/cps')
            
            # Plot the spectrum of the PHASE calibrator
            outPS = outPlotDir + '/' + phCalName + IFext + \
                    "_P%d" % passNum + '_PcalSpec.ps'
            if os.path.exists(outPS):
                os.remove(outPS)
            print mir.uvspec(vis=phCalSplit, stokes='xx,yy', interval=1e3,
                             nxy='1,1', axis='frequency,amplitude',
                             device=outPS + '/cps')
            
    # Flag and calibrate the Flux Calibrator(s) -----------------------------#

    # Split out the flux calibrator
    fluxCalSplit = workingDir + '/' + pDict['fluxcal'] + '.FLUX' + IFext
    if passNum==1:
        log_wr(LF, "\n> Splitting out the flux cal '%s' ..." %  \
               pDict['fluxcal'])
        if os.path.exists(fluxCalSplit):
            log_wr(LF, ">>> Removing previous version ...")
            shutil.rmtree(fluxCalSplit, ignore_errors=True)
        print mir.uvcat(vis=uvDataMS, out=fluxCalSplit,
                        options='nocal,nopol,nopass,unflagged',
                        select='source(%s)' % pDict['fluxcal'])

        # Copy over the bandpass table
        print mir.gpcopy(vis=bpCalSplit, out=fluxCalSplit, mode='copy',
                         options='nocal,nopol')
    else:
        log_wr(LF, "\n> Using previously split flux cal '%s' ..." % \
               pDict['fluxcal'])

    # Run the stored FLAG commands in the config file [with label: fluxcal1]
    run_stored_flags(flgDict, 'fluxcal' + str(passNum), fluxCalSplit, LF)

    # Calibrate the primary and bootstrap
    log_wr(LF, "\n> Bootstrapping the flux scale.")
    print mir.gpcal(vis=fluxCalSplit, interval=pDict['intervalFC'],
                    refant=pDict['refant'], options='xyvary',
                    nfbin=pDict['nfbin'])
    print mir.gpboot(vis=phCalSplitLst[0], cal=fluxCalSplit)
    print mir.mfboot(vis="%s,%s" % (fluxCalSplit, phCalSplitLst[0]),
                     select="source(%s)" % pDict['fluxcal'])

    # Split out a test calibrator and plot the spectral shape ----------------#
    if 'testcal' in pDict:        
        testCalSplit = workingDir + '/' + pDict['testcal'] + '.TEST' + IFext
        if passNum==1:
            log_wr(LF, "\n> Splitting out the test cal '%s' ..." %  \
                   pDict['testcal'])
            if os.path.exists(testCalSplit):
                log_wr(LF, ">>> Removing previous version ...")
                shutil.rmtree(testCalSplit, ignore_errors=True)
            print mir.uvcat(vis=uvDataMS, out=testCalSplit,
                            options='nocal,nopol,nopass,unflagged',
                            select='source(%s)' % pDict['testcal'])
        else:
            log_wr(LF, "\n> Using previously split test cal '%s' ..." % \
                   pDict['fluxcal'])            

        # Copy over the bandpass table
        print mir.gpcopy(vis=bpCalSplit, out=testCalSplit, mode='copy')

        # Plot the spectrum of the TEST calibrator
        outPS = outPlotDir + '/' + pDict['testcal'] + IFext + \
                "_P%d" % passNum + '_TESTspec.ps'
        if os.path.exists(outPS):
            os.remove(outPS)
        print mir.uvspec(vis=testCalSplit, stokes='xx,yy', interval=1e3,
                         nxy='1,1', axis='frequency,amplitude',
                         device=outPS + '/cps')

    # Return the locations of the bandpass, phase and flux calibrators
    # Calibrations tables may be copied from the phCalSplitLst[0]
    # Flag tables may be copied from all calibrators
    dataDict = {}
    dataDict['bpCal'] = bpCalSplit
    dataDict['phaseCals'] = phCalSplitLst
    dataDict['fluxCal'] = fluxCalSplit
    return dataDict


#-----------------------------------------------------------------------------#
def print_usage():
    print "\n\tUSAGE: ./calObsDay.py [-p] <YYYY-MM-DD.freq>\n"
    sys.exit()
    

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
if __name__ == '__main__':
    main()
