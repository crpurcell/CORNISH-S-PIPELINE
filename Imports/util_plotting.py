#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plotting.py                                                  #
#                                                                             #
# PURPOSE:  Helper functions and classes for plotting in the ATCA data        #
#           reduction pipeline.                                               #
#                                                                             #
# MODIFIED: 30-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#           parse_gainlog        ... parse the log produced by GPPLT          #
#           parse_BPgainlog      ... parse the time BP log produced by GPPLT  #
#           plot_gains_vs_time   ... plot the gain soulutions from GPPLT      #
#           plot_BPwaterfallse   ... plot the BP-time soulutions from GPPLT   #
#           make_nice_ax         ... format axis label padding and ticks      #
#                                                                             #
#=============================================================================#
import os
import sys
import re
import numpy as np
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import LineCollection


#-----------------------------------------------------------------------------#
def parse_gainlog(inFileName, verbose=False):
    """
    Parse the log produced by GPPLT. Return 3D arrays of gain for XX and YY
    feeds with structure [ant, fbin, gain] and an array of times in hours. If
    the XY phase has been plotted return XY phase in the XX array.
    """
    
    gainType = None
    polType = None

    # List to hold the time axis
    timeLst_hrs = []

    # Lists to hold the Gains
    gainXX1Lst = []
    gainXX2Lst = []
    gainXX3Lst = []
    gainXX4Lst = []
    gainXX5Lst = []
    gainXX6Lst = []
    gainYY1Lst = []
    gainYY2Lst = []
    gainYY3Lst = []
    gainYY4Lst = []
    gainYY5Lst = []
    gainYY6Lst = []
    
    # temporary arrays
    xx1Lst = []
    yy1Lst = []
    xx2Lst = []
    yy2Lst = []
    xx3Lst = []
    yy3Lst = []
    xx4Lst = []
    yy4Lst = []
    xx5Lst = []
    yy5Lst = []
    xx6Lst = []
    yy6Lst = []
    
    # Define regular expressions to match new FBins & the two entry lines
    newNFbinRe = re.compile("^# Listing of the (\S+) of the gains for (\S+)")
    ent1Re = re.compile("^\s(\d)\s(\d{2}:\d{2}:\d{2})\s+(\S+)\s+(\S+)" +
                        "\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
    ent2Re = re.compile("^(\s{5})(\s{5}\s+)(\S+)\s+(\S+)\s+(\S+)\s+(\S+)" +
                        "\s+(\S+)\s+(\S+)")

    # Loop through the lines in the file
    freqBinCnt = 0
    FH = open(inFileName, 'r')
    for line in FH:        
        line = line.rstrip("\r\n")

        # Trigger a match for new freq bin
        mch = newNFbinRe.match(line)
        if mch:
            gainType = mch.group(1)
            polType = mch.group(2)            
            freqBinCnt += 1
            print "Reading frequency bin %d" % freqBinCnt

            # Append the previous bin's reading to the aggregate list
            if freqBinCnt>1:                
                gainXX1Lst.append(xx1Lst)
                gainXX2Lst.append(xx2Lst)
                gainXX3Lst.append(xx3Lst)
                gainXX4Lst.append(xx4Lst)
                gainXX5Lst.append(xx5Lst)
                gainXX6Lst.append(xx6Lst)
                
                if polType!='XY':
                    gainYY1Lst.append(yy1Lst)
                    gainYY2Lst.append(yy2Lst)
                    gainYY3Lst.append(yy3Lst)
                    gainYY4Lst.append(yy4Lst)
                    gainYY5Lst.append(yy5Lst)
                    gainYY6Lst.append(yy6Lst)

            # Starting a new bin, so reset the temporary lists
            xx1Lst = []
            yy1Lst = []
            xx2Lst = []
            yy2Lst = []
            xx3Lst = []
            yy3Lst = []
            xx4Lst = []
            yy4Lst = []
            xx5Lst = []
            yy5Lst = []
            xx6Lst = []
            yy6Lst = []
            
        # Trigger a match on the first line of an entry
        mch = ent1Re.match(line)
        if mch:
            if polType=='XY':
                xx1Lst.append(float(mch.group(3)))
                xx2Lst.append(float(mch.group(4)))
                xx3Lst.append(float(mch.group(5)))
                xx4Lst.append(float(mch.group(6)))
                xx5Lst.append(float(mch.group(7)))
                xx6Lst.append(float(mch.group(8)))                
            else:
                xx1Lst.append(float(mch.group(3)))
                yy1Lst.append(float(mch.group(4)))
                xx2Lst.append(float(mch.group(5)))
                yy2Lst.append(float(mch.group(6)))
                xx3Lst.append(float(mch.group(7)))
                yy3Lst.append(float(mch.group(8)))

            # Populate the time list from the first freq bin
            if freqBinCnt==1:
                days = float(mch.group(1))
                time = mch.group(2).split(':')
                time = [float(x) for x in time]
                time_hrs = days*24.0 + time[0] + time[1]/60. + time[2]/3600.
                timeLst_hrs.append(time_hrs)

        # Trigger a match on the second line of an entry
        mch = ent2Re.match(line)
        if mch:
            xx4Lst.append(float(mch.group(3)))
            yy4Lst.append(float(mch.group(4)))
            xx5Lst.append(float(mch.group(5)))
            yy5Lst.append(float(mch.group(6)))
            xx6Lst.append(float(mch.group(7)))
            yy6Lst.append(float(mch.group(8)))

    # Populate the last frequency bin (not triggered)
    if polType=='XY':
        gainXX1Lst.append(xx1Lst)
        gainXX2Lst.append(xx2Lst)
        gainXX3Lst.append(xx3Lst)
        gainXX4Lst.append(xx4Lst)
        gainXX5Lst.append(xx5Lst)
        gainXX6Lst.append(xx6Lst)
    else:
        gainXX1Lst.append(xx1Lst)
        gainYY1Lst.append(yy1Lst)
        gainXX2Lst.append(xx2Lst)
        gainYY2Lst.append(yy2Lst)
        gainXX3Lst.append(xx3Lst)
        gainYY3Lst.append(yy3Lst)
        gainXX4Lst.append(xx4Lst)
        gainYY4Lst.append(yy4Lst)
        gainXX5Lst.append(xx5Lst)
        gainYY5Lst.append(yy5Lst)
        gainXX6Lst.append(xx6Lst)
        gainYY6Lst.append(yy6Lst)
                
    # Clean up
    FH.close()

    # Format the lists into arrays [antenna,bin,time]
    gainsXXArr =  np.array([gainXX1Lst, gainXX2Lst, gainXX3Lst,
                            gainXX4Lst, gainXX5Lst, gainXX6Lst])
    gainsYYArr =  np.array([gainYY1Lst, gainYY2Lst, gainYY3Lst,
                            gainYY4Lst, gainYY5Lst, gainYY6Lst])
    timeArr_hrs = np.array(timeLst_hrs)
    
    return timeArr_hrs, gainsXXArr, gainsYYArr, gainType, polType


#-----------------------------------------------------------------------------#
def parse_BPgainlog(inFileName, verbose=False):
    """
    Parse the log produced by GPPLT of time-varying bandpass. Return 3D arrays
    of passband gains for XX and YY feeds with structure [ant, soln#, gain]
    and an array of frequencies in Hz.
    """
    
    gainType = None

    # List to hold the time axis
    freqLst_Hz = []

    # Lists to hold the Gains
    gainXX1Lst = []
    gainXX2Lst = []
    gainXX3Lst = []
    gainXX4Lst = []
    gainXX5Lst = []
    gainXX6Lst = []
    gainYY1Lst = []
    gainYY2Lst = []
    gainYY3Lst = []
    gainYY4Lst = []
    gainYY5Lst = []
    gainYY6Lst = []
    
    # temporary arrays
    xx1Lst = []
    yy1Lst = []
    xx2Lst = []
    yy2Lst = []
    xx3Lst = []
    yy3Lst = []
    xx4Lst = []
    yy4Lst = []
    xx5Lst = []
    yy5Lst = []
    xx6Lst = []
    yy6Lst = []
    
    # Define regular expressions to match new FBins & the two entry lines
    newSolnRe = re.compile("^# Listing of the (\S+)\s+")
    ent1Re = re.compile("^\s+(\d+\.\d+)\s+(\S+)\s+(\S+)" +
                        "\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
    ent2Re = re.compile("^\s{10,13}(\S+)\s+(\S+)\s+(\S+)\s+(\S+)" +
                        "\s+(\S+)\s+(\S+)")

    # Loop through the lines in the file
    solnCnt = 0
    FH = open(inFileName, 'r')
    for line in FH:        
        line = line.rstrip("\r\n")

        # Trigger a match for new solution
        mch = newSolnRe.match(line)
        if mch:
            gainType = mch.group(1)
            solnCnt += 1
            print "Reading solution %d" % solnCnt
            
            # Append the previous bin's reading to the aggregate list
            if solnCnt>1:
                gainXX1Lst.append(xx1Lst)
                gainYY1Lst.append(yy1Lst)
                gainXX2Lst.append(xx2Lst)
                gainYY2Lst.append(yy2Lst)
                gainXX3Lst.append(xx3Lst)
                gainYY3Lst.append(yy3Lst)
                gainXX4Lst.append(xx4Lst)
                gainYY4Lst.append(yy4Lst)
                gainXX5Lst.append(xx5Lst)
                gainYY5Lst.append(yy5Lst)
                gainXX6Lst.append(xx6Lst)
                gainYY6Lst.append(yy6Lst)

            # Starting a new bin, so reset the temporary lists
            xx1Lst = []
            yy1Lst = []
            xx2Lst = []
            yy2Lst = []
            xx3Lst = []
            yy3Lst = []
            xx4Lst = []
            yy4Lst = []
            xx5Lst = []
            yy5Lst = []
            xx6Lst = []
            yy6Lst = []
            
        # Trigger a match on the first line of an entry
        mch = ent1Re.match(line)
        if mch:
            xx1Lst.append(float(mch.group(2)))
            yy1Lst.append(float(mch.group(3)))
            xx2Lst.append(float(mch.group(4)))
            yy2Lst.append(float(mch.group(5)))
            xx3Lst.append(float(mch.group(6)))
            yy3Lst.append(float(mch.group(7)))

            # Populate the frequnecy list from the first solution
            if solnCnt==1:
                freqLst_Hz.append(float(mch.group(1)))

        # Trigger a match on the second line of an entry
        mch = ent2Re.match(line)
        if mch:
            xx4Lst.append(float(mch.group(1)))
            yy4Lst.append(float(mch.group(2)))
            xx5Lst.append(float(mch.group(3)))
            yy5Lst.append(float(mch.group(4)))
            xx6Lst.append(float(mch.group(5)))
            yy6Lst.append(float(mch.group(6)))

    # Populate the last solution (not triggered)
    gainXX1Lst.append(xx1Lst)
    gainYY1Lst.append(yy1Lst)
    gainXX2Lst.append(xx2Lst)
    gainYY2Lst.append(yy2Lst)
    gainXX3Lst.append(xx3Lst)
    gainYY3Lst.append(yy3Lst)
    gainXX4Lst.append(xx4Lst)
    gainYY4Lst.append(yy4Lst)
    gainXX5Lst.append(xx5Lst)
    gainYY5Lst.append(yy5Lst)
    gainXX6Lst.append(xx6Lst)
    gainYY6Lst.append(yy6Lst)
                
    # Clean up
    FH.close()

    # Format the lists into arrays [antenna,sol,chan]
    gainsXXArr = np.array([gainXX1Lst, gainXX2Lst, gainXX3Lst,
                            gainXX4Lst, gainXX5Lst, gainXX6Lst])
    #print gainsXXArr.shape

    gainsXXArr = np.where(gainsXXArr==0.0, np.nan, gainsXXArr)
    gainsYYArr = np.array([gainYY1Lst, gainYY2Lst, gainYY3Lst,
                            gainYY4Lst, gainYY5Lst, gainYY6Lst])
    gainsYYArr = np.where(gainsYYArr==0.0, np.nan, gainsYYArr)
    freqArr_Hz = np.array(freqLst_Hz)

    
    return freqArr_Hz, gainsXXArr, gainsYYArr, gainType


#-----------------------------------------------------------------------------#
def plot_gains_vs_time(timeArr_hrs, gainsArr, gainType, feed="XX"):
    """
    Create a figure showing Amplitude or Phase gains versus time from the
    logfile output by GPPLT. Returns a MPL figure object.
    """

    # Setup the figure page
    fig = plt.figure(figsize=(12.0, 8.0))
    
    # One plot per antenna
    for antN in range(gainsArr.shape[0]):

        # Axis placement
        xi = antN%3
        yi = antN/3
        ax1 = fig.add_subplot(2,3,antN+1)

        # Assemble the linecollection
        lineLst = []
        for binN in range(gainsArr.shape[1]):
            lineLst.append(gainsArr[antN, binN, :])
        lineSegs = LineCollection([list(zip(timeArr_hrs, y))
                                   for y in lineLst],
                                  linewidths = 2,
                                  linestyles = 'solid')
        lineSegs.set_array(np.linspace(0,gainsArr.shape[1]-1,
                                       gainsArr.shape[1]))
        ax1.add_collection(lineSegs)
        if xi==2:
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.00)
            cbar = fig.colorbar(lineSegs, cax=cax)
            cbar.set_label('Frequency Bin')

        # Formatting
        if gainType=='Phase':
            ax1.set_ylim(-180.0, 180.0)
        else:
            gainType='Amplitude'
            ax1.set_ylim(np.nanmin(gainsArr), np.nanmax(gainsArr))
        ax1.set_xlim(np.nanmin(timeArr_hrs), np.nanmax(timeArr_hrs))
        ax1 = make_nice_ax(ax1, pad=7)            
        ax1.set_title('Antenna %d - %s' % (antN+1, feed))        
        if (yi==0):
            ax1.axes.xaxis.set_ticklabels([])
        else:
            ax1.set_xlabel("Time (hrs)")
        if not xi==0:
            ax1.axes.yaxis.set_ticklabels([])
        else:
            ax1.set_ylabel(gainType)
            
        fig.subplots_adjust(left=0.08, right=0.94, bottom=0.08, top=0.94,
                            wspace=0.08,hspace=0.11)
        
    return fig


#-----------------------------------------------------------------------------#
def plot_BPwaterfalls(freqArr_Hz, gainsArr, gainType, feed="XX",
                      saveFits=False, outDir="."):
    """
    Create a figure showing Amplitude or Phase gains versus time from the
    logfile output by GPPLT. Returns a MPL figure object.
    """

    # Setup the figure page
    fig = plt.figure(figsize=(8.0, 10.0))

    # Determine colourscale range
    med = nanmedian(gainsArr)
    mad = MAD(gainsArr.flatten())
    vMin = med-mad*3
    vMax = med+mad*3

    # One plot per antenna
    for antN in range(gainsArr.shape[0]):
        
        data = gainsArr[antN,:,:]
        aspect = float(data.shape[0]) / float(data.shape[1])
        
        # Axis placement
        ax = fig.add_subplot(6,1,antN+1)
        im = ax.imshow(data, interpolation="Nearest", origin='lower',
                       aspect = aspect, vmax=vMax, vmin=vMin,
                       extent=[freqArr_Hz[0], freqArr_Hz[-1],
                               1.0, gainsArr.shape[-2]],
                       cmap=plt.get_cmap('jet'))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.set_ylabel("Ant %d" % (antN+1))
        if antN==0:
            ax.set_title("Time-varying Bandpass Solutions (%s) - %s feed" \
                         % (gainType, feed))
        if antN==5:
            ax.set_xlabel("Frequency (GHz)")
            
        # Add a colourbar
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax)

        # Save a FITS file
        if not saveFits is None:
            outFile = outDir + "/BPsoln_%s_" % gainType
            outFile += "Time_Ant%d_%s.fits" % ((antN+1),feed)
            pf.writeto(outFile, data, clobber=True)

    return fig


#-----------------------------------------------------------------------------#
def make_nice_ax(ax, pad=10):
    """
    Do some formatting on a MPL axis
    """
    
    ax.tick_params(pad=pad)
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markeredgewidth(1)
    return ax
