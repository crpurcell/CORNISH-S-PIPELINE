#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     4_genEquTileGrid.py                                               #
#                                                                             #
# USAGE:    ./4_genEquTileGrid.py                                             #
#                                                                             #
# PURPOSE:  Read the pointing coordinates from the mosaic files and lay down  #
#           grid of Equatorial tiles for image testing purposes.              #
#                                                                             #
# MODIFIED: 05-Jan-2017 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Hardcoded paths
dataRootDir = "../DATA"

# CORNISH-South border out to edge of data in Galactic coordinates
# Make this large enough to account for data out to the primary beam FWHM
bMax_deg = +1.2
bMin_deg = -1.2
lMax_deg = 350.2
lMin_deg = 294.8

# Tile parameters 
imSize_px = [4000, 4000]              # pixels [x, y] tile size
pixScale_deg = [0.30/3600, 0.30/3600]   # pixel scale [dx, dy]
overlap_deg  = [60.0/3600, 60.0/3600]   # overlap between tiles [x, y] (deg)

#-----------------------------------------------------------------------------#

import os
import sys
import copy
import glob
import re
import math as m
import numpy as np
from pyslalib import slalib
import pylab as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Ellipse, RegularPolygon, Polygon, Patch
from matplotlib.collections import PatchCollection
import sqlite3

from Imports.util_ATCA_PIPE import sort_nicely

# Constants
C = 2.998e8

#-----------------------------------------------------------------------------#
def main():

    # Create a polygon describing the Galactic border of the survey
    # Oversample each edge and combine into an ordered set of vertices
    lBorLst_deg = np.linspace(lMin_deg, lMax_deg, 5500).tolist()
    bBorLst_deg = np.linspace(bMin_deg, bMax_deg, 220).tolist()
    borderPolyGalLst = zip(lBorLst_deg, [bMin_deg]*len(lBorLst_deg))
    borderPolyGalLst += zip([lMax_deg]*len(bBorLst_deg), bBorLst_deg,)[1:]
    borderPolyGalLst += zip(lBorLst_deg[::-1], [bMax_deg]*len(lBorLst_deg))[1:]
    borderPolyGalLst += zip([lMin_deg]*len(bBorLst_deg), bBorLst_deg[::-1])
    borderPolyGalArr = np.array(borderPolyGalLst)
    lRange_deg = lMax_deg - lMin_deg
    bRange_deg = bMax_deg - bMin_deg
    
    # Convert the Galactic polygon vertices into Equatorial coordinates and
    # determine the maximum and minimum RA and Dec. limits
    borderPolyEquLst = []
    for e in borderPolyGalLst:
        ra_rad, dec_rad = slalib.sla_galeq(m.radians(e[0]), m.radians(e[1]))
        borderPolyEquLst.append( (m.degrees(ra_rad), m.degrees(dec_rad)) )
    borderPolyEquArr = np.array(borderPolyEquLst)
    raMax_deg = np.max(borderPolyEquArr[:,0])
    raMin_deg = np.min(borderPolyEquArr[:,0])
    decMax_deg = np.max(borderPolyEquArr[:,1])
    decMin_deg = np.min(borderPolyEquArr[:,1])
    raRange_deg = raMax_deg - raMin_deg
    decRange_deg = decMax_deg - decMin_deg

    # Calculate the constant Dec (y) increment between tile centres
    yIncr_deg = imSize_px[1] * pixScale_deg[1] - overlap_deg[1]

    #------------------------------------------------------------------------#
    # NOTE:
    # Start at the bottom-left of the Equ grid and advance along a Dec. line
    # setting down tiles. Starting coordinate = decMin_deg, raMin_deg.
    # Note: Placing tiles on lines of constant Dec does not take into account
    # the curvature of the Dec lines as we approach the equatorial pole,
    # however, it should be good enough if the overlap between the tiles is
    # enough and the cos(Dec) factor is calculated at the most negative Dec.
    #------------------------------------------------------------------------#
    raCentTileLst_deg = []
    decCentTileLst_deg = []
    vertexTileEquLst_deg = []
    vertexTileGalLst_deg = []
    
    # Loop through Dec rows until decMax reached
    i = 0
    while True:
        
    	# Calculate the Dec at the centre top and bottom of the current row
	decTileCent_deg = decMin_deg + (yIncr_deg - 2 * overlap_deg[1]) * i
        decTileTop_deg = decTileCent_deg - yIncr_deg/2.0
        decTileBot_deg = decTileCent_deg + yIncr_deg/2.0
        
	# Calculate the RA increment for this row
	cosDecCent = m.cos(m.radians(decTileCent_deg))
	cosDecTop = m.cos(m.radians(decTileTop_deg))
	cosDecBot = m.cos(m.radians(decTileBot_deg))
        cosDec = min(cosDecCent, cosDecTop, cosDecBot)
	xIncr_deg = (imSize_px[0] * pixScale_deg[0] - 2*overlap_deg[0])/cosDec
        i += 1
        
	# Loop through the RAs until raMax reached
        j = 0
        while True:

            # Calculate RA for this tile
            raTileCent_deg = raMin_deg + xIncr_deg * j
            raCentTileLst_deg.append(raTileCent_deg)
            decCentTileLst_deg.append(decTileCent_deg)
            j += 1

            # Calculate the tile corner coorinates in Equ
            xIncrTop_deg = (imSize_px[0] * pixScale_deg[0])/cosDecTop
            xIncrBot_deg = (imSize_px[0] * pixScale_deg[0])/cosDecBot
            x1y2 = (raTileCent_deg + xIncrTop_deg / 2.0, decTileTop_deg)
            x2y2 = (raTileCent_deg - xIncrTop_deg / 2.0, decTileTop_deg)
            x2y1 = (raTileCent_deg - xIncrBot_deg / 2.0, decTileBot_deg)
            x1y1 = (raTileCent_deg + xIncrBot_deg / 2.0, decTileBot_deg)
            vertexTileEquLst_deg.append(np.array([x1y1,x1y2,x2y2,x2y1]))

            # Calculate the tile corner coordinates in Gal
            lV_rad, bV_rad = slalib.sla_eqgal(m.radians(x1y2[0]),
                                              m.radians(x1y2[1]))
            x1y2 = (m.degrees(lV_rad), m.degrees(bV_rad))
            lV_rad, bV_rad = slalib.sla_eqgal(m.radians(x2y2[0]),
                                              m.radians(x2y2[1]))
            x2y2 = (m.degrees(lV_rad), m.degrees(bV_rad))
            lV_rad, bV_rad = slalib.sla_eqgal(m.radians(x2y1[0]),
                                              m.radians(x2y1[1]))
            x2y1 = (m.degrees(lV_rad), m.degrees(bV_rad))
            lV_rad, bV_rad = slalib.sla_eqgal(m.radians(x1y1[0]),
                                              m.radians(x1y1[1]))
            x1y1 = (m.degrees(lV_rad), m.degrees(bV_rad))
            vertexTileGalLst_deg.append(np.array([x1y1,x1y2,x2y2,x2y1]))

            # End of RA While loop
            if raTileCent_deg>=raMax_deg:
                break

        # End of Dec While loop
        if decTileCent_deg>=decMax_deg:
            break

    # Convert the tile centre coordinates to Galactic
    lCentTileLst_deg = []
    bCentTileLst_deg = []
    for i in range(len(raCentTileLst_deg)):
        l_rad, b_rad = slalib.sla_eqgal(m.radians(raCentTileLst_deg[i]),
                                        m.radians(decCentTileLst_deg[i]))
        lCentTileLst_deg.append(m.degrees(l_rad))
        bCentTileLst_deg.append(m.degrees(b_rad))

    # Filter both Equ and Gal lists for tiles outside the survey area
    # Must iterate from highest index when using 'pop' function
    for i in range(len(raCentTileLst_deg)-1, -1, -1):
        
        if not (lCentTileLst_deg[i]>=lMin_deg and
                lCentTileLst_deg[i]<=lMax_deg and
                bCentTileLst_deg[i]>=bMin_deg and
                bCentTileLst_deg[i]<=bMax_deg):
            lCentTileLst_deg.pop(i)
            bCentTileLst_deg.pop(i)
            raCentTileLst_deg.pop(i)
            decCentTileLst_deg.pop(i)
            vertexTileEquLst_deg.pop(i)
            vertexTileGalLst_deg.pop(i)
            
    # Sort the list of tiles into increasing RA
    multiLst = zip(raCentTileLst_deg,
                   decCentTileLst_deg,
                   lCentTileLst_deg,
                   bCentTileLst_deg)
    multiLst.sort()
    (raCentTileLst_deg,
     decCentTileLst_deg,
     lCentTileLst_deg,
     bCentTileLst_deg) = zip(*multiLst)

    # Create the remaining columns (ID, pixel-scale, num-pixels
    tileIDLst = range(1, len(lCentTileLst_deg)+1)
    pixScaleXLst_asec = [pixScale_deg[0]*3600.0] * len(lCentTileLst_deg)
    pixScaleYLst_asec = [pixScale_deg[1]*3600.0] * len(lCentTileLst_deg)
    nPixXLst = [imSize_px[0]]* len(lCentTileLst_deg)
    nPixYLst = [imSize_px[1]]* len(lCentTileLst_deg)

    #------------------------------------------------------------------------#

    # Upload the tile parameters into each database file in the data directory
    dbFileLst = glob.glob(dataRootDir + '/*.sqlite')

    # Loop through the database files
    for dbFile in dbFileLst:
        print ">> Writing tile_coords table to %s ..." % dbFile
        
        # Connect to the database
        conn = sqlite3.connect(dbFile)
        cursor = conn.cursor()
        
        # Drop old tile_coords table and create a new one
        sql = "DROP TABLE IF EXISTS tile_coords"
        cursor.execute(sql)
        sql = """
        CREATE TABLE tile_coords (
        tileID INTEGER PRIMARY KEY, 
        RA_deg DOUBLE,
        Dec_deg DOUBLE,
        l_deg DOUBLE,
        b_deg DOUBLE,
        pixScaleX_asec DOUBLE,
        pixscaleY_asec DOUBLE,
        nPixX INTEGER,
        nPixY INTEGER );
        """
        cursor.execute(sql)

        # Insert the entries into the table
        for i in range(len(raCentTileLst_deg)):
        
            sql = """
            INSERT INTO tile_coords
            (tileID,
            RA_deg,
            Dec_deg,
            l_deg,
            b_deg,
            pixScaleX_asec,
            pixscaleY_asec,
            nPixX,
            nPixY)
            """        
            sql += 'VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?) '
            vals = [tileIDLst[i],
                    raCentTileLst_deg[i],
                    decCentTileLst_deg[i],
                    lCentTileLst_deg[i],
                    bCentTileLst_deg[i],
                    pixScaleXLst_asec[i],
                    pixScaleYLst_asec[i],
                    nPixXLst[i],
                    nPixYLst[i]]
            cursor.execute(sql, vals)
        
        # Commit changed to the database and close connection
        conn.commit()
        cursor.close()
        conn.close()
    
    #------------------------------------------------------------------------#
    #------------------------------------------------------------------------#
    
    # Plot the tile centres over the survey border
    fig = plt.figure(figsize=(18.0, 10.0))    

    # EQUATORIAL PLOT -------------------------------------------------------#
    ax1 = fig.add_axes([0.08, 0.4, 0.88, 0.58])

    # Plot the tile centres and vertices
    ax1.scatter(np.array(raCentTileLst_deg)/15.0, decCentTileLst_deg, s=2,
                zorder=2)
    squares = []
    for vertex in vertexTileEquLst_deg:
        vertex = np.array(vertex)
        vertex[:,0]/=15.0
        square = Polygon(xy=np.array(vertex), closed=True)
        squares.append(square)
        s = PatchCollection(squares, alpha=1.0, edgecolor='black',
                            facecolor='none', zorder=3)
    ax1.add_collection(s)

    # Plot the border and format the axis
    ax1.plot(borderPolyEquArr[:,0]/15.0, borderPolyEquArr[:,1])
    ax1.yaxis.grid(True, which='major')
    ax1.xaxis.grid(True, which='major')
    ax1.set_xlim((raMax_deg+0.01*raRange_deg)/15,
                 (raMin_deg-0.01*raRange_deg)/15)
    ax1.set_ylim(decMin_deg-0.05*decRange_deg, decMax_deg+0.05*decRange_deg)
    ax1.set_aspect(1.0/15.0/cosDec)
    ax1.set_ylabel('Dec. (deg)')
    ax1.set_xlabel('R.A. (hrs)')

    # Annotate the Equatorial plot with tile numbers
#    for i in range(len(raCentTileLst_deg)):
#        ax1.annotate(str(tileIDLst[i]),
#                     xy=(raCentTileLst_deg[i]/15.0, decCentTileLst_deg[i]),
#                     horizontalalignment='center',
#                     verticalalignment='center',
#                     fontsize=6,
#                     textcoords='data',
#                     clip_on=True, backgroundcolor='w')
    
    # GALACTIC PLOT ---------------------------------------------------------#
    ax2 = fig.add_axes([0.08, 0.05, 0.88, 0.30])
    ax2.plot(borderPolyGalArr[:,0], borderPolyGalArr[:,1])

    # Plot the tile centres and vertices
    #ax2.scatter(np.array(lCentTileLst_deg), bCentTileLst_deg, s=2,
    #            zorder=2)
    squares = []
    for vertex in vertexTileGalLst_deg:
        square = Polygon(xy=np.array(vertex), closed=True)
        squares.append(square)
        s = PatchCollection(squares, alpha=1.0, edgecolor='black',
                            facecolor='none', zorder=3)
    ax2.add_collection(s)
    
    # Plot the border and format the axis
    ax2.yaxis.set_major_locator(MaxNLocator(4))
    ax2.yaxis.grid(True, which='major')
    ax2.xaxis.grid(True, which='major')
    ax2.set_xlim(lMax_deg+0.02*lRange_deg, lMin_deg-0.02*lRange_deg)
    ax2.set_ylim(bMin_deg-0.19*bRange_deg, bMax_deg+0.19*bRange_deg)
    ax2.set_aspect(1.0)
    ax2.set_ylabel('Glong. (deg)')
    ax2.set_xlabel('Glat. (deg)')
    
    # Annotate the Galactic plot with tile numbers
    for i in range(len(lCentTileLst_deg)):
        ax2.annotate(str(tileIDLst[i]),
                     xy=(lCentTileLst_deg[i], bCentTileLst_deg[i]),
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=8,
                     textcoords='data',
                     clip_on=True)

    fig.show()
    fig.savefig('tile_layout.pdf')
    
    print "Press <RETURN> to exit ..."
    raw_input()
    
    
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    main()



