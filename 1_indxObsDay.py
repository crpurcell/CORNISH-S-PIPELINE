#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     1_indxObsDay.py                                                   #
#                                                                             #
# USAGE:    $>./1_indxObsDay.py <YYYY-MM-DD.freq>                             #
#             ./1_indxObsDay.py 2010-12-22.5500                               #
#                                                                             #
# PURPOSE:  Index uv-data from a single day of ATCA observations and update   #
#           the database with the observation parameters. The database is     #
#           created if necesary using the sqlite3 module from the schema      #
#           defined in 'Imports/DBSchema.sql'. Apply any initial gross        #
#           flagging commands in 'YYYY-MM-DD.freq' file, if it exists.        #
#           These flagging commands should have a label 'onindex'.            #
#                                                                             #
# MODIFIED: 03-Jan-2017 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Hardcoded paths
dataRootDir = "../DATA"

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
import numpy.lib.recfunctions as rec
import sqlite3

from Imports.util_ATCA_PIPE import *
from Imports.util_tables import *
register_sqlite3_numpy_dtypes()   # Map numpy data-types to sqlite3 data-types

# Setup MIRpy 
print "Setting up the MIRpy environment ...",
sys.stdout.flush()
from Imports.mirpy import miriad as mir
print "done."

#-----------------------------------------------------------------------------#
def main():
    
    startTime = time.time()

    # Get the dayName from the command line
    args = sys.argv[1:]
    if len(args)<1:
        print_usage()
    dayName, IFext = args[0].rsplit('.', 1)

    # Open a log file
    if not os.path.exists('LOGS'):
        os.mkdir('LOGS')
    LF = open('LOGS/' + dayName + '.' + IFext + '_indx.log', 'w')
    
    # Check for required files and directories
    startDir = os.getcwd()
    fail_not_exists(dataRootDir, 'directory', LF)
    atlodDir = dataRootDir + "/uvdata_atlod"
    fail_not_exists(atlodDir, 'directory', LF)
    
    # Parse the schema file for the table column definitions and
    # their python equivalent dtypes
    try:
        schemaFile = 'Imports/DBSchema.sql'
        tableDtypeDict, tableSQLdict = schema_to_dtypes(schemaFile)
    except Exception:
        log_fail(LF, "Err: Failed to parse the DB schema '%s'." % schemaFile)
    
    # Check that the ATLOD data exists
    uvDataAtlod = atlodDir + '/' + dayName + '.' + IFext
    fail_not_exists(uvDataAtlod, 'file', LF)

    # Parse any initial flag commands
    flgDict = {}
    configFile = dayName + '.' + IFext + '.config'
    if os.path.exists(configFile):
        dummy, flgDict = config_read_flag(configFile)

    # Backup the original ATLOD flag table, header, vartable and history
    backupDir = uvDataAtlod + '/tabsOriginal'
    if not os.path.exists(backupDir):
        log_wr(LF, "Backing up original ATLOD tables ...")
        os.mkdir(backupDir)
        shutil.copy(uvDataAtlod + '/flags', backupDir + '/flags')
        shutil.copy(uvDataAtlod + '/header', backupDir + '/header')
        shutil.copy(uvDataAtlod + '/history', backupDir + '/history')
        shutil.copy(uvDataAtlod + '/vartable', backupDir + '/vartable')

    # Apply any global flags (not frequency specific) in .config file
    # Flag commands must be tagged with <flagcommand.>onindex='...'
    run_stored_flags(flgDict, 'onindex', uvDataAtlod, LF)

    # Create or connect to a database for each frequency extension
    dbFile = dataRootDir + '/DB_' + str(IFext) + '.sqlite'
    if os.path.exists(dbFile):
        try:
            log_wr(LF, "Connecting to existing DB file '%s' ..." % dbFile)
            conn = sqlite3.connect(dbFile)
            cursor = conn.cursor()
        except Exception:
            log_fail(LF, "Err: Failed to connect to DB '%s'" % dbFile)
    else:
        try:
            log_wr(LF, "Creating DB file '%s' using '%s' ..." % (dbFile,
                                                                 schemaFile))
            conn = sqlite3.connect(dbFile)
            cursor = conn.cursor()
            for tableName, sql in tableSQLdict.iteritems():
                cursor.execute(sql)
            conn.commit()
        except Exception:
            os.remove(dbFile)
            log_fail(LF, "Err: Failed to create database from schema!")

    # Run UVINDEX and parse to find observation parameters
    log_wr(LF, "Indexing %s" % (dayName + '.' + IFext))
    uvindexLog = atlodDir + '/' +  dayName + '.' + IFext + '.uvindex'
    print mir.uvindex(vis=uvDataAtlod, log=uvindexLog)
    log_wr(LF, "Parsing UVINDEX log ...")
    timeTab, pntTab, freqTab, fieldTab = parse_uvindex_new(uvindexLog)

    # Add the dayName column to the recarray
    timeTab['dayName'] = dayName
        
    # Update/insert the DB entries
    colNameLst = ['timeStamp', 'dayName', 'pntName', 'calCode', 'nAnts',
                  'fieldName', 'RA_deg', 'Dec_deg']
    insert_arr_db(cursor, timeTab, 'observations', colNameLst)
    colNameLst = ['fieldName', 'RA_deg', 'Dec_deg', 'calCode']
    insert_arr_db(cursor, fieldTab, 'field_coords', colNameLst)
    colNameLst = ['pntName', 'calCode', 'RA_deg', 'Dec_deg', 'fieldName']
    insert_arr_db(cursor, pntTab, 'pointings', colNameLst)
    conn.commit()

    # Query the distinct pointing entries to fill the mosaic_obs table
    sql = """
    SELECT DISTINCT dayName || '.' || pntName as dayPntID,
    dayName,
    pntName,
    calCode,
    RA_deg,
    Dec_deg,
    0 as splitFlag,
    0 as useFlag
    FROM observations
    WHERE dayName = ?;
    """
    mosaicTab = select_into_arr(cursor, sql, args=([dayName]))
    insert_arr_db(cursor, mosaicTab, 'mosaic_obs')
    conn.commit()

    endTime = time.time()
    log_wr(LF, "Duration = %.2f min" % ((endTime-startTime)/60.0))
    
    # Clean up
    cursor.close()
    conn.close()


#-----------------------------------------------------------------------------#
def print_usage():
    print "\n\tUSAGE: 1_indexObsDay.py <YYYY-MM-DD.freq>\n"
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
