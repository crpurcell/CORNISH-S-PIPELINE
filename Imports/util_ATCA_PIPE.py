#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_ATCA_PIPE.py                                                 #
#                                                                             #
# PURPOSE:  Helper functions for the ATCA data reduction pipeline.            #
#                                                                             #
# REQUIRED: Requires numpy.                                                   #
#                                                                             #
# MODIFIED: 21-May-2015 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  run_command          ... run a system command                              #
#  argstr2dict          ... convert space-delimited key=val string to dict    #
#  config_read          ... read a key=value format text file                 #
#  config_read_flag     ... read a key=value format text file (flag version)  #
#  deg2dms              ... convert decimal degrees to dms string             #
#  fail_not_exists      ... check for dir/file existence, exit if missing     #
#  log_fail             ... record a fatal error in the log file and exit     #
#  log_wr               ... record a message in the log file                  #
#  nanmedian            ... np.median ignoring NaNs                           #
#  MAD                  ... calculate the madfm                               #
#  calc_stats           ... calculate the statistics of an array              #
#  calc_clipped_stats   ... calculate the stats after sigma-clipping          #
#  parse_uvindex        ... parse the UVINDEX log file                        #
#  query_obs_in_field   ... query the uv-data belonging to a field            #
#  query_obs_in_pointing .. query the uv-data belonging to a mosaic pointing  #
#  query_tile_parameters .. query the parameters of a tile a field            #
#  query_fields_in_tile ... query the fields overlapping a tile               #
#  query_pointings_in_tile  query the mosaic pointings overlapping a tile     #
#  parse_uvlist_spectral .. parse the UVLIST log for channel information      #
#  parse_uvlist_array   ... parse the UVLIST log for baseline information     #
#  parse_imstat         ... parse the IMSTAT log                              #
#  make_lin_axis        ... create linear array from FITS axis key values     #
#  sort_nicely          ... sort a list in the order a human would            #
#  merge_dicts          ... marge any number of dictionaries                  #
#  pause                ... pause using raw input                             #
#                                                                             #
#=============================================================================#

import os
import sys
import re
import time
import datetime
import math as m
import numpy as np
import numpy.ma as ma

from util_tables import *


#-----------------------------------------------------------------------------#
def run_command(command):
    """
    Execute a system command in the $PATH.
    """
    
    # Clean up the spaces in the command
    spaces = re.compile('\s+')
    command = command.strip()
    command = spaces.sub(' ', command)
    
    # Print the command to screen and execute
    print "-"*80
    print command
    print "-"*80
    sys.stdout.flush()
    return os.system(command)

    
#-----------------------------------------------------------------------------#
def argstr2dict(argStr, delim=" "):
    """
    Parse a argument string into a dictionary, e.g.,
    argStr = 'in=myfilein select=source(g19.5) device=/xs'
    """
    argStrLst = argStr.split(delim)
    argDict = {}
    for arg in argStrLst:
        keyVal = arg.split('=')
        argDict[keyVal[0]] = keyVal[1]

    return argDict

    
#-----------------------------------------------------------------------------#
def config_read(filename, delim='=', doValueSplit=True):
    """
    Read a configuration file and output a 'KEY=VALUE' dictionary.
    """

    configTable = {}
    CONFIGFILE = open(filename, "r")
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    commaOrSpace = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    keyVal = re.compile('^.+' + delim + '.+')
    key_ext = re.compile('^.+\.\d+')
    uvflg = re.compile('')

    # Read in the input file, line by line
    for line in CONFIGFILE:

        valueLst=[]
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line) and keyVal.match(line):

            # Weed out internal comments & split on 1st space
            line = comment.sub('',line)
            (keyword, value) = line.split(delim,1)

            # If the line contains a value
            keyword = keyword.strip()              # kill external whitespace
            keyword = spaces.sub('', keyword)      # kill internal whitespaces
            value = value.strip()                  # kill external whitespace
            value = spaces.sub(' ', value)         # shrink internal whitespace
            value = value.replace("'", '')         # kill quotes
            value = commaAndSpaces.sub(',', value) # kill ambiguous spaces

            # Split comma/space delimited value strings
            if doValueSplit:
                valueLst = commaOrSpace.split(value)
                if len(valueLst)<=1:
                    valueLst = valueLst[0]
                configTable[keyword] = valueLst
            else:
                configTable[keyword] = value

    return configTable


#-----------------------------------------------------------------------------#
def config_read_flag(filename, delim='=', doValueSplit=True):
    """
    Read a configuration file and output a 'KEY=VALUE' dictionary.
    Parse flagging commands (uvflag, pgflag & mirflag) separately.
    """

    configTable = {}
    flgTable = {}
    CONFIGFILE = open(filename, "r")
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    commaOrSpace = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    keyVal = re.compile('^.+' + delim + '.+')
    key_ext = re.compile('^.+\.(\w+)')
    uvflag = re.compile('^uvflag')
    pgflag = re.compile('^pgflag')
    mirflag = re.compile('^mirflag')

    # Read in the input file, line by line
    for line in CONFIGFILE:

        valueLst=[]
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line) and keyVal.match(line):

            # Weed out internal comments & split on 1st space
            line = comment.sub('',line)
            (keyword, value) = line.split(delim,1)

            # If the line contains a value
            keyword = keyword.strip()              # kill external whitespace
            keyword = spaces.sub('', keyword)      # kill internal whitespaces
            value = value.strip()                  # kill external whitespace
            value = spaces.sub(' ', value)         # shrink internal whitespace
            value = commaAndSpaces.sub(',', value) # kill ambiguous spaces
            value = value.replace("'", '')         # kill quotes
            
            # Check for flag commands
            flagCmd = None
            if uvflag.match(keyword):
                flagCmd = 'uvflag'
            if pgflag.match(keyword):
                flagCmd = 'pgflag'
            if mirflag.match(keyword):
                flagCmd = 'mirflag'
            if not flagCmd is None:
                mch = key_ext.match(keyword)
                if mch:
                    flgExt = mch.group(1)
                else:
                    flgExt = '0'
                value = "cmd=%s %s" % (flagCmd, value)
                if flgExt in flgTable:
                    flgTable[flgExt].append(value)
                else:
                    flgTable[flgExt] = [value]
                continue
                
            # Split comma/space delimited value strings
            if doValueSplit:
                valueLst = commaOrSpace.split(value)
                if len(valueLst)<=1:
                    valueLst = valueLst[0]
                configTable[keyword] = valueLst
            else:
                configTable[keyword] = value

    return configTable, flgTable


#-----------------------------------------------------------------------------#
def dms2deg(dms):
    
    try:
        delim = re.compile('[,| |:|h|d|m|s]')
        dmsLst = delim.split(dms)
        d = float(dmsLst[0])
        m = float(dmsLst[1])
        s = float(dmsLst[2])
        sign = 1
        if d != 0.0:
            sign = int( d/abs(d) )
        return sign * ( abs(d) + m / 60.0 + s / 3600.0 )
    except Exception:
        return None


#-----------------------------------------------------------------------------#
def deg2dms(deg, delim=':', doSign=False, nPlaces=2):
    """
    Convert a float in degrees to 'dd mm ss' format.
    """

    try:
        angle = abs(deg)
        sign=1
        if angle!=0: sign = angle/deg
        
        # Calcuate the degrees, min and sec
        dd = int(angle)
        rmndr = 60.0*(angle - dd)
        mm = int(rmndr)
        ss = 60.0*(rmndr-mm)

        # If rounding up to 60, carry to the next term
        if float("%05.2f" % ss) >=60.0:
            mm+=1.0
            ss = ss - 60.0
        if float("%02d" % mm) >=60.0:
            dd+=1.0
            mm = mm -60.0
        if nPlaces> 0:
            formatCode = "%0" + "%s.%sf" % (str(2 + nPlaces + 1), str(nPlaces))
        else:
            formatCode = "%02.0f"
        if sign>0:
            if doSign:
                formatCode = "+%02d%s%02d%s" + formatCode
            else:
                formatCode = "%02d%s%02d%s" + formatCode
        else:
            formatCode = "-%02d%s%02d%s" + formatCode
        return formatCode % (dd, delim, mm, delim, ss)
        
    except Exception:
        return None

    
#-----------------------------------------------------------------------------#
def fail_not_exists(item, type='file', LF=None):
    """
    Check for a file or directory, exit if not found. Messages to STDOUT or a
    log file via a file handle.
    """
    if LF is None:
        LF = sys.stdout
    if not os.path.exists(item):
        log_fail(LF, "Err: The %s '%s' does not exist." % (type, item))
    else:
        log_wr(LF, "Found the %s '%s'." % (type, item))


#-----------------------------------------------------------------------------#
def log_fail(LF, errStr):
    """
    Record a fatal error in the log file, echo to STDOUT and exit.
    """
    
    if LF is None:
        LF = sys.stdout
    else:
        print errStr
    if not LF==sys.stdout:
        LF.write(errStr + '\n')
        LF.close()
    sys.exit(1)


#-----------------------------------------------------------------------------#
def log_wr(LF, message):
    """
    Record a message in a log file and echo to STDOUT.
    """
    
    if LF is None:
        LF = sys.stdout
    else:
        print message
    if not LF==sys.stdout:
        LF.write(message + '\n')
    sys.stdout.flush()
    


#-----------------------------------------------------------------------------#
def nanmedian(arr, **kwargs):
    """
    Returns median ignoring NANs.
    """
    
    return ma.median( ma.masked_where(arr!=arr, arr), **kwargs )


#-----------------------------------------------------------------------------#
def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std

    """
    
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m


#-----------------------------------------------------------------------------#
def calc_stats(a, maskzero=False):
    """
    Calculate the statistics of an array (masked version).
    """
    
    statsDict = {}
    a = np.array(a)

    # Mask off bad values and count valid pixels
    if maskzero:
        a = np.where( np.equal(a, 0.0), np.nan, a)
    am = ma.masked_invalid(a)
    statsDict['npix'] = np.sum(~am.mask)
    
    if statsDict['npix']>=2:
        statsDict['stdev'] = float(np.std(am))
        statsDict['mean'] = float(np.mean(am))
        statsDict['median'] = float(nanmedian(am))
        statsDict['max'] = float(np.max(am))
        statsDict['min'] = float(np.min(am))
        statsDict['centmax'] = list(np.unravel_index(np.argmax(am),
                                                     a.shape))
        statsDict['madfm'] = float(MAD(am.flatten()))
        statsDict['success'] = True
        
    else:
        statsDict['npix'] == 0
        statsDict['stdev']   = 0.0
        statsDict['mean']    = 0.0
        statsDict['median']  = 0.0
        statsDict['max']     = 0.0
        statsDict['min']     = 0.0
        statsDict['centmax'] = (0.0, 0.0)
        statsDict['madfm']   = 0.0
        statsDict['success'] = False
        
    return statsDict


#-----------------------------------------------------------------------------#
def calc_clipped_stats(data, clip=3.0, nIter=10, maskzero=False):
    """
    Calculate the statistics of an array given a sigma clip.
    """
    
    if maskzero:
        data = np.where( np.equal(data, 0.0), np.nan, data)
        
    ms = calc_stats(data)

    if clip>0 and nIter>0:
        convergeFlg = 0
        itCnt = 0
        while convergeFlg==0 and itCnt<nIter:
            meanOld, stdOld, madOld = ms['mean'], ms['stdev'], ms['madfm']
            minVal = ms['mean'] - (clip * ms['madfm'])
            maxVal = ms['mean'] + (clip * ms['madfm'])
            
            # Blank values outside the 3-sigma range
            dataMsk = np.where(np.greater(data, maxVal), np.nan, data)
            dataMsk = np.where(np.less(data, minVal), np.nan, dataMsk)
            
            # Measure the statistics
            ms = calc_stats(dataMsk)
            dataMsk = []
    
            if ms['mean'] == meanOld and ms['madfm'] == madOld:
                convergFlg = 1
            itCnt += 1

    return ms


#-----------------------------------------------------------------------------#
def parse_uvindex(fileName):
    """
    Parse the UVINDEX log file for frequency, source and time information.
    """
    
    nChan = []
    sFreq = []
    dFreq = []
    sources = []
    freqSrch = False
    sourceSrch =False
    FH = open(fileName, "r")
    
    # Regular expressions to match the freq and source tables
    freqMarkerRe = re.compile('Frequency Configuration \d+')
    freqEntryRe = re.compile('^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+GHz')
    srcMarkerRe = re.compile('.*following pointings:')
    srcEntryRe = re.compile('(\S+)\s+\S{0,1}\s+\d+:\d{2}:\d{2}')
    
    # Read in the input file, line by line
    for line in FH:
        line = line.rstrip("\n\r")

        # Search for the frequency settups
        if freqMarkerRe.match(line):
            freqSrch=True
        if freqSrch:
            mch = freqEntryRe.match(line)
            if mch:
                freqSrch=False
                nChan.append(float(mch.group(1)))
                sFreq.append(float(mch.group(2)))
                dFreq.append(float(mch.group(3)))

        # Search for the source names
        if srcMarkerRe.match(line):
            sourceSrch=True
        if sourceSrch:
            mch = srcEntryRe.match(line)
            if mch:
                sources.append(mch.group(1))
            #elif not mch and len(sources)>1:
            #    source_srch=False

    FH.close()

    # Consolidate the frequency information
    freqParams=[]
    for i in range(len(nChan)):
        minMaxMid = []
        
        # Work out the mid and end frequencies
        minMaxMid.append(sFreq[i])
        minMaxMid.append(sFreq[i]+(dFreq[i]*(nChan[i]-1)))
        minMaxMid.append(sFreq[i]+(dFreq[i]*(nChan[i]-1)/2))
        freqParams.append([minMaxMid, dFreq[i], nChan[i]])
    
    return [freqParams, sources]


#-----------------------------------------------------------------------------#
def parse_uvindex_new(fileName):
    """
    Parse the UVINDEX log file for frequency, source and time information.
    """

    # Dictionaries of lists to hold the tables
    timeDol = {'timeStamp':[],
               'pntName':[],
               'calCode':[],
               'nAnts':[],
               'freqNum':[],
               'fieldName':[],
               'RA_deg':[],
               'Dec_deg':[]}
    freqDol = {'nChan':[],
               'sFreq':[],
               'dFreq':[],
               'minFreq':[],
               'midFreq':[],
               'maxFreq':[]}
    pntDol = {'pntName':[],
              'calCode':[],
              'RA_hms':[],
              'Dec_dms':[],
              'RA_deg':[],
              'Dec_deg':[],
              'fieldName':[]}
        
    # Regular expressions to match the start of the freq, source and time
    # tables, and the entries in each table.
    timeMarkerRe = re.compile('.*Config   No.')
    timeEntryRe = re.compile('^(\d{2}\w{3}\d{2}:\S+)\s+(\S+)\s+(\w)\s' +
                             '{2}(\d)\s+\d+\s+\d+\s+(\d)')
    freqMarkerRe = re.compile('Frequency Configuration \d+')
    freqEntryRe = re.compile('^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+GHz')
    srcMarkerRe = re.compile('.*following pointings:')
    srcEntryRe = re.compile('(\S+)\s+(\w)\s+(\d{2}:\d{2}:\d{2}\.\d+)\s+' +
                            '([+|-]\d{2}:\d{2}:\d{2}\.\d+)\s+')
    
    # Read in the input file, line by line
    doTimeSrch = False
    doFreqSrch = False
    doSrcSrch = False
    FH = open(fileName, "r")
    for line in FH:
        line = line.rstrip("\n\r")

        # Search for the time entries
        if timeMarkerRe.match(line):
            doTimeSrch = True
            doFreqSrch = False
            doSrcSrch = False
        if doTimeSrch:
            mch = timeEntryRe.match(line)
            if mch:
                timeDol['timeStamp'].append(mch.group(1))
                timeDol['pntName'].append(mch.group(2))
                timeDol['calCode'].append(mch.group(3))
                timeDol['nAnts'].append(mch.group(4))
                timeDol['freqNum'].append(mch.group(5))
                
        # Search for the frequency setups
        if freqMarkerRe.match(line):
            doTimeSrch = False
            doFreqSrch = True
            odSrcSrch = False
        if doFreqSrch:
            mch = freqEntryRe.match(line)
            if mch:
                freqDol['nChan'].append(float(mch.group(1)))
                freqDol['sFreq'].append(float(mch.group(2)))
                freqDol['dFreq'].append(float(mch.group(3)))

        # Search for the source names
        if srcMarkerRe.match(line):
            doTimeSrch = False
            doFreqSrch = False
            doSrcSrch = True
        if doSrcSrch:
            mch = srcEntryRe.match(line)
            if mch:
                pntDol['pntName'].append(mch.group(1))
                pntDol['calCode'].append(mch.group(2))
                pntDol['RA_hms'].append(mch.group(3))
                pntDol['Dec_dms'].append(mch.group(4))    
    FH.close()

    # Create a field table keyed with a unique field name and a
    # pointing table keyed with a unique pointing name
    # Duplicate keys are automatically replaced a dictionary
    fieldRA_hms = {}
    fieldDec_dms = {}
    fieldRA_deg = {}
    fieldDec_deg = {}
    fieldCalcode = {}
    fieldPosNames = {}
    pntFieldNames = {}
    pntRA_deg = {}
    pntDec_deg = {}
    for i in range(len(pntDol['pntName'])):
        fieldName = pntDol['RA_hms'][i].replace(':', '') + \
                    pntDol['Dec_dms'][i].replace(':', '')
        fieldRA_hms[fieldName] = pntDol['RA_hms'][i]
        fieldDec_dms[fieldName] = pntDol['Dec_dms'][i]
        RA_deg = dms2deg(pntDol['RA_hms'][i]) * 15.0
        Dec_deg = dms2deg(pntDol['Dec_dms'][i])
        fieldRA_deg[fieldName] = RA_deg
        fieldDec_deg[fieldName] = Dec_deg
        fieldCalcode[fieldName] = pntDol['calCode'][i]
        if fieldName in fieldPosNames:
            fieldPosNames[fieldName].append(pntDol['pntName'][i])
        else:
            fieldPosNames[fieldName] = [pntDol['pntName'][i]]
        pntFieldNames[pntDol['pntName'][i]] = fieldName
        pntRA_deg[pntDol['pntName'][i]] = RA_deg
        pntDec_deg[pntDol['pntName'][i]] = Dec_deg

        # Fill the name and position columns in the pointing table
        pntDol['fieldName'].append(fieldName)
        pntDol['RA_deg'].append(fieldRA_deg[fieldName])
        pntDol['Dec_deg'].append(fieldDec_deg[fieldName])

    # Calculate the frequency information
    for i in range(len(freqDol['nChan'])):
        minMidMax = []
        minMidMax.append(freqDol['sFreq'][i])
        minMidMax.append(freqDol['sFreq'][i] +
                         (freqDol['dFreq'][i]*(freqDol['nChan'][i]-1)/2))
        minMidMax.append(freqDol['sFreq'][i] +
                         (freqDol['dFreq'][i]*(freqDol['nChan'][i]-1)))
        minMidMax.sort()
        freqDol['minFreq'].append(minMidMax[0])
        freqDol['midFreq'].append(minMidMax[1])
        freqDol['maxFreq'].append(minMidMax[2])

    # Poplulate the fieldName, RA_deg and Dec_deg fields in the time table
    for i in range(len(timeDol['timeStamp'])):
        timeDol['fieldName'].append(pntFieldNames[timeDol['pntName'][i]])
        timeDol['RA_deg'].append(pntRA_deg[timeDol['pntName'][i]])
        timeDol['Dec_deg'].append(pntDec_deg[timeDol['pntName'][i]])

    # Write the information in each table to numpy record arrays
    dayNameLst = ['0'*20] * len(timeDol['timeStamp'])
    timeTab = np.rec.fromrecords(zip(timeDol['timeStamp'],
                                     dayNameLst,
                                     timeDol['pntName'],
                                     timeDol['calCode'],
                                     timeDol['nAnts'],
                                     timeDol['freqNum'],
                                     timeDol['fieldName'],
                                     timeDol['RA_deg'],
                                     timeDol['Dec_deg']),
                                 names=['timeStamp',
                                        'dayName',
                                        'pntName',
                                        'calCode',
                                        'nAnts',
                                        'freqNum',
                                        'fieldName',
                                        'RA_deg',
                                        'Dec_deg'])
    pntTab = np.rec.fromrecords(zip(pntDol['pntName'],
                                    pntDol['calCode'],
                                    pntDol['RA_hms'],
                                    pntDol['Dec_dms'],
                                    pntDol['RA_deg'],
                                    pntDol['Dec_deg'],
                                    pntDol['fieldName']),
                                names=['pntName',
                                       'calCode',
                                       'RA_hms',
                                       'Dec_dms',
                                       'RA_deg',
                                       'Dec_deg',
                                       'fieldName'])
    freqTab = np.rec.fromrecords(zip(freqDol['nChan'],
                                     freqDol['sFreq'],
                                     freqDol['dFreq'],
                                     freqDol['minFreq'],
                                     freqDol['midFreq'],
                                     freqDol['maxFreq'],),
                                 names=['nChan',
                                        'sFreq',
                                        'dFreq',
                                        'minFreq',
                                        'midFreq',
                                        'maxFreq'])

    # Assemble the field coordinate table
    fieldNameLst = fieldRA_deg.keys()
    fieldNameLst.sort()
    RALst_deg = []
    DecLst_deg = []
    fieldCalcodeLst = []
    for field in fieldNameLst:
        RALst_deg.append(fieldRA_deg[field])
        DecLst_deg.append(fieldDec_deg[field])
        fieldCalcodeLst.append(fieldCalcode[field])
    fieldTab = np.rec.fromrecords(zip(fieldNameLst,
                                      RALst_deg,
                                      DecLst_deg,
                                      fieldCalcodeLst),
                                  names=['fieldName',
                                         'RA_deg',
                                         'Dec_deg',
                                         'calCode'])
                                      
    return timeTab, pntTab, freqTab, fieldTab


#-----------------------------------------------------------------------------#
def timestamp_to_datetime(timeStampStr):
    """
    Convert a MIRIAD timestamp to a datetime object.
    """

    stampRe = re.compile('(\d{2})(\w{3})(\d{2}):(\d{2}):(\d{2}):(\d{2})\.(\d)')
    monthDict = {'JAN': 1,
                 'FEB': 2,
                 'MAR': 3,
                 'APR': 4,
                 'MAY': 5,
                 'JUN': 6,
                 'JUL': 7,
                 'AUG': 8,
                 'SEP': 9,
                 'OCT': 10,
                 'NOV': 11,
                 'DEC': 12}

    mch = stampRe.match(timeStampStr)
    if not mch:
        return None
    else:

        # MIRIAD assumes the old century turnover is at 2040
        year = int(mch.group(1))
        if year>40:
            year +=  1900
        else:
            year +=  2000
        month = monthDict[mch.group(2)]
        day = int(mch.group(3))
        hour = int(mch.group(4))
        minute = int(mch.group(5))
        sec = int(mch.group(6))
        uSec = int(mch.group(7)) * 100000

    return datetime.datetime(year, month, day, hour, minute, sec, uSec)


#-----------------------------------------------------------------------------#
def datetime_to_timestamp(dtObj):
    """
    Convert a datetime object to a MIRIAD timestamp.
    """

    monthDict = {1: 'JAN',
                 2: 'FEB',
                 3: 'MAR',
                 4: 'APR',
                 5: 'MAY',
                 6: 'JUN',
                 7: 'JUL',
                 8: 'AUG',
                 9: 'SEP',
                 10: 'OCT',
                 11: 'NOV',
                 12: 'DEC'}
    
    year = str(dtObj.year)[2:]
    month = monthDict[dtObj.month]
    day = "%02d" % dtObj.day
    hour = ":%02d" % dtObj.hour
    minute = ":%02d" % dtObj.minute
    sec = ":%02d" % dtObj.second
    uSec = ".%01d" % int(dtObj.microsecond/1e5)

    return year + month + day + hour + minute + sec + uSec
        
    
#-----------------------------------------------------------------------------#
def query_obs_in_field(cursor, fieldName, splitFlag=True, useFlag=True,
                       dayName=None):
    """
    Query the observations (i.e., timestamps) corresponding to each field.
    By default only return those timestamps which have been successfully
    split out (observations.splitFlag=1 in the DB) and are flagged to be
    used in the imaging (observations.useFlag=1). By default both these flags
    are set to 1 by the splitting script.
    """
    
    # Query the uv-data contributing to each field
    sql = "SELECT * FROM observations WHERE fieldName = ?"
    if splitFlag:
        sql += " AND splitFlag = 1"
    if useFlag:
        sql += " AND useFlag = 1 "
    args = [fieldName]
    if dayName is not None:
        if len(dayName)>0:
            sql += " AND dayname REGEXP ? "
            args = [fieldName, dayName]

    # Fetch the rows
    obsTab = select_into_arr(cursor, sql, args=(args))
    
    return obsTab

#-----------------------------------------------------------------------------#
def query_obs_in_pointing(cursor, pntName, splitFlag=True, useFlag=True,
                       dayName=None):
    """
    Query the split mosaics uv-data corresponding to each pointing (unique
    mosaic field coordinates given given by pntName).
    """
    
    # Query the uv-data contributing to each field
    sql = "SELECT * FROM mosaic_obs WHERE pntName = ?"
    if splitFlag:
        sql += " AND splitFlag = 1"
    if useFlag:
        sql += " AND useFlag = 1 "
    args = [pntName]
    if dayName is not None:
        if len(dayName)>0:
            sql += " AND dayname REGEXP ? "
            args = [pntName, dayName]

    # Fetch the rows
    obsTab = select_into_arr(cursor, sql, args=(args))
    
    return obsTab


#-----------------------------------------------------------------------------#
def query_tile_parameters(cursor, tileID):
    """
    Query the parameters of a tile given the ID number.
    """
    
    sql = "SELECT * FROM tile_coords WHERE tileID = ?"
    tileTab = select_into_arr(cursor, sql, args=([tileID]))
    
    return tileTab


#-----------------------------------------------------------------------------#
def query_fields_in_tile(cursor, tileID, radFld_deg):    
    """
    Query the fields in a tile given the tileID and the radius-of-view.
    """
    
    # Dictionary to store found field params
    fldParmDict = dict()
    
    # Query the tile coordinates and parameters
    # ('tileID', 'RA_deg', 'Dec_deg', 'l_deg', 'b_deg',
    #  'pixScaleX_asec', 'pixscaleY_asec', 'nPixX', 'nPixY')
    tlTab = query_tile_parameters(cursor, tileID)[0]
    
    # Calculate the high and low tile Dec boundaries
    dDec_deg = abs(tlTab['nPixY'] * tlTab['pixscaleY_asec'] / 3600.0)
    decHigh_deg = tlTab['Dec_deg'] + dDec_deg / 2.0
    decLow_deg = tlTab['Dec_deg'] - dDec_deg / 2.0

    # First search for fields within the DEC range
    # Inclusive boundary is a further radFld_deg outside tile boundary
    sql = """
    SELECT
    field_coords.fieldName,
    field_coords.RA_deg,
    field_coords.Dec_deg
    FROM field_coords
    WHERE field_coords.Dec_deg < ? + (?)
    AND field_coords.Dec_deg > ? - (?)
    AND field_coords.calcode = 'n'
    """
    fldTab = select_into_arr(cursor, sql, args=([decHigh_deg, radFld_deg,
                                                 decLow_deg, radFld_deg]))
    
    if len(fldTab)>0:
            
        # Loop through the fields if they overlap (in RA) with the tile
        #-----------------------------------------------------------------#
        # The blocky procedure below is needed to properly                #
        # account for the cos(Dec) factor in the coordinates. We generate #
        # the RA coordinates of the current field, and the RA coordinates #
        # of the tile border at the upper and lower Dec of the current    #
        # field. If any of the field corners is inside the tile, then the #
        # field is accepted for mosaicing.                                #
        #-----------------------------------------------------------------#
        for e in fldTab:
        
            # Calculate the coordinates at the field corners
            fldDecTop_deg = e['Dec_deg'] + radFld_deg
            fldDecBot_deg = e['Dec_deg'] - radFld_deg
            dRAtop_deg = radFld_deg / (m.cos(m.radians(fldDecTop_deg)))
            dRAbot_deg = radFld_deg / (m.cos(m.radians(fldDecBot_deg)))
            fldRAtopL_deg = e['RA_deg'] + dRAtop_deg # RA, top-left    
            fldRAtopR_deg = e['RA_deg'] - dRAtop_deg # RA, top-right   
            fldRAbotL_deg = e['RA_deg'] + dRAbot_deg # RA, bottom-left 
            fldRAbotR_deg = e['RA_deg'] - dRAbot_deg # RA, bottom-right

            # Calculate the left and right RA boundaries of the *tile*
            # at the top and bottom DEC of the *current field*
            dRA_deg = tlTab['nPixX'] * tlTab['pixScaleX_asec'] /3600.0/2.0
            dRAtop_deg = dRA_deg/(m.cos(m.radians(fldDecTop_deg)))
            dRAbot_deg = dRA_deg/(m.cos(m.radians(fldDecBot_deg)))
            tlRAtopL_deg = tlTab['RA_deg'] + dRAtop_deg # RA, top-left
            tlRAtopR_deg = tlTab['RA_deg'] - dRAtop_deg # RA, top-right
            tlRAbotL_deg = tlTab['RA_deg'] + dRAbot_deg # RA, bottom-left
            tlRAbotR_deg = tlTab['RA_deg'] - dRAbot_deg # RA, bottom-right

            # Measure the field against the tile
            trTest = fldRAtopR_deg <= tlRAtopL_deg # top-right is low enough
            tlTest = fldRAtopL_deg >= tlRAtopR_deg # top-left is high enough
            brTest = fldRAbotR_deg <= tlRAbotL_deg # bottom-right is low enough
            blTest = fldRAbotL_deg >= tlRAbotR_deg # bottom-left is high enough
                                  
            # Include the field in the list if any corner overlaps the tile
            if (tlTest and trTest) or (blTest and brTest):
                fldParmDict[e['fieldName']] = [e['RA_deg'], e['Dec_deg']]

    # Return a dictionary of field parameters {name: [ra, dec]}
    return fldParmDict


#-----------------------------------------------------------------------------#
def query_pointings_in_tile(cursor, tileID, radFld_deg):    
    """
    Query the mosaic pointings in a tile given the tileID and the
    radius-of-view.
    """
    
    # Dictionary to store found field params
    pntParmDict = dict()
    
    # Query the tile coordinates and parameters
    # ('tileID', 'RA_deg', 'Dec_deg', 'l_deg', 'b_deg',
    #  'pixScaleX_asec', 'pixscaleY_asec', 'nPixX', 'nPixY')
    tlTab = query_tile_parameters(cursor, tileID)[0]
    
    # Calculate the high and low tile Dec boundaries
    dDec_deg = abs(tlTab['nPixY'] * tlTab['pixscaleY_asec'] / 3600.0)
    decHigh_deg = tlTab['Dec_deg'] + dDec_deg / 2.0
    decLow_deg = tlTab['Dec_deg'] - dDec_deg / 2.0
    
    # First search for fields within the DEC range
    # Inclusive boundary is a further radFld_deg outside tile boundary
    sql = """
    SELECT
    pointings.pntName,
    pointings.RA_deg,
    pointings.Dec_deg
    FROM pointings
    WHERE pointings.Dec_deg < ? + (?)
    AND pointings.Dec_deg > ? - (?)
    AND pointings.calcode = 'n'
    """
    pntTab = select_into_arr(cursor, sql, args=([decHigh_deg, radFld_deg,
                                                 decLow_deg, radFld_deg]))
    
    if len(pntTab)>0:
            
        # Loop through the pointings if they overlap (in RA) with the tile
        #-----------------------------------------------------------------#
        # The blocky procedure below is needed to properly                #
        # account for the cos(Dec) factor in the coordinates. We generate #
        # the RA coordinates of the current field, and the RA coordinates #
        # of the tile border at the upper and lower Dec of the current    #
        # field. If any of the field corners is inside the tile, then the #
        # field is accepted for mosaicing.                                #
        #-----------------------------------------------------------------#
        for e in pntTab:
        
            # Calculate the coordinates at the field corners
            pntDecTop_deg = e['Dec_deg'] + radFld_deg
            pntDecBot_deg = e['Dec_deg'] - radFld_deg
            dRAtop_deg = radFld_deg / (m.cos(m.radians(pntDecTop_deg)))
            dRAbot_deg = radFld_deg / (m.cos(m.radians(pntDecBot_deg)))
            pntRAtopL_deg = e['RA_deg'] + dRAtop_deg # RA, top-left    
            pntRAtopR_deg = e['RA_deg'] - dRAtop_deg # RA, top-right   
            pntRAbotL_deg = e['RA_deg'] + dRAbot_deg # RA, bottom-left 
            pntRAbotR_deg = e['RA_deg'] - dRAbot_deg # RA, bottom-right

            # Calculate the left and right RA boundaries of the *tile*
            # at the top and bottom DEC of the *current field*
            dRA_deg = tlTab['nPixX'] * tlTab['pixScaleX_asec'] /3600.0/2.0
            dRAtop_deg = dRA_deg/(m.cos(m.radians(pntDecTop_deg)))
            dRAbot_deg = dRA_deg/(m.cos(m.radians(pntDecBot_deg)))
            tlRAtopL_deg = tlTab['RA_deg'] + dRAtop_deg # RA, top-left
            tlRAtopR_deg = tlTab['RA_deg'] - dRAtop_deg # RA, top-right
            tlRAbotL_deg = tlTab['RA_deg'] + dRAbot_deg # RA, bottom-left
            tlRAbotR_deg = tlTab['RA_deg'] - dRAbot_deg # RA, bottom-right

            # Measure the field against the tile
            trTest = pntRAtopR_deg <= tlRAtopL_deg # top-right is low enough
            tlTest = pntRAtopL_deg >= tlRAtopR_deg # top-left is high enough
            brTest = pntRAbotR_deg <= tlRAbotL_deg # bottom-right is low enough
            blTest = pntRAbotL_deg >= tlRAbotR_deg # bottom-left is high enough
                                  
            # Include the field in the list if any corner overlaps the tile
            if (tlTest and trTest) or (blTest and brTest):
                pntParmDict[e['pntName']] = [e['RA_deg'], e['Dec_deg']]

    # Return a dictionary of field parameters {name: [ra, dec]}
    return pntParmDict


#-----------------------------------------------------------------------------#
def parse_uvlist_spectral(fileName):
    """
    Parse the UVLIST log file for channel information. Must have been run
    with option=spectral.
    """
    
    outDict = {}
    FH = open(fileName, "r")
    
    # Regular expressions to match velocity tables
    restFreqRe = re.compile('^Rest frequency\s+\(GHz\) :\s+(.*)')
    startChanRe = re.compile('^Start channel\s+:\s+(\d+)')
    nChanRe = re.compile('^Number of channels\s+:\s+(.*)')
    startFreqRe = re.compile('^Start frequency\s+\(GHz\) :\s+(.*)')
    freqIncRe =  re.compile('^Frequency increment \(GHz\) :\s+(.*)')
    
    radioVelRe = re.compile('^Radio Velocities:')
    opticalVelRe = re.compile('^Optical Velocities:')
    startVelRe = re.compile('\s+Start velocity\s+\(km\/s\):\s+(.*)')
    endVelRe = re.compile('\s+End velocity\s+\(km\/s\):\s+(.*)')
    velIncRe = re.compile('\s+Velocity increment\(km\/s\):\s+(.*)')
    
    # Read in the input file, line by line
    dosrch = False
    for line in FH:
        line = line.rstrip("\n\r")

        # Match the frequency parameters and number of channels
        mch = restFreqRe.match(line)
        if mch:
            outDict['restFreq_GHz'] = float(mch.group(1))
        mch = startChanRe.match(line)
        if mch:
            outDict['startChan'] = float(mch.group(1))
        mch = nChanRe.match(line)
        if mch:
            outDict['nChan'] = int(mch.group(1))
        mch = startFreqRe.match(line)
        if mch:
            outDict['startFreq_GHz'] = float(mch.group(1))
        mch = freqIncRe.match(line)
        if mch:
            outDict['freqInc_GHz'] = float(mch.group(1))

        # Toggle on/off velocity matching
        if radioVelRe.match(line):
            dosrch = True
        if opticalVelRe.match(line):
            dosrch = False

        # Match the radio velocities
        if dosrch:
            mch = startVelRe.match(line)
            if mch:
                outDict['startVel_kms'] = float(mch.group(1))
            mch = endVelRe.match(line)
            if mch:
                outDict['endVel_kms'] = float(mch.group(1))
            mch = velIncRe.match(line)
            if mch:
                outDict['velInc_kms'] = float(mch.group(1))

    FH.close()
    
    return outDict


#-----------------------------------------------------------------------------#
def parse_uvlist_array(fileName, flagAntLst=[]):
    """
    Parse the UVLIST log file for baseline lengths. Must have been run
    with option=array.    
    """

    if not isinstance(flagAntLst, list):
        flagAntLst = [flagAntLst]

    
    outDict = {}
    FH = open(fileName, "r")
    
    # Regular expressions to match velocity tables
    antposRe = re.compile('^\s+(\d)\s+(\S+)\s+(\S+)\s+(\S+)')
    
    # Read in the input file, line by line
    posDict = {}
    for line in FH:
        line = line.rstrip("\n\r")

        # Match the antenna positions
        mch = antposRe.match(line)
        if mch:
            if not (int(mch.group(1)) in flagAntLst):
                posDict[int(mch.group(1))] = np.array([float(mch.group(2)),
                                                       float(mch.group(3)),
                                                       float(mch.group(4))])

    # Determin the baseline lengths
    baselineDict = {}
    ants = posDict.keys()
    ants.sort()
    for i in range(len(ants)):
        ant1 = ants[i]
        for k in range(i+1, len(ants)):
            ant2 = ants[k]
            length_m = np.sum((posDict[ant1]-posDict[ant2])**2.0)**0.5
            baselineDict[length_m] = "%d,%d" % (ant1, ant2)
    FH.close()

    return baselineDict

    
#-----------------------------------------------------------------------------#
def parse_imstat(fileName):
    """
    Parse the IMSTAT log file for information on the image properties.
    """

    outDict = {}
    FH = open(fileName, "r")

    # Regular expressions 
    totalRe = re.compile('^\sTotal')

    # Read in the input file, line by line
    dosrch = False
    for line in FH:
        line = line.rstrip("\n\r")

        # Toggle on/off on 'Total" line
        if totalRe.match(line):            
            dosrch = True
            continue
        
        # Split up the last line (each field is 10 characters)
        #      Sum      Mean      rms     Maximum   Minimum    Npoints
        #    0.100     2.316E-08 1.728E-04 9.142E-04-8.482E-04 4318084
        if dosrch:
            dosrch = False
            outDict['sum'] = float(line[21:31])
            outDict['mean'] = float(line[31:41])
            outDict['rms'] = float(line[41:51])
            outDict['max'] = float(line[51:61])
            outDict['min'] = float(line[61:71])
            outDict['npix'] = int(line[71:81])
            
    return outDict


#-----------------------------------------------------------------------------#
def make_lin_axis(naxis, crpix, crval, cdelt):
    """
    Create an array containing the axis values, assuming a simple linear
    projection scheme. A
    """
    
    start = crval + (1 - crpix) * cdelt
    stop = (crval + (naxis + 1 - crpix) *  cdelt)
    
    return np.arange(start, stop, cdelt)


#-----------------------------------------------------------------------------#
def sort_nicely(l):
    """
    Sort a list in the order a human would.
    """
    
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    l.sort( key=alphanum_key ) 


#-----------------------------------------------------------------------------#
def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


#-----------------------------------------------------------------------------#
def pause():
    print "\nPress <RETURN> to continue ...",
    raw_input()
    
