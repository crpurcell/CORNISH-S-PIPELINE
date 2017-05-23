#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_tables.py                                                    #
#                                                                             #
# PURPOSE:  Helper functions and classes for manipulating numpy recarrays     #
#           and sqlite3 databases.                                            #
#                                                                             #
# MODIFIED: 23-Mar-2017 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  register_sqlite3_numpy_dtypes ... setup sqlite3 to read numpy arrays       #
#  regexp                   ... regualr expression function for sqlite3       #
#  schema_to_dtypes         ... parse the schema file for numpy dtypes        #
#  columns_view             ... return a recarray view filtered by column     #
#  insert_arr_db            ... insert all recarray entries into the database #
#  update_arr_db            ... update DB entries using a recarray            #
#  select_into_arr          ... run a SQL query and return a numpy recarray   #
#  select_into_lod          ... run a SQL query and return a list-of-dicts    #
#  rowlst_to_rowdict        ... convert a LOD to a dictionary of dictionaries #
#                                                                             #
#=============================================================================#
import os
import sys
import time
import re
import numpy as np
import sqlite3


#-----------------------------------------------------------------------------#
def register_sqlite3_numpy_dtypes():
    """
    Map numpy data-types to the limited sqlite data-types. This must be called
    before using the sqlite3 database or INSERT statements will fail.
    """
    
    for t in (np.int8, np.int16, np.int32, np.int64,
          np.uint8, np.uint16, np.uint32, np.uint64):
        sqlite3.register_adapter(t, long)
    for t in (np.float16, np.float32, np.float64,
          np.float128, np.double):
        sqlite3.register_adapter(t, float)

        
#-----------------------------------------------------------------------------#
def regexp(expr, item):
    reg = re.compile(expr)    
    return reg.search(item) is not None

    
#-----------------------------------------------------------------------------#
def schema_to_dtypes(schemaFile):
    """
    Parse the SQL file containing the CREATE TABLE definitions and convert to
    a dtype format that python can understand. This can be used to initialise
    numpy record arrays for storing the tables in memory. Assumes the
    create statement has the form:
    
    'CREATE TABLE myTableName (entry1 double [args+], entry2 int(2), ... );'

    The 'CREATE TABLE' statement must have the same case. WARNING: recognises
    only , 'DOUBLE', 'INTEGER' and 'VARCHAR(n)' types in the SQL statement.
    """

    # Return these
    tableNameLst = []
    tableDefDict = {}
    tableSQLdict = {}
    
    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    commaAndSpaces = re.compile(',\s+')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    createRe =  re.compile('^(CREATE TABLE|create table) (\w+)\s*\((.+)\)\s*$')
    
    # Translation function for data types: SQL to recarray
    # RecArray dtypes:
    #   bytes                = b<n>, e.g. b1
    #   ints                 = i<n>, e.g. i2, i4, i8,
    #   unsigned ints        = u<n>, e.g. u1, u2, u4, u8
    #   floats               = f<n>, e.g. f2, f4, f8
    #   complex              = c<n>, e.g. c8, c16
    #   fixed length strings = a<n>, e.g. a10, a100
    # where <n> is the number of bytes / chars, so float32=f4, float64=f8
    def trtype(dtype):
        floatRe = re.compile('^FLOAT')
        doubleRe = re.compile('^DOUBLE')
        intRe = re.compile('^INT')
        charRe = re.compile('^VARCHAR\((\d+)\)')
        if floatRe.match(dtype):
            return 'f4'
        if doubleRe.match(dtype):
            return 'f8'
        if intRe.match(dtype):
            return 'i8'
        mch = charRe.match(dtype)
        if mch:
             return 'a' + mch.group(1)
        return 'f8'    # default to float64
        
    # Loop through the SQL statements
    sqlLst = open(schemaFile).read().split(';')
    for sql in sqlLst:

        # Simplify the SQL statement
        sql = sql.replace('\r', ' ')        # kill carriage-return
        sql = sql.replace('\n', ' ')        # kill newlines
        sql = sql.strip()                   # kill external whitespace
        sql = spaces.sub(' ', sql)          # shrink internal whitespaces
        sql = commaAndSpaces.sub(',', sql)  # kill ambiguous spaces

        mch = createRe.match(sql)
        if mch:
            tableName = mch.group(2)
            tableSQLdict[tableName] = mch.group(0)
            tableNameLst.append(tableName)
            colDefLst = mch.group(3).strip().split(',')
            colDefLst = [x.split(' ')[:2] for x in colDefLst]

            # Translate the data types into a python recarray dtype list
            for i in range(len(colDefLst)):
                colDefLst[i][1] = trtype(colDefLst[i][1])
                colDefLst[i] = tuple(colDefLst[i])

            # Add to the table definition dictionary
            tableDefDict[tableName] = colDefLst

    return tableDefDict, tableSQLdict

        
#-----------------------------------------------------------------------------#
def columns_view(arr, colNameLst=None):
    """
    Return a view of a numpy record array containing only the column names in
    the colNameLst argument. 'colNameLst' should be a list of column names.
    """

    # Default to all column
    if not colNameLst:
        colNameLst = arr.dtype.names
    
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in colNameLst})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)


#-----------------------------------------------------------------------------#
def insert_arr_db(cursor, recArr, tabName, colNameLst=None,
                     insertSQL="INSERT OR REPLACE"):
    """
    Insert a numpy recarray into a database via a cursor object. It is assumed
    that the columns in the recarray and database have the same names and
    compatable datatypes. If the recarray contains columns NOT in the database
    the user must provide a list of existing columns to be to be inserted.
    """
    
    # Default to all cols
    if not colNameLst:
        colNameLst = recArr.dtype.names    
    
    sql = '%s INTO %s (%s) ' % (insertSQL, tabName, ', '.join(colNameLst))
    sql += 'VALUES(%s) ' % (', '.join(['?']*len(colNameLst)))

    #print sql
    #print columns_view(recArr, colNameLst)
    #raw_input()
    cursor.executemany(sql, columns_view(recArr, colNameLst))


#-----------------------------------------------------------------------------#
def update_arr_db(cursor, recArr, tabName, keyName, colNameLst=None):
    """
    Do an UPDATE on existing rows. Unlike the INSERT OR REPLACE function
    above, this only updates the columns given in colNameLst.
    """
    
    # Default to all fields
    if not colNameLst:
        colNameLst = list(recArr.dtype.names)

    # Key must exist in list of field names
    if not keyName in colNameLst:
        print "ERR: Key '%s' not in column list" % keyName
        return

    # Remove the key from the list and format the SQL
    colNameLst.remove(keyName)        
    sql = 'UPDATE %s SET ' % tabName
    sql += '=?, '.join(colNameLst) + '=? '
    sql += 'WHERE %s = ?' % keyName

    # Attach the key to the end of the field list and fetch a view
    # Use fancy indexing to force the key to the last column
    colNameLst.append(keyName)
    a = fields_view(recArr, colNameLst)[colNameLst]
    cursor.executemany(sql, a)

    
#-----------------------------------------------------------------------------#
def select_into_arr(cursor, sql, args=[]):
    """
    Run a SQL query and return a numpy recordarray.
    """

    if args == []:
        cursor.execute(sql)
    else:
        cursor.execute(sql, tuple(args))
    try:
        rows = cursor.fetchall()
        if len(rows) == 0:
            rows = np.array([], dtype='i4')
        else:
            columnNameLst = zip(*cursor.description)[0]
            rows = np.rec.fromrecords(rows, names=columnNameLst)
        return rows
    except Exception:
        print "WARNING: failed to convert SQL result to a recarray!"
        return None


#-----------------------------------------------------------------------------#
def select_into_lod(cursor, sql, args=[]):
    """
    Run a SQL query and return a List of Dictionaries.
    """

    rowLst = []
    if args == []:
        cursor.execute(sql)
    else:
        cursor.execute(sql, tuple(args))
    columnNameLst = zip(*cursor.description)[0]
    rows = cursor.fetchall()
    for row in rows:
        e = {}
        for i in range(len(columnNameLst)):
            e[columnNameLst[i]] = row[i]
        rowLst.append(e)
    return rowLst


#-----------------------------------------------------------------------------#
def rowlst_to_rowdict(rowLst, keyName):
    """
    Convert a list of dictionaries into a dictionary of dictionaries.
    """
    
    rowDict = {}
    for e in rowLst:
        key = e.pop(keyName)
        rowDict[key] = e

    return rowDict
