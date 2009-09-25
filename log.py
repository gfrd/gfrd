import logging
import logging.handlers
import os
import sys

'''
This logging faciliy can be used for info en debug messages and such.

Usage:
    * from log import *
    * log.info( 'foo' )
    * log.debug( 'bar' )
'''
def setupLogging():
    global log 

    log = logging.getLogger( 'pdpp' )

    if 'LOGFILE' in os.environ:
        if 'LOGSIZE' in os.environ and int( os.environ[ 'LOGSIZE' ] ) != 0:
            handler = logging.handlers.\
                RotatingFileHandler( os.environ[ 'LOGFILE' ], mode='w',
                                     maxBytes=int( os.environ[ 'LOGSIZE' ] ) )
        else:
            handler = logging.FileHandler( os.environ[ 'LOGFILE' ], 'w', )
            
    else:
        handler = logging.StreamHandler( sys.stdout )

    formatter = logging.Formatter( '%(message)s' )
    handler.setFormatter( formatter )
    log.addHandler( handler )
        
    LOGLEVELS = { 'CRITICAL': logging.CRITICAL,
                  'ERROR': logging.ERROR,
                  'WARNING': logging.WARNING,
                  'INFO': logging.INFO, 
                  'DEBUG': logging.DEBUG, 
                  'NOTSET': logging.NOTSET }

    if 'LOGLEVEL' in os.environ:
        log.setLevel( LOGLEVELS[ os.environ[ 'LOGLEVEL' ] ] )
    else:
        log.setLevel( logging.INFO )


setupLogging()


