#!/usr/bin/env python

#python mapk/plot_mean.py 09/data/mapk3_1e-15_0.25_fixed_0_normal_ALL_tc.dat 09/data/mapk3_1e-15_0.5_fixed_0_normal_ALL_tc.dat 09/data/mapk3_1e-15_1_fixed_0_normal_ALL_tc.dat 09/data/mapk3_1e-15_2_fixed_0_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_0_normal_ALL_tc.dat

# python mapk/plot_mean.py 09/data/mapk3_1e-15_0.25_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_0.5_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_1_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_2_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-2_normal_ALL_tc.dat
# python mapk/plot_mean.py 09/data/mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_0.5_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_1_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_2_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-6_normal_ALL_tc.dat
# python mapk/plot_mean.py 09/data/mapk3_1e-15_4_fixed_0_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-6_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-5_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-4_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-3_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-2_normal_ALL_tc.dat 09/data/mapk3_1e-15_4_fixed_1e-1_normal_ALL_tc.dat


import sys
import os
import glob

import numpy
import scipy.io

from matplotlib.pylab import *

def load_header( filename ):
    file = open( filename )
    header = []
    for line in file.readlines():
        if line[0:2] == '#@':
            hline = line[2:].lstrip()
            header.append( hline )

    return header

def resample( x, y, newx ):

    indices = numpy.searchsorted( x, newx )

    indices = indices.clip( 0, len(y) - 1 )
    #print indices, len(y)

    return y.take( indices )
    

def add_columns( data, ycolumns ):

    y = numpy.array([ data[:,col] for col in ycolumns ]) 

    y = y.sum(0)

    return y


def load_data( filename ):
    ycolumns = [1,]
    #ycolumns = [2,6]
    #ycolumns = [3,5]
    #ycolumns = [2,6,3,5]

    header = load_header( filename )
    print header
    for l in header:
        exec( l )

    #data = numpy.loadtxt( filename )
    data = load( filename )
    x = data[:,0]
    y = add_columns( data, ycolumns )

    return x, y


def plot_file( filename, lp='-' ):

    x, y = load_data( filename )

    #plot_theory( N_K, N_P, Keq, x[-1] )
    plot( x, y, lp )

    #psd( y )
    #ylim( 1, 5e4 )

def plot_mean( filelist, end, l='' ):

    start = 0.
    interval = (end-start) / 1000.
    rx = numpy.mgrid[start:end:interval]

    data = []

    assert filelist

    for filename in filelist:
        print 'file ', filename
        x, y = load_data( filename )
        print x,y
        ry = resample( x, y, rx )
        print ry.shape
        data.append( ry )

        mry = numpy.array( data ).mean( 0 )


    plot( rx, mry, label=l )


def plot_mean_pattern( pattern, end ):
    globpattern = pattern.replace('ALL','*')
    
    l = os.path.basename( os.path.splitext( pattern )[0] )
    print 'pattern ', l

    filelist = glob.glob( globpattern )

    plot_mean( filelist, end )


if __name__ == '__main__':


    import glob
    import os

    xmax = 60

    for pattern in sys.argv[1:]:
        plot_mean_pattern( pattern, xmax )


    #plot_file('/home/shafi/wrk/brown/samples/mapk/Kpp_ODE_0.ecd', 'k-' )

    xticks( size=20 )
    yticks( size=20 )

    xlabel( r'time [s]', size=22 )
    ylabel( r'Kpp', size=22 )



#title( figtitle )

#savefig( 'figs/' + figtitle + '.png', dpi=80 )

    show()
