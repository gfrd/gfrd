#!/usr/bin/env python

'''
 PYTHONPATH=../.. python -O run.py single.0.out 5e-4 3e-8 100000
'''



from egfrd import *
import numpy

def run( outfilename, T, S, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T, S )
        outfile.write( '%g\n' % d )

        assert t == T

    outfile.close()



def singlerun( T, S ):

    s = EGFRDSimulator( worldSize = 1e-3 )

    s.setUserMaxShellSize( S )

    A = Species( 'A', 1e-12, 5e-9 )
    s.addSpecies( A )
    
    particleA = s.placeParticle( A, [0,0,0] )

    endTime = T
    s.step()

    while 1:
        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()

    distance = s.distance( numpy.array([0,0,0]), particleA.pos )

    return distance, s.t
    
if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), float( sys.argv[3] ),
         int( sys.argv[4] ) )
