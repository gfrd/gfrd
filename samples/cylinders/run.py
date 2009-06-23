#!/usr/bin/env python

from egfrd import *
from vtklogger import VTKLogger

def run( ):
    s = EGFRDSimulator()
    L = 1e-6
    s.setWorldSize( L )

    D = 1e-12
    sigma = 2.5e-8
    A = Species( 'A', D, sigma )
    s.addSpecies( A )

    dnaR = sigma


    box1 = CuboidalSurface( [0,0,0], [L,L,L], 'world' )

    dna = CylindricalSurface( [L/2,L/2,L/2], sigma, [0,0,1], L, 'dna')
    s.addSurface( dna )

    s.throwInParticles( A, 1, dna )

    # Initialize before shellMatrix.addCylinder.
    s.initialize()

    #s.shellMatrix.addCylinder( (dna, 0), dna.outside )

    l = VTKLogger(s, 'run1')
    for i in range(100):
        l.log()
        s.step()

    l.stop()
    
if __name__ == '__main__':
    run( )
