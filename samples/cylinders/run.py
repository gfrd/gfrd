#!/usr/bin/env python

from egfrd import *
from vtklogger import VTKLogger

def run( ):
    s = EGFRDSimulator('run')
    factor = 1e9
    L = factor * 1e-6
    s.setWorldSize( L )

    D = factor * factor * 1e-12
    sigma = factor * 2.5e-8
    O = Species( 'O', 0, sigma )
    s.addSpecies( O )
    A = Species( 'A', D, sigma )
    s.addSpecies( A )

    dnaR = sigma

    box1 = CuboidalSurface( [0,0,0], [L,L,L], 'world' )

    dna = CylindricalSurface( [L/2,L/2,L/2], sigma, [0,1,0], L/2, 'dna')
    s.addSurface( dna )

    membrane = PlanarSurface( [L/2,L/2,0], [1,0,0], [0,1,0], L/2, L/2, sigma, 'membrane', )
    s.addSurface( membrane )

    #s.throwInParticles( A, 2, box1 )
    s.throwInParticles( A, 2, dna )
    s.throwInParticles( A, 2, membrane )
    #s.placeParticle( O, [L/2,L/2,L/2], dna )
    #s.placeParticle( A, [100,100,0], membrane )
    #s.placeParticle( A, [100,200,0], membrane )

    # Initialize before s.cylinderMatrix.add.
    s.initialize()

    #s.cylinderMatrix.add( (dna, 0), dna.outside )

    for i in range(50):
        s.step()

    s.stop( s.t )
    
if __name__ == '__main__':
    run( )
