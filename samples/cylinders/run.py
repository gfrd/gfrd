#!/usr/bin/env python

from egfrd import *
from vtklogger import VTKLogger

def run( ):
    s = EGFRDSimulator('run')
    #s = EGFRDSimulator()
    factor = 1e9
    L = factor * 1e-6
    s.setWorldSize( L )

    D = factor * factor * 1e-12
    sigma = factor * 2.5e-8
    O = Species( 'O', 0, sigma )
    s.addSpecies( O )
    A = Species( 'A', D, sigma )
    s.addSpecies( A )
    B = Species( 'B', D, sigma )
    s.addSpecies( B )
    C = Species( 'C', D, sigma )
    s.addSpecies( C )

    r1 = BindingReactionType( A, B, C, 1e18 )
    s.addReactionType( r1 )
    r2 = UnbindingReactionType( C, A, B, 1e0 )
    s.addReactionType( r2 )

    box1 = CuboidalSurface( [0,0,0], [L,L,L], 'world' )
    dna = CylindricalSurface( [L/2,L/2,L/2], 2*sigma, [0,1,0], L/2, 'dna')
    s.addSurface( dna )
    membrane = PlanarSurface( [L/2,L/2,0], [1,0,0], [0,1,0], L/2, L/2, sigma, 'membrane', )
    s.addSurface( membrane )

    # Test periodic boundary condition with pair.
    s.placeParticle( A, [500,500,975], membrane )
    s.placeParticle( B, [500,500,25.1], membrane )





    #s.throwInParticles( A, 2, box1 )
    #s.throwInParticles( A, 2, dna )
    #s.throwInParticles( B, 2, dna )

    #s.throwInParticles( C, 2, dna )
    #s.throwInParticles( A, 2, membrane )
    #s.throwInParticles( B, 2, membrane )

    #s.throwInParticles( C, 2, membrane )
    #s.placeParticle( O, [L/2,L/2,L/2], dna )





    s.placeParticle( C, [709,748,3], membrane )

    # Initialize before s.cylinderMatrix.add.
    s.initialize()

    #s.cylinderMatrix.add( (dna, 0), dna.outside )

    for i in range(20):
        s.step()

    s.stop( s.t )
    
if __name__ == '__main__':
    run( )
