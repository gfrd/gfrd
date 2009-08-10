#!/usr/bin/env python

from egfrd import *
from vtklogger import VTKLogger

def run( ):
    '''
    Simulator.
    '''
    s = EGFRDSimulator('run')
    #s = EGFRDSimulator()

    factor = 1e9
    L = factor * 1e-6
    s.setWorldSize( L )

    D = factor * factor * 1e-12
    sigma = factor * 2.5e-8


    '''
    Species.

    Todo: somehow make sure a species can only be on 1 surface.
    '''
    O = Species( 'O', 0, sigma )
    s.addSpecies( O )

    membraneA = Species( 'membraneA', D, sigma )
    s.addSpecies( membraneA )
    membraneB = Species( 'membraneB', D, sigma )
    s.addSpecies( membraneB )
    membraneC = Species( 'membraneC', D, sigma )
    s.addSpecies( membraneC )

    membrane2A = Species( 'membrane2A', D, sigma )
    s.addSpecies( membrane2A )
    membrane2B = Species( 'membrane2B', D, sigma )
    s.addSpecies( membrane2B )
    membrane2C = Species( 'membrane2C', D, sigma )
    s.addSpecies( membrane2C )

    dnaA = Species( 'dnaA', D, sigma )
    s.addSpecies( dnaA )
    dnaB = Species( 'dnaB', D, sigma )
    s.addSpecies( dnaB )
    dnaC = Species( 'dnaC', D, sigma )
    s.addSpecies( dnaC )

    worldA = Species( 'worldA', D, sigma )
    s.addSpecies( worldA )
    worldB = Species( 'worldB', D, sigma )
    s.addSpecies( worldB )
    worldC = Species( 'worldC', D, sigma )
    s.addSpecies( worldC )


    '''
    Surfaces.
    '''
    box1 = CuboidalSurface( [0,0,0], [L,L,L], 'world' )
    s.cuboidalSurfaces = [ box1 ]

    dna = CylindricalSurface( [L/2,L/2,L/2], 2*sigma, [0,1,0], L/2, 'dna')
    #s.addSurface( dna )

    membrane = PlanarSurface( [L/2,L/2,2*L/10], [1,0,0], [0,1,0], L/2, L/2, sigma, 'membrane', )
    s.addSurface( membrane )

    membrane2 = PlanarSurface( [L/2,L/2,8*L/10], [1,0,0], [0,1,0], L/2, L/2, sigma, 'membrane2', )
    s.addSurface( membrane2 )


    '''
    Reactions.
    '''
    '''
    r1 = BindingReactionType( A, B, C, 1e18 )
    s.addReactionType( r1 )
    r2 = UnbindingReactionType( C, A, B, 1e0 )
    s.addReactionType( r2 )
    '''

    i1 = SurfaceBindingInteractionType( worldA, dna, dnaA, 1e8 )
    s.addInteractionType( i1 )
    r1 = SurfaceUnbindingReactionType( dnaA, worldA, 1e10 )
    s.addReactionType( r1 )

    i2 = SurfaceBindingInteractionType( worldA, membrane, membraneA, 1 )
    s.addInteractionType( i2 )
    r2 = SurfaceUnbindingReactionType( membraneA, worldA, 1e10 )
    s.addReactionType( r2 )

    i2 = SurfaceBindingInteractionType( worldA, membrane2, membrane2A, 1 )
    s.addInteractionType( i2 )
    r2 = SurfaceUnbindingReactionType( membrane2A, worldA, 1e10 )
    s.addReactionType( r2 )


    '''
    Particles.
    '''
    #s.throwInParticles( worldA, 5, box1 )
    #s.throwInParticles( worldB, 1, box1 )
    #s.throwInParticles( worldC, 1, box1 )

    #s.throwInParticles( dnaA, 5, dna )
    #s.throwInParticles( dnaB, 1, dna )
    #s.throwInParticles( dnaC, 1, dna )

    s.throwInParticles( membraneA, 5, membrane )
    #s.throwInParticles( membraneB, 1, membrane )
    #s.throwInParticles( membraneC, 1, membrane )

    s.throwInParticles( membrane2A, 5, membrane2 )
    #s.throwInParticles( membrane2B, 1, membrane2 )
    #s.throwInParticles( membrane2C, 1, membrane2 )

    #s.placeParticle( O, [L/2,L/2,L/2], dna )


    '''
    Simulation.
    '''
    s.initialize()
    for i in range(200):
        try:
            s.step()
        except Stop, message:
            print message
            break
    s.stop( s.t )
    
if __name__ == '__main__':
    run( )
