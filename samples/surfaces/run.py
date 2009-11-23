#!/usr/bin/env python

from egfrd import *
from vtklogger import VTKLogger

def run( ):
    '''
    Simulator.
    '''
    factor = 1 #1e9
    L = factor * 1e-6

    s = EGFRDSimulator( worldSize=L )

    D = factor * factor * 1e-12
    sigma = factor * 2.5e-8


    '''
    Surfaces.
    '''
    dna = CylindricalSurface( [L/2,L/2,L/2], 2*sigma, [0,1,0], L/2, 'dna')
    s.addSurface( dna )

    membrane = PlanarSurface( [L/2,L/2,2*L/10], [1,0,0], [0,1,0], L/2, L/2, sigma, 'membrane', )
    s.addSurface( membrane )

    membrane2 = PlanarSurface( [L/2,L/2,8*L/10], [1,0,0], [0,1,0], L/2, L/2, sigma, 'membrane2', )
    s.addSurface( membrane2 )


    '''
    Species.

    Todo: somehow make sure a species can only be on 1 surface.
    '''
    #O = Species( 'O', 0, sigma )
    #s.addSpecies( O )

    Am = Species( 'A', D, sigma )
    s.addSpecies( Am, membrane )
    Bm = Species( 'B', D, sigma)
    s.addSpecies( Bm, membrane  )
    Cm = Species( 'C', D, sigma )
    s.addSpecies( Cm, membrane  )

    Am2 = Species( 'A', D, sigma )
    s.addSpecies( Am2, membrane2 )
    Bm2 = Species( 'B', D, sigma )
    s.addSpecies( Bm2, membrane2 )
    Cm2 = Species( 'C', D, sigma )
    s.addSpecies( Cm2, membrane2 )

    '''
    Ad = Species( 'A', D, sigma )
    s.addSpecies( Ad, dna )
    Bd = Species( 'B', D, sigma )
    s.addSpecies( Bd, dna )
    Cd = Species( 'C', D, sigma )
    s.addSpecies( Cd, dna )
    '''

    Aw = Species( 'A', D, sigma )
    s.addSpecies( Aw )
    Bw = Species( 'B', D, sigma )
    s.addSpecies( Bw )
    Cw = Species( 'C', D, sigma )
    s.addSpecies( Cw )


    '''
    Reactions.
    '''
    '''
    r1 = BindingReactionType( A, B, C, 1e18 )
    s.addReactionType( r1 )
    r2 = UnbindingReactionType( C, A, B, 1e0 )
    s.addReactionType( r2 )
    '''

    '''
    i1 = SurfaceBindingInteractionType( Aw, Ad, 1e-5 )
    s.addInteractionType( i1 )
    r1 = SurfaceUnbindingReactionType( Ad, Aw, 1e-5 )
    s.addReactionType( r1 )
    '''

    i2 = SurfaceBindingInteractionType( Aw, Am, 1e-5 )
    s.addInteractionType( i2 )
    r2 = SurfaceUnbindingReactionType( Am, Aw, 1e-5 )
    s.addReactionType( r2 )

    i2 = SurfaceBindingInteractionType( Aw, Am2, 1e-5 )
    s.addInteractionType( i2 )
    r2 = SurfaceUnbindingReactionType( Am2, Aw, 1e-5 )
    s.addReactionType( r2 )


    '''
    Particles.
    '''
    s.throwInParticles( Aw, 5 )
    #s.throwInParticles( Bw, 1 )
    #s.throwInParticles( Cw, 1 )

    #s.throwInParticles( Ad, 2 )
    #s.throwInParticles( Bd, 1 )
    #s.throwInParticles( Cd, 1 )

    s.throwInParticles( Am, 5 )
    #s.throwInParticles( Bm, 1 )
    #s.throwInParticles( Cm, 1 )

    #s.throwInParticles( Am2, 5 )
    #s.throwInParticles( Bm2, 1 )
    #s.throwInParticles( Cm2, 1 )

    #s.placeParticle( O, [L/2,L/2,L/2], dna )


    '''
    Simulation.
    '''
    s.initialize()
    vtklogger = VTKLogger( s, 'run' )
    for i in range(100):
        try:
            vtklogger.log()
            s.step()
        except Exception, message:
            print message
            break
    vtklogger.stop()
    s.stop( s.t )
    
if __name__ == '__main__':
    run( )
