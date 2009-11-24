#!/usr/bin/env python

from egfrd import *
#from vtklogger import VTKLogger

'''
A + A <-> B
'''

def run( ):
    ''' Dimensions.

    '''
    L = 1e-6            # Size of simulation box.
    D = 1e-12           # Diffusion constant.
    radius = 2.5e-8     # Radius of particles.


    ''' Simulator.

    A cube of dimension (L x L x L) is created, with periodic boundary 
    conditions.

    '''
    s = EGFRDSimulator( worldSize=L )


    ''' Surfaces. 

    Define the surfaces to add to the simulator.

    Usage:
        s.addPlanarSurface( origin, vectorX, vectorY, Lx, Ly, Lz, name )
        s.addCylindricalSurface( origin, radius, orientation, size, name )
    See the docstrings. 

    Note that surface are not allowed to touch or overlap. Surfaces should 
    extend over the whole simulation box, so finite surfaces are not 
    supported.

    Both methods return the surface they create. Assign it 
    to some variable, so it can be used when adding species.

    '''
    WORLD = True
    MEMBRANE = True
    MEMBRANE2 = True
    DNA = True

    if MEMBRANE:
        membrane = s.addPlanarSurface( [L/2,L/2,2*L/10], [1,0,0], [0,1,0],
                                       L/2, L/2, radius, 'membrane' )

    if MEMBRANE2:
        membrane2 = s.addPlanarSurface( [L/2,L/2,8*L/10], [1,0,0], [0,1,0],
                                        L/2, L/2, radius, 'membrane2' )

    if DNA:
        dna = s.addCylindricalSurface( [L/2,L/2,L/2], 2*radius, [0,1,0], L/2,
                                       'dna')


    ''' Species.

    Specify which type of particles can exist on which surface.
    
    Usage:
        s.addSpecies( id, D, radius )
    See the docstring.

    If particles of type 'A' should be able to exist in the world as well as 
    on one of the previously added surfaces, then 2 species should be defined. 
    The diffusion constants and the radii of the 2 species can be different.

    Both species should be created using addSpecies, the latter however requires 
    an explicit 2nd argument to addSpecies referencing the surface it can 
    exist on. They should have different ids.

    A species can not exist on more than 1 surface.

    There is a 3D surface called 'world' by default. Although it is not really 
    a surface, all species and particles that have not been explicitly 
    assigned a surface are assumed to be in the 'world'.

    The method returns the surface it creates. Assign each species to a 
    variable, so it can be used when defining reactions in the next steps.

    '''
    if WORLD:
        Aw = s.addSpecies( 'Aw', D, radius )
        Bw = s.addSpecies( 'Bw', D, radius )

    if MEMBRANE:
        Am = s.addSpecies( 'Am', D, radius, membrane )
        Bm = s.addSpecies( 'Bm', D, radius, membrane )

    if MEMBRANE2:
        Am2 = s.addSpecies( 'Am2', D, radius, membrane2 )
        Bm2 = s.addSpecies( 'Bm2', D, radius, membrane2 )

    if DNA:
        Ad = s.addSpecies( 'Ad', D, radius, dna )
        Bd = s.addSpecies( 'Bd', D, radius, dna )


    ''' Reactions.

    Usage:
        s.addReactionType( SomeReactionType( arguments ) )

    Options (see ReactionTypes in gfrdbase.py):
        - UnimolecularReactionType( reactant, product, k )
        - DecayReactionType( reactant, k )
        - BindingReactionType( reactant1, reactant2, product, k )
        - UnbindingReactionType( reactant, product1, product2, k )

    The reactant and product species for every ReactionType should be on the 
    same surface.

    Combinations of species for which no BindingReactionType is defined are 
    repulsive by default.

    '''
    kf = 1e2
    kb = 1e2

    if WORLD:
        s.addReactionType( BindingReactionType( Aw, Aw, Bw, kf ) )
        s.addReactionType( UnbindingReactionType( Bw, Aw, Aw, kb ) )

    if MEMBRANE:
        s.addReactionType( BindingReactionType( Am, Am, Bm, kf ) )
        s.addReactionType( UnbindingReactionType( Bm, Am, Am, kb ) )

    if MEMBRANE2:
        s.addReactionType( BindingReactionType( Am2, Am2, Bm2, kf ) )
        s.addReactionType( UnbindingReactionType( Bm2, Am2, Am2, kb ) )

    if DNA:
        s.addReactionType( BindingReactionType( Ad, Ad, Bd, kf ) )
        s.addReactionType( UnbindingReactionType( Bd, Ad, Ad, kb ) )


    ''' Surface binding interactions.

    Usage:
        s.addReactionType( SurfaceBindingReactionType( reactant, product, k ) )

    The reactant species for every SurfaceBindingInteractionType should be a 
    species that can only exist in the 'world'. The product species should 
    have been added to the simulator with a second argument defining the 
    surface it can live on.

    When no SurfaceBindingInteractionType is defined for a combination of 
    species and surface, they are repulsive by default.

    '''
    kon = 1e2
    koff = 1e2

    if MEMBRANE and WORLD:
        s.addReactionType( SurfaceBindingReactionType( Aw, Am, kon ) )
        s.addReactionType( SurfaceBindingReactionType( Bw, Bm, kon ) )

    if MEMBRANE2 and WORLD:
        s.addReactionType( SurfaceBindingReactionType( Aw, Am2, kon ) )
        s.addReactionType( SurfaceBindingReactionType( Bw, Bm2, kon ) )

    if DNA and WORLD:
        s.addReactionType( SurfaceBindingReactionType( Aw, Ad, kon ) )
        s.addReactionType( SurfaceBindingReactionType( Bw, Bd, kon ) )


    ''' Surface unbinding reactions.

    Usage:
        s.addReactionType( SurfaceUnBindingReactionType( reactant, product, k ) )

    Unbinding is a Poissonian reaction, from a surface to the 'world'.

    The reactant species should have been added using s.addSpecies with a 2nd 
    argument defining the surface it can exist on. The product species should  
    always be a species that can only exist in the 'world'.

    '''
    if MEMBRANE and WORLD:
        s.addReactionType( SurfaceUnbindingReactionType( Am, Aw, koff ) )
        s.addReactionType( SurfaceUnbindingReactionType( Bm, Bw, koff ) )
    if MEMBRANE2 and WORLD:
        s.addReactionType( SurfaceUnbindingReactionType( Am2, Aw, koff ) )
        s.addReactionType( SurfaceUnbindingReactionType( Bm2, Bw, koff ) )
    if DNA and WORLD:
        s.addReactionType( SurfaceUnbindingReactionType( Ad, Aw, koff ) )
        s.addReactionType( SurfaceUnbindingReactionType( Bd, Bw, koff ) )


    '''
    Particles.

    '''
    if WORLD:
        if MEMBRANE and MEMBRANE2:
            # Add world particles inside the two planes.
            # Note that a CuboidalRegion is defined by 2 corners.
            box1 = CuboidalRegion( [0,0,L*2/10], [L,L,L*8/10] )
            s.throwInParticles( Bw, 2, box1 )
        else:
            # Particles are added to world (3D) by default.
            s.throwInParticles( Bw, 2 )

    if MEMBRANE: s.throwInParticles( Bm, 2, membrane )
    if MEMBRANE2: s.throwInParticles( Bm2, 2, membrane2 )
    if DNA: s.throwInParticles( Bd, 2, dna )


    '''
    Simulation.

    '''
    s.initialize()
    print s.dumpReactions()

    # Write vtk files for last 100 steps only. Use VTK or Paraview to
    # visualize.
    vtklogger = VTKLogger( s, 'run', 100 )

    for i in range( 200 ):
        try:
            vtklogger.log()
            print i
            s.step()
        except RuntimeError, message:
            print message
            break

    vtklogger.stop()
    s.stop( s.t )
    

if __name__ == '__main__':
    run( )

