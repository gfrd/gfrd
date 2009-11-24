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
        membrane = s.addPlanarSurface( origin=[ L / 2, L / 2, 2 * L / 10 ],
                                       vectorX=[ 1, 0, 0 ],
                                       vectorY=[ 0, 1, 0 ],
                                       Lx=(L / 2),
                                       Ly=(L / 2),
                                       Lz=radius,
                                       name='membrane' )

    if MEMBRANE2:
        membrane2 = s.addPlanarSurface( origin=[ L / 2, L / 2, 8 * L / 10 ],
                                        vectorX=[ 1, 0, 0 ], 
                                        vectorY=[ 0, 1, 0 ],
                                        Lx=(L / 2),
                                        Ly=(L / 2),
                                        Lz=radius,
                                        name='membrane2' )

    if DNA:
        dna = s.addCylindricalSurface( origin=[ L / 2, L / 2, L / 2 ],
                                       radius=radius,
                                       orientation=[ 0, 1, 0 ],
                                       size=(L / 2),
                                       name='dna' )

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
        s.addReaction( [reactants], [products], rate )

    The reactant and product species for every Reaction should be on the same 
    surface.

    Combinations of species which do not appear together as reactants in any 
    reaction are made repulsive by default.

    '''
    kf = 1e2
    kb = 1e2

    if WORLD:
        s.addReaction( [ Aw, Aw ],   [ Bw ],       kf )
        s.addReaction( [ Bw ],       [ Aw, Aw ],   kb )

    if MEMBRANE:
        s.addReaction( [ Am, Am ],   [ Bm ],       kf )
        s.addReaction( [ Bm ],       [ Am, Am ],   kb )

    if MEMBRANE2:
        s.addReaction( [ Am2, Am2 ], [ Bm2 ],      kf )
        s.addReaction( [ Bm2 ],      [ Am2, Am2 ], kb )

    if DNA:
        s.addReaction( [ Ad, Ad ],   [ Bd ],       kf )
        s.addReaction( [ Bd ],       [ Ad, Ad ],   kb )


    ''' Surface binding reactions.

    Usage:
        s.addReaction( [reactant], [product], rate ) )

    The reactant species for every surface binding reaction should be a 
    species that can only exist in the 'world'. The product species should 
    have been added to the simulator with a second argument defining the 
    surface it can live on.

    When no surface binding reaction is defined for a combination of a species 
    and a surface, they are made repulsive by default.

    '''
    kon = 1e2
    koff = 1e2

    if MEMBRANE and WORLD:
        s.addReaction( [ Aw ],  [ Am ],  kon )
        s.addReaction( [ Bw ],  [ Bm ],  kon )

    if MEMBRANE2 and WORLD:
        s.addReaction( [ Am2 ], [ Am2 ], kon )
        s.addReaction( [ Bm2 ], [ Bm2 ], kon )

    if DNA and WORLD:
        s.addReaction( [ Aw ],  [ Ad ],  kon )
        s.addReaction( [ Bw ],  [ Bd ],  kon )


    ''' Surface unbinding reactions.

    Usage:
        s.addReaction( [reactant], [product], k ) )

    Unbinding is a Poissonian reaction, from a surface to the 'world'.

    The reactant species should have been added using s.addSpecies with a 2nd 
    argument defining the surface it can exist on. The product species should  
    always be a species that can only exist in the 'world'.

    '''
    if MEMBRANE and WORLD:
        s.addReaction( [ Am ], [ Aw ], koff )
        s.addReaction( [ Bm ], [ Bw ], koff )
    if MEMBRANE2 and WORLD:
        s.addReaction( [ Am2 ], [ Aw ], koff )
        s.addReaction( [ Bm2 ], [ Bw ], koff )
    if DNA and WORLD:
        s.addReaction( [ Ad ], [ Aw ], koff )
        s.addReaction( [ Bd ], [ Bw ], koff )


    '''
    Particles.

    '''
    if WORLD:
        if MEMBRANE and MEMBRANE2:
            # Add world particles inside the two planes.
            # Note that a CuboidalRegion is defined by 2 corners.
            box1 = CuboidalRegion( [ 0, 0, 2 * L / 10 ], [ L, L, 8 * L / 10 ] )
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
            s.step()
        except RuntimeError, message:
            print message
            break

    vtklogger.stop()
    s.stop( s.t )
    

if __name__ == '__main__':
    run( )

