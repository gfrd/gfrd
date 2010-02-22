#!/usr/bin/env python

'''
To run this script:

PYTHONPATH=../../ python example.py
'''

from egfrd import *
from vtklogger import VTKLogger

'''
The reaction network:

A     <-> B
A + B <-> C
C      -> 0
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


    ''' Adding surfaces. 

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
    MEMBRANE1 = True
    MEMBRANE2 = True
    DNA = True

    if MEMBRANE1:
        m1 = s.addPlanarSurface(origin=[ L / 2, L / 2, 2 * L / 10 ],
                                vectorX=[ 1, 0, 0 ],
                                vectorY=[ 0, 1, 0 ],
                                Lx=(L / 2),
                                Ly=(L / 2),
                                Lz=radius,
                                name='m1')

    if MEMBRANE2:
        m2 = s.addPlanarSurface( origin=[ L / 2, L / 2, 8 * L / 10 ],
                                 vectorX=[ 1, 0, 0 ], 
                                 vectorY=[ 0, 1, 0 ],
                                 Lx=(L / 2),
                                 Ly=(L / 2),
                                 Lz=radius,
                                 name='m2' )

    if DNA:
        d = s.addCylindricalSurface( origin=[ L / 2, L / 2, L / 2 ],
                                     radius=radius,
                                     orientation=[ 0, 1, 0 ],
                                     size=(L / 2),
                                     name='d' )

    ''' Defining species.

    Define the type of particles that can exist.

    Usage:
        species = Species( 'name', D, radius )

    If no D or radius is specified, it has to be specified when adding it to a 
    specific surface.

    '''
    A = Species( 'A', D, radius )
    B = Species( 'B', D, radius )
    C = Species( 'C', D, radius )


    ''' Adding species to 3D space.

    Usage:
        s.addSpecies( species )

    When adding a species to the system without explicitly assigning the 
    surface it can live on, it is assumed to be in the 'world'. The 'world' is 
    a 3D space, also referred to as the 'bulk' or 'cytoplasm'.

    '''
    if WORLD:
        s.addSpecies( A )
        s.addSpecies( B )
        s.addSpecies( C )


    ''' Adding species to surfaces.

    Specify which species can exist on which surface.

    Usage:
        s.addSpecies( species, surface, D, radius )
    See the docstring.

    Per surface a different diffusion constant D and radius can be specified. 
    By default the ones for the 'world' are used. 

    Note: if particles of species 'A' should be able to exist in the world as 
    well as on one of the previously added surfaces, then it should be added 
    twice. Ones with and ones without an explicit surface argument.

    '''
    if MEMBRANE1:
        s.addSpecies( A, m1 )
        s.addSpecies( B, m1 )
        s.addSpecies( C, m1 )

    if MEMBRANE2:
        # No species can live on membrane 2.
        pass

    if DNA:
        s.addSpecies( A, d )
        s.addSpecies( B, d, D * 2 )
        s.addSpecies( C, d, D * 2, radius * 2 )


    ''' Adding reactions in 3D.

    Usage:
        s.addReaction( [reactants], [products], rate )

    For now: a bimolecular reaction can only have 1 product species.

    Units for the rate of a bimolecular reaction:
        [kf] = meters^3 / (molecules * second)

    (Order of magnitude kf: 1e-18)

    '''
    sigma = 2 * radius
    D_tot = D
    tau = sigma * sigma / D_tot
    # Bimolecular. 
    kf_2 = 100 * sigma * D_tot
    kb_2 = 0.1 / tau
    # Unimolecular.
    kf_1 = 0.1 / tau
    kb_1 = 0.1 / tau

    if WORLD:
        # A     <-> B
        # A + B <-> C
        # C      -> 0
        s.addReaction( [ A    ], [ B    ], kf_1 )
        s.addReaction( [ B    ], [ A    ], kb_1 )
        s.addReaction( [ A, B ], [ C    ], kf_2 )
        s.addReaction( [ C    ], [ A, B ], kb_2 )
        s.addReaction( [ C    ], [      ], kf_1 )


    ''' Adding reactions on surfaces.

    Usage:
        s.addReaction( [reactants], [products], rate )
    , where one reactant or product is a tuple: (species, surface).

    Combinations of species which do not appear together as reactants in any 
    reaction are made repulsive by default on every surface.

    '''
    if MEMBRANE1:
        # A     <-> B
        # A + B <-> C
        # C      -> 0
        s.addReaction( [ (A, m1),          ], [ (B, m1)           ], kf_1 )
        s.addReaction( [ (B, m1),          ], [ (A, m1)           ], kb_1 )
        s.addReaction( [ (A, m1),  (B, m1) ], [ (C, m1)           ], kf_2 )
        s.addReaction( [ (C, m1)           ], [ (A, m1),  (B, m1) ], kb_2 )
        s.addReaction( [ (C, m1),          ], [                   ], kf_1 )

    if MEMBRANE2:
        # No particles can live on membrane2.
        pass

    if DNA:
        # A     <-> B
        # A + B <-> C
        # C      -> 0
        s.addReaction( [ (A, d),          ], [ (B, d)           ], kf_1 )
        s.addReaction( [ (B, d),          ], [ (A, d)           ], kb_1 )
        s.addReaction( [ (A, d),  (B, d)  ], [ (C, d)           ], kf_2 )
        s.addReaction( [ (C, d)           ], [ (A, d),  (B, d)  ], kb_2 )
        s.addReaction( [ (C, d),          ], [                  ], kf_1 )


    ''' Surface binding reactions.

    Usage:
        s.addReaction( [reactant], [product], rate ) )
    , where product is a tuple: (species, surface).

    The reactant species for every surface binding reaction is always a 
    species that can only exist in the 'world', so no surface has to be 
    specified there.

    If a surface should be absorbive, specify (0, surface) as the product.

    When no surface binding reaction is defined for a combination of a species 
    and a surface, they are made repulsive by default.

    '''
    # Molecule + surface.
    kon = kf_2

    if MEMBRANE1 and WORLD:
        # Species C can bind to the membrane. The membrane is reflective, by 
        # default, to species A and B.
        s.addReaction( [ C ],  [ (C, m1)  ], kon )

    if MEMBRANE2 and WORLD:
        # Membrane 2 absorbs all particles.
        s.addReaction( [ A ],  [ (0, m2) ], kon )
        s.addReaction( [ B ],  [ (0, m2) ], kon )
        s.addReaction( [ C ],  [ (0, m2) ], kon )
        pass

    if DNA and WORLD:
        # Species C can bind to the dna. The dna is reflective, by default, to 
        # species A and B.
        s.addReaction( [ C ],  [ (C, d)  ], kon )


    ''' Surface unbinding reactions.

    Usage:
        s.addReaction( [reactant], [product], k ) )
    , where reactant is a tuple: (species, surface).

    Unbinding is a Poissonian reaction, from a surface to the 'world'.

    The product species for every surface binding reaction is always a 
    species that can only exist in the 'world', so no surface has to be 
    specified there.

    '''
    # Unimolecular.
    koff = kb_1

    if MEMBRANE1 and WORLD:
        # Species C can unbind to the 'world'.
        s.addReaction( [ (C, m1)  ], [ C ], koff )

    if MEMBRANE2 and WORLD:
        # No particles can live on membrane2.
        pass

    if DNA and WORLD:
        # Species C can unbind to the 'world'.
        s.addReaction( [ (C, d)  ], [ C ], koff )


    '''
    Particles.

    '''
    if WORLD:
        if MEMBRANE1 and MEMBRANE2:
            # Add world particles inside the two planes.
            # Note that a CuboidalRegion is defined by 2 corners.
            box1 = CuboidalRegion( [ 0, 0, 2 * L / 10 ], [ L, L, 8 * L / 10 ] )
            s.throwInParticles( C, 2, box1 )
        else:
            # Particles are added to world (3D) by default.
            s.throwInParticles( C, 2 )

    if MEMBRANE1: s.throwInParticles( C, 2, m1 )
    if MEMBRANE2: pass
    if DNA: s.throwInParticles( C, 2, d )


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

