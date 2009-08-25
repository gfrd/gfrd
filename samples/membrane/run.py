#!/usr/bin/env python

from egfrd import *
from vtklogger import VTKLogger
from logger import *
import time

def run( ):
    '''
    Setting up the simulator.
    '''
    factor = 1#e9
    L = factor * 1e-6	# The world is a micrometer large

    s = EGFRDSimulator(L)#, 'run')

    # s.setWorldSize( L )

    #box1 = CuboidalSurface( [0,0,0], [L,L,L], 'world' )
    #s.cuboidalSurfaces = [ box1 ]


    '''
    Adding the two membranes. The cytosol is 1 micrometer long and 0.6 micrometer in diameter
    '''
    membrane1 = PlanarSurface( [L/2,L/2,2*L/10], [1,0,0], [0,1,0], L/2, L/2, 2.5e-9, 'membrane1' )
    membrane2 = PlanarSurface( [L/2,L/2,8*L/10], [1,0,0], [0,1,0], L/2, L/2, 2.5e-9, 'membrane2' )
    # membrane2 isn't used in any reaction, but is just there for compartmentalization

    s.addSurface( membrane1 )
    s.addSurface( membrane2 )	


    '''
    Species.

    Todo: somehow make sure a species can only be on 1 surface.
    '''
#    signal_mol = Species( 'signal_molecule', 1e-11, 1e-8)	# the signal molecule
#    signal_mol_mem = Species( 'signal_molecule_mem', 1e-12, 1e-8)# the signal molecule on the membrane
    receptor = Species( 'Receptor', 1e-13, 1e-8)		# the membrane bound recepter protein
    receptor_a = Species( 'Receptor_a', 1e-13, 1e-8)		# activated receptor protein (by ligand)
    receptor_a_2 = Species( 'Receptor_a_dimer', 1e-13, 2e-8)	# a dimer of the activated receptor proteins
    carrier = Species( 'Carrier', 2e-12, 1e-8)			# cytosolic carrier protein (TF?)
    carrier_mem = Species( 'Carrier_mem', 2e-13, 1e-8)		# carrier protein associated with the membrane
    carrier_a = Species( 'Carrier_a', 2e-12, 1e-8)		# activated cytosolic carrier protein
    carrier_a_mem = Species( 'Carrier_a_mem', 2e-13, 1e-8)	# active cytosolic carrier protein on membrane
    receptor_carrier_C = Species( 'receptor_carrier_C', 1e-13, 3e-8)	# receptor-carrier complex

#    s.addSpecies( signal_mol )
#    s.addSpecies( signal_mol_mem , membrane1)
    s.addSpecies( receptor , membrane1)
    s.addSpecies( receptor_a , membrane1)
    s.addSpecies( receptor_a_2 , membrane1)
    s.addSpecies( carrier )
    s.addSpecies( carrier_mem , membrane1)
    s.addSpecies( carrier_a )
    s.addSpecies( carrier_a_mem ,membrane1)
    s.addSpecies( receptor_carrier_C ,membrane1)




    '''
    Reactions.
    '''
    '''
    r1 = BindingReactionType( A, B, C, 1e18 )
    r2 = UnbindingReactionType( C, A, B, 1e0 )
    s.addReactionType( r1 )
    s.addReactionType( r2 )
    '''

    # signal_mol + membrane -> association
    # association -> signal_mol + membrane
#    i1 = SurfaceBindingInteractionType( signal_mol, signal_mol_mem, 1e6)
#    r1 = SurfaceUnbindingReactionType ( signal_mol_mem, signal_mol, 1e1)

    # receptor + signal_mol -> activated receptor
    # activated receptor -> receptor
    # SO WE LOSE SIGNAL MOLECULES!
#    r2 = BindingReactionType ( receptor, signal_mol_mem, receptor_a, 1e6)
    r2 = UnimolecularReactionType ( receptor, receptor_a, 1e4)
    r3 = UnimolecularReactionType ( receptor_a, receptor, 1e0)

    # receptor_a + receptor_a -> receptor dimer
    # receptor dimer -> receptor_a + receptor_a						
    r4 = BindingReactionType ( receptor_a, receptor_a, receptor_a_2, 1e4)
    r5 = UnbindingReactionType ( receptor_a_2, receptor_a, receptor_a, 1e1)

    # carrier protein + membrane -> membrane associated carrier
    # membrane associated carrier -> carrier protein
    i2 = SurfaceBindingInteractionType( carrier, carrier_mem, 1e6)
    r6 = SurfaceUnbindingReactionType ( carrier_mem, carrier, 1e1)

    # membrane associated carrier + receptor dimer -> receptor-carrier complex
    # receptor-carrier complex -> receptor dimer + active membrane associated carrier
    r7 = BindingReactionType( carrier_mem, receptor_a_2, receptor_carrier_C, 1e6)
    r8 = UnbindingReactionType ( receptor_carrier_C, receptor_a_2, carrier_a_mem, 1e2)
    r8_2 = UnbindingReactionType ( receptor_carrier_C, receptor_a_2, carrier_mem, 1e1)

    # active carrier in membrane -> active carrier protein
    # active carrier protein + membrane -> association
    r9 = SurfaceUnbindingReactionType ( carrier_a_mem, carrier_a, 1e6)
    i3 = SurfaceBindingInteractionType( carrier_a, carrier_a_mem, 1e1)

    # decay of activated carrier: active carrier -> carrier
    r10 = UnimolecularReactionType ( carrier_a, carrier, 1e-2)

#    s.addInteractionType( i1 )
    s.addInteractionType( i2 )
    s.addInteractionType( i3 )
#    s.addReactionType( r1 )
    s.addReactionType( r2 )
    s.addReactionType( r3 )
    s.addReactionType( r4 )
    s.addReactionType( r5 )
    s.addReactionType( r6 )
    s.addReactionType( r7 )
    s.addReactionType( r8 )
    s.addReactionType( r8_2 )
    s.addReactionType( r9 )
    s.addReactionType( r10 )

    '''
    Putting in particles in the system.
    '''
    box_outside = CuboidalSurface( [0, 0, 0], [L, L, 2*L/10 - 1e-8], 'box_outside' )
    box_inside = CuboidalSurface( [0, 0, 2*L/10 + 1e-8], [L, L, (6*L/10) - 1e-8], 'box_inside' )
    #s.addSurface( box_outside )
    #s.addSurface( box_inside )
    
#    s.throwInParticles( signal_mol, 10, box_outside )
    s.throwInParticles( carrier, 12, box_inside )	# 2 particle is 1nM in this volume
    s.throwInParticles( receptor, 6)

    #s.placeParticle( O, [L/2,L/2,L/2], dna )

    l = Logger( s, 'membrane_stuff' , 'data', 'eerste test')
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') ) ## functie bestaat helemaal niet!
    l.setInterval( 1e-4 )
    l.log()


    '''
    Simulation.
    '''
    s.initialize()
    old = time.time()
    startTime = old
    for i in range(1000000):
        new = time.time()
#        print '\ttime = ', new - old
#        print '\ttotaltime = ', new - startTime
        old = new
        try:
            s.step()
#	    print s.dumpPopulation()
        except Stop, message:
            print message
	    s.vtklogger = VTKLogger(s, 'run_20090824')
            s.vtklogger.log()
            break
	l.log()		## print the number of particle for certain times
    s.stop( s.t )
    
if __name__ == '__main__':
    run( )
