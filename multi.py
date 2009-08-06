from bd import *
from shell import *


'''
Used internally by Multi.
'''

### This is the simulator that is used for multis.
class MultiBDCore( BDSimulatorCoreBase ):
    
    def __init__( self, main, multi ):

        BDSimulatorCoreBase.__init__( self, main )

        # this has to be ref, not proxy, since it is used for comparison.
        self.multiref = weakref.ref( multi )

        self.particleMatrix = SphereMatrix()
        self.particleMatrix.setWorldSize( self.main.worldSize )

        self.shellMatrix = SphereMatrix()
        self.shellMatrix.setWorldSize( self.main.worldSize )

        self.escaped = False


        
    def updateParticle( self, particle, pos ):

        self.particleMatrix.update( particle, pos, particle.radius )
        # Todo. Is it ok to use moveParticle here?
        self.main.moveParticle( particle, pos )

    def initialize( self ):

        BDSimulatorCoreBase.initialize( self )
        self.updateShellMatrix()


    def step( self ):

        self.escaped = False
        BDSimulatorCoreBase.step( self )


    '''
    def sync( self ):
        for particle in self.particleList:
            self.main.updateOnParticleMatrix( particle, particle.pos )
'''

    def updateShellMatrix( self ):

        self.shellMatrix.clear()
        for shell in self.multiref().shellList:
            self.shellMatrix.add( shell, shell.origin, shell.radius )

    def addParticle( self, particle ):

        self.addToParticleList( particle )
        self.particleMatrix.add( particle,
                                 particle.pos, particle.radius )

    def removeParticle( self, particle ):

        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )
        self.particleMatrix.remove( particle )


    def createParticle( self, species, pos, surface ):

        #if not self.withinShell( pos, species.radius ):
        #    self.escaped = True
        #    self.clearOuterVolume( pos, species.radius )

        particle = self.main.createParticle( species, pos, surface )
        self.addParticle( particle )

        return particle


    def moveParticle( self, particle, pos ):

        particle.pos = pos
        self.updateParticle( particle, pos )

        
    def clearVolume( self, pos, radius, ignore=[] ):

        if not self.withinShell( pos, radius ):
            self.escaped = True
            self.clearOuterVolume( pos, radius, ignore )


    def clearOuterVolume( self, pos, radius, ignore=[] ):

        self.main.clearVolume( pos, radius, ignore=[self.multiref(),] )
        if not self.main.checkOverlap( pos, radius, ignore ):
            raise NoSpace()

 
    def withinShell( self, pos, radius ):

        n, _ = self.shellMatrix.getNeighborsWithinRadiusNoSort( pos, - radius )
        return n

        
    def checkOverlap( self, pos, radius, ignore=[] ):
        
        n, _ = self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        n = list( n )
        for particle in ignore:
            if particle in n:
                n.remove( particle )

        return not n

    '''
    def getNeighborParticles( self, pos, n=None ):

        n, d = self.particleMatrix.getNeighbors( pos, n )
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d
    '''

    def getParticlesWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        neighbors, _ = \
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
        return [ n for n in neighbors if n not in ignore ]


    ### Why don't we borrow this stuff from ParticleSimulatorBase?
    def getClosestParticle( self, pos, ignore=[] ):

        neighbors, distances =\
            self.particleMatrix.getNeighbors( pos, len( ignore ) + 1 )

        for i, neighbor in enumerate( neighbors ):
            particle = Particle( neighbor[0], neighbor[1] )
            if particle not in ignore:
                return particle, distances[i]

        # default case: none left.
        return None, INF


    def check( self ):

        BDSimulatorCoreBase.check( self )

        # shellMatrix consistency
        for shell in self.multiref().shellList:
            pos, radius = self.shellMatrix.get( shell )
            assert not ( pos - shell.origin ).any()
            assert radius == shell.radius


        # shells are contiguous
        for shell in self.multiref().shellList:
            n, d = self.shellMatrix.getNeighbors( shell.origin )
            assert d[1] - shell.radius < 0.0, 'shells are not contiguous.'

        # all particles within the shell.
                
        for p in self.particleList:
            assert self.withinShell( p.pos, p.species.radius ),\
                'not all particles within the shell.'


### How can this be so much shorter than the Pair class?
class Multi( object ):

    def __init__( self, main ):

        ### A multi does not have a shell of itself.
        self.shellList = []
        self.eventID = None

        ### Each multi has it's own BD simulator!
        self.sim = MultiBDCore( main, self )



    def initialize( self, t ):

        self.lastTime = t
        self.startTime = t

        self.sim.initialize() # ??


    def getDt( self ):
        return self.sim.dt

    dt = property( getDt )


    def getMultiplicity( self ):
        return len( self.sim.particleList )

    multiplicity = property( getMultiplicity )

    ### So a multi contains particles+their shells, but not singles.
    def addParticle( self, particle ):
        self.sim.addParticle( particle )

    def addShell( self, origin, radius ):
        self.shellList.append( Sphere( origin, radius ) )


    def check( self ):

        self.sim.check()


    def __str__( self ):

        if len( self.sim.particleList ) == 0:
            return 'Multi()'

        buf = 'Multi( '
        buf += str( self.sim.particleList[0] )
        for particle in self.sim.particleList[1:]:
            buf += ', ' + str( particle )
        buf += ' )'

        return buf
