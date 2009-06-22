
import os
#import re
#import logging
import numpy

INF = numpy.inf
PARTICLES = 0
SHELLS = 1
CYLINDERS = 2

BROWNIAN = False
from vtk_xml_serial_unstructured import VTK_XML_Serial_Unstructured

# Logger that can be used to visualize data with Kitware ParaView.
class VTKLogger:

    def __init__(self, sim, dir):
        self.sim = sim
        self.vtk_writer = VTK_XML_Serial_Unstructured()

        # Filename stuff.
        datadir = 'data/' + dir + '/steps/'
        if not os.path.exists(datadir):
            os.makedirs(datadir)

        self.dir = 'data/' + dir + '/'
        self.shellsPrefix = 'shells'
        self.cylindersPrefix = 'cylinders'
        self.particlesPrefix = 'particles'
        self.fileNamesShells = []    # For .pvd file.
        self.fileNamesCylinders = []    # For .pvd file.
        self.fileNamesParticles = [] # For .pvd file.

        self.i = 0          # Counter.
        self.times = []     # Store times for .pvd file.
        self.deltaT = 1e-13 # Needed for hack.
        self.lastTime = 0   # Needed for hack.


    # Write data to file.
    def makeSnapshot( self, prefix, fileList, index, time, doc ):
        fileName = 'steps/' + prefix + str(index) + '.vtu'
        self.vtk_writer.writeDoc(doc, self.dir + fileName, time)
        fileList.append(fileName)


    def log( self ):
        time = self.sim.t
        if (abs(time - self.lastTime) < 1e-10):
            # Hack to make Paraview understand this is a different step but
            # with the same time.
            if BROWNIAN != True:
                time = self.lastTime + self.deltaT

        # Get and process data.
        particlesDoc = self.getParticleData( )
        if BROWNIAN != True:
            shellsDoc    = self.getShellData( )
            cylindersDoc = self.getCylinderData( )

            if self.i == 0:
                self.previousShellsDoc = shellsDoc
                self.previousCylindersDoc = cylindersDoc

            # First a snapshot with only the particles updated.
            # Don't use new docs, but previousShellsDoc and previousCylindersDoc.
            self.times.append(time)
            self.makeSnapshot( self.particlesPrefix, self.fileNamesParticles, self.i, time, particlesDoc )
            self.makeSnapshot( self.shellsPrefix, self.fileNamesShells, self.i, time, self.previousShellsDoc )
            self.makeSnapshot( self.cylindersPrefix, self.fileNamesCylinders, self.i, time, self.previousCylindersDoc )
            self.i += 1
            time += self.deltaT  # Hack.
            self.previousShellsDoc = shellsDoc
            self.previousCylindersDoc = cylindersDoc


        # Then a normal snapshot.
        self.times.append(time)
        self.makeSnapshot( self.particlesPrefix, self.fileNamesParticles, self.i, time, particlesDoc )
        if BROWNIAN != True:
            self.makeSnapshot( self.shellsPrefix, self.fileNamesShells, self.i, time, shellsDoc )
            self.makeSnapshot( self.cylindersPrefix, self.fileNamesCylinders, self.i, time, cylindersDoc )
        self.i += 1
        self.lastTime = time


    def stop( self ):
        self.vtk_writer.writePVD(self.dir + self.particlesPrefix + '.pvd', self.fileNamesParticles, self.times)
        self.vtk_writer.writePVD(self.dir + self.shellsPrefix + '.pvd', self.fileNamesShells, self.times)
        self.vtk_writer.writePVD(self.dir + self.cylindersPrefix + '.pvd', self.fileNamesCylinders, self.times)

    
    # Would be nicer to merge this in with getShellData, but I need it in a
    # different .vtu file.
    def getCylinderData( self ):
        posList, radiusList, typeList, lengthList, orientationList = [], [], [], [], []

        topTime = self.sim.scheduler.getTopTime()

        # Add DNA by default.
        dnaL = 1e-6
        sigma = 2.5e-8
        dnaR = sigma
        # Paraview wants to know the middle of the cylinder.
        pos = numpy.array( [dnaL/2, dnaL/2, dnaL/2] )
        # Note: when scaling, setting radius has no effect.
        # Note: set scalemode to vector components!
        # Note: set radius to 1 in Paraview, also for particles.
        scale = numpy.array( [dnaR, dnaR, dnaL])
        # Orientationlist should contain the number of degrees [0-360] of 
        # rotation in each direction. Not used. Commented out in 
        # vtk_xml_serial_unstructured.

        # Not by default at the moment.
        #self.appendLists( posList, pos, radiusList, 0 , \
        #        typeList, 0, lengthList, scale, orientationList, [0,0,0] )



        # Get data from object matrix.
        cylinders = self.sim.shellMatrix.getCylinders( [0,0,0] )
        for cylinder in cylinders:
            type = 1
            # Multiply cylinder.size by 2 because we are storing halfLengths 
            # and Paraview wants full length.
            scale = numpy.array([cylinder.radius, cylinder.radius, cylinder.size * 2])
            self.appendLists( posList, cylinder.pos, radiusList, cylinder.radius, typeList, type, lengthList, scale, orientationList, [0,0,0] )


        # Get data from scheduler.
        """
        for eventIndex in range( self.sim.scheduler.getSize() ):
            # Get event
            event = self.sim.scheduler.getEventByIndex( eventIndex )
            object = event.getArg()

            type = 1 # Color

            # Highlight topEvent.
            if event.getTime() == topTime:
                type = 0

            try:
                length = object.halfLength  # Only cylinders have length.
                scale = numpy.array([object.radius_b, object.radius_b, length])
                self.appendLists( posList, object.pos, radiusList, object.radius_b , \
                        typeList, type, lengthList, scale, orientationList, object.orientation )
            except AttributeError:
                # So this is not a cylinder.
                pass
        """

        # Initialize XML doc.
        doc, grid = self.vtk_writer.createDoc()
        # Add data as XML child to grid.
        return self.vtk_writer.addPiece(doc, grid, posList, orientations=orientationList, \
                lengths=lengthList, radii=radiusList, colors=typeList )


    def getShellData( self ):
        posList, radiusList, typeList, lengthList = [], [], [], []

        topTime = self.sim.scheduler.getTopTime()

        for eventIndex in range( self.sim.scheduler.getSize() ):
            # Get event
            event = self.sim.scheduler.getEventByIndex( eventIndex )
            object = event.getArg()

            # Give different color to singles/paires/multis.
            type = object.multiplicity
            # Highlight topEvent.
            if event.getTime() == topTime:
                pass #type = 0

            try:
                # Don't show sphere for cylinder.
                length = object.halfLength  # Only cylinders have length.
            except:
                # Fixme. How do I find if this object is a single/pair/multi?
                # Single or Pair.
                try:
                    if len(object.pos) ==3: # Ugliest hack ever. Todo.
                        self.appendLists( posList, object.pos, radiusList, object.radius, typeList, type )
                except:
                    # This is one of those 1dim shells that I have to think
                    # about still.
                    pass
                try:
                    # Multi.
                    for single in object.shellList:
                        self.appendLists( posList, single.pos, radiusList, single.size, typeList, type )
                except:
                    pass

        # Initialize XML doc.
        doc, grid = self.vtk_writer.createDoc()
        # Add data as XML child to grid.
        return self.vtk_writer.addPiece(doc, grid, posList, radii=radiusList, colors=typeList )


    def getParticleData( self ):
        posList, radiusList, typeList, lengthList = [], [], [], []

        for speciesIndex,species in enumerate( self.sim.speciesList.values() ):
            for particlePos in species.pool.positions:
                self.appendLists( posList, particlePos, radiusList, species.radius, typeList, speciesIndex * 2)

        # Initialize XML doc.
        doc, grid = self.vtk_writer.createDoc()
        # Add data as XML child to grid.
        return self.vtk_writer.addPiece(doc, grid, posList, radii=radiusList, colors=typeList )


    # Helper.
    def appendLists( self, posList, pos, radiusList, radius, typeList, type, lengthList=[], scale=None, orientationList=[], orientation=None ):
        # Convert all lengths to nanometers.
        factor = 1e9

        posList.append( pos * factor )
        radiusList.append( radius * factor )
        typeList.append( type )
        if scale != None:
            scale *= factor
            lengthList.append( scale )
        #orientationList.append( orientation )
        orientationList.append( [0,0,0] )


