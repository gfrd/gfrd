
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
        self.spheresPrefix = 'spheres'
        self.cylindersPrefix = 'cylinders'
        self.particlesPrefix = 'particles'
        self.fileNamesShells = []    # For .pvd file.
        self.fileNamesCylinders = []    # For .pvd file.
        self.fileNamesParticles = [] # For .pvd file.

        self.i = 0          # Counter.
        self.times = []     # Store times for .pvd file.
        self.deltaT = 1e-11 # Needed for hack.
        # Note: don't make too small, should be relative to max time.
        self.lastTime = 0   # Needed for hack.


    # Write data to file.
    def makeSnapshot( self, prefix, fileList, index, time, doc ):
        fileName = 'steps/' + prefix + str(index) + '.vtu'
        self.vtk_writer.writeDoc(doc, self.dir + fileName, time)
        fileList.append(fileName)


    def log( self ):
        time = self.sim.t
        if (abs(time - self.lastTime) < 1e-9):
            # Hack to make Paraview understand this is a different step but
            # with the same time.
            # 1. During multi global time is not updated.
            # 2. During initialization time is 0.
            # And I want every step recorded, so I can find out what happened 
            # in which step (divide by 2 actually).
            if BROWNIAN != True:
                time = self.lastTime + self.deltaT

        # Get and process data.
        particlesDoc = self.getParticleData( )
        if BROWNIAN != True:
            spheresDoc    = self.getShellData( )
            cylindersDoc = self.getCylinderData( )

            if self.i == 0:
                self.previousShellsDoc = spheresDoc
                self.previousCylindersDoc = cylindersDoc

            # First a snapshot with only the particles updated.
            # Don't use new docs, but previousShellsDoc and previousCylindersDoc.
            self.times.append(time)
            self.makeSnapshot( self.particlesPrefix, self.fileNamesParticles, self.i, time, particlesDoc )
            self.makeSnapshot( self.spheresPrefix, self.fileNamesShells, self.i, time, self.previousShellsDoc )
            self.makeSnapshot( self.cylindersPrefix, self.fileNamesCylinders, self.i, time, self.previousCylindersDoc )
            self.i += 1
            time += self.deltaT  # Hack.
            self.previousShellsDoc = spheresDoc
            self.previousCylindersDoc = cylindersDoc


        # Then a normal snapshot.
        self.times.append(time)
        self.makeSnapshot( self.particlesPrefix, self.fileNamesParticles, self.i, time, particlesDoc )
        if BROWNIAN != True:
            self.makeSnapshot( self.spheresPrefix, self.fileNamesShells, self.i, time, spheresDoc )
            self.makeSnapshot( self.cylindersPrefix, self.fileNamesCylinders, self.i, time, cylindersDoc )
        self.i += 1
        self.lastTime = time


    def stop( self ):
        self.vtk_writer.writePVD(self.dir + self.particlesPrefix + '.pvd', self.fileNamesParticles, self.times)
        self.vtk_writer.writePVD(self.dir + self.spheresPrefix + '.pvd', self.fileNamesShells, self.times)
        self.vtk_writer.writePVD(self.dir + self.cylindersPrefix + '.pvd', self.fileNamesCylinders, self.times)

    
    # Would be nicer to merge this in with getShellData, but I need it in a
    # different .vtu file.
    def getCylinderData( self ):
        posList, radiusList, typeList, lengthList, orientationList = [], [], [], [], []


        # If DNA is along the z-axis in userscript, then in Paraview it should 
        # be drawn along the y-asis. No rotation needed, just set origin to 
        # [500,500,500], height to 10000 and radius to 25.

        """
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
        """



        # Get data from object matrix.
        keys = self.sim.cylinderMatrix.getAll( )
        for key in keys:
            single = key[0]
            #try:
            cylinder = single.shellList[0]
            #except: 
            #    cylinder = single.outside
            type = 1
            scale = numpy.array([cylinder.radius, cylinder.radius, cylinder.size])
            self.appendLists( posList, cylinder.origin, radiusList, cylinder.radius, typeList, type, lengthList, scale, orientationList, [0,0,0] )


        # Get data from scheduler.
        """
        topTime = self.sim.scheduler.getTopTime()
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

        # Get data from object matrix.
        keys = self.sim.sphereMatrix.getAll( )
        for key in keys:
            single = key[0]
            sphere = single.shellList[0]
            type = 1
            self.appendLists( posList, sphere.origin, radiusList, sphere.size, typeList, type )

        """
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
                    if len(object.origin) ==3: # Ugliest hack ever. Todo.
                        self.appendLists( posList, object.origin, radiusList, object.radius, typeList, type )
                except:
                    # This is one of those 1dim spheres that I have to think
                    # about still.
                    pass
                try:
                    # Multi.
                    for single in object.shellList:
                        self.appendLists( posList, single.origin, radiusList, single.size, typeList, type )
                except:
                    pass
        """

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
        factor = 1

        # Multiply radii and lenghts by by 2 because Paraview set radius to 
        # 0.5 by default, and it wants full lengths for cylinders and we are 
        # storing half lengths.
        posList.append( pos * factor )
        radiusList.append( radius * 2 * factor )
        typeList.append( type )
        if scale != None:
            scale *= 2 * factor
            lengthList.append( scale )
        orientationList.append( orientation )
        #orientationList.append( [0,0,0] )


