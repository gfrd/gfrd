
import os
#import re
#import logging
import numpy
from vtk_xml_serial_unstructured import VTK_XML_Serial_Unstructured
from shape import DummyCylinder, DummyBox, Cylinder, Box

INF = numpy.inf
BROWNIAN = False


# Logger that can be used to visualize data with Kitware ParaView.
# log():
#   (get.......Data = createDoc() + addPiece()) + writeDoc()
#
# stop():
#   writePvd(name.pvd)
class VTKLogger:
    def __init__(self, sim, name):
        self.sim = sim
        self.vtk_writer = VTK_XML_Serial_Unstructured()

        # Filename stuff.
        self.name = name
        if not os.path.exists('data/' + self.name + '/files'):
            os.makedirs('data/' + self.name + '/files')

        self.fileList = {'particles':[], 'spheres':[], 'cylinders':[], 
                'cylindricalSurfaces':[], 'planarSurfaces':[]}  # For .pvd file.

        self.i = 0          # Counter.
        self.deltaT = 1e-11 # Needed for hack. Note: don't make too small, 
                            # should be relative to max time.
        self.lastTime = 0   # Needed for hack.


    def log( self ):
        time = self.sim.t
        #print 'log called, time = ', time
        if (abs(time - self.lastTime) < 1e-9):
            # Hack to make Paraview understand this is a different simulator 
            # step but with the same time.
            # 1. During multi global time is not updated.
            # 2. During initialization time is 0.
            # 3. Sometimes shells have mobilityRadius = 0 --> dt=0.
            # And I want every step recorded, so I can find out what happened 
            # in which step (divide by 2 actually).
            if BROWNIAN != True:
                # Now Paraview should perceive a state change.
                time = self.lastTime + self.deltaT

        # Get and process data.
        particlesDoc = self.getParticleData()
        cylindricalSurfacesDoc = self.getCylindricalSurfaceData()
        planarSurfacesDoc = self.getPlanarSurfaceData()
        if BROWNIAN != True:
            spheresDoc   = self.getShellData( )
            cylindersDoc = self.getCylinderData( )

            if self.i == 0:
                self.previousShellsDoc = spheresDoc
                self.previousCylindersDoc = cylindersDoc

            # First a snapshot with only the particles updated. Don't use new 
            # docs, but previousShellsDoc and previousCylindersDoc.
            self.makeSnapshot( 'particles', self.i, time, particlesDoc )
            self.makeSnapshot( 'spheres', self.i, time, self.previousShellsDoc )
            self.makeSnapshot( 'cylinders', self.i, time, self.previousCylindersDoc )
            self.makeSnapshot( 'cylindricalSurfaces', self.i, time, cylindricalSurfacesDoc )
            self.makeSnapshot( 'planarSurfaces', self.i, time, planarSurfacesDoc )
            self.i += 1
            time += self.deltaT  # Hack.
            self.previousShellsDoc = spheresDoc
            self.previousCylindersDoc = cylindersDoc


        # Then a normal snapshot.
        self.makeSnapshot( 'particles', self.i, time, particlesDoc )
        self.makeSnapshot( 'cylindricalSurfaces', self.i, time, cylindricalSurfacesDoc )
        self.makeSnapshot( 'planarSurfaces', self.i, time, planarSurfacesDoc )
        if BROWNIAN != True:
            self.makeSnapshot( 'spheres', self.i, time, self.previousShellsDoc )
            self.makeSnapshot( 'cylinders', self.i, time, self.previousCylindersDoc )
        self.i += 1
        self.lastTime = time


    # Write data to file.
    # Store filename and time in fileList.
    def makeSnapshot( self, type, index, time, doc ):
        fileName = 'files/' + type + str(index) + '.vtu'
        self.vtk_writer.writeDoc(doc, 'data/' + self.name + '/' + fileName)
        self.fileList[type].append( (fileName, time) )


    def stop( self ):
        self.vtk_writer.writePVD('data/' + self.name + '/' + 'files.pvd', self.fileList)


    def getParticleData( self ):
        posList, radiusList, typeList = [], [], []

        for speciesIndex,species in enumerate(self.sim.speciesList.values()):
            for particlePos in species.pool.positions:
                self.appendLists( posList, particlePos, radiusList,
                        species.radius, typeList, speciesIndex * 2)

        return self.vtk_writer.createDoc(posList, radii=radiusList, 
                colors=typeList )


    def getShellData( self ):
        posList, radiusList, typeList = [], [], []

        # Get data from object matrix.
        keys = self.sim.sphereMatrix.getAll( )
        for key in keys:
            single = key[0]
            sphere = single.shellList[0]
            type = 1
            self.appendLists( posList, sphere.origin, radiusList, sphere.size, 
                    typeList, type )

        return self.vtk_writer.createDoc(posList, radii=radiusList, 
                colors=typeList )


    def getCylinderData( self ):
        # Get data from object matrix.
        keys = self.sim.cylinderMatrix.getAll( )
        cylinders = [ key[0].shellList[0] for key in keys ]
        return self.processCylinders( cylinders )


    def getCylindricalSurfaceData( self ):
        cylinders = [ surface.outside for surface in self.sim.surfaceList if isinstance(surface.outside, Cylinder) ]
        return self.processCylinders( cylinders )

    def getPlanarSurfaceData( self ):
        posList, radiusList, typeList, scaleList, orientationList, tensorList = \
                [], [], [], [], [], []

        boxes = [ surface.outside for surface in self.sim.surfaceList if isinstance(surface.outside, Box) ]
        if len(boxes) == 0:
            # Add dummy box to stop tensorGlyph from complaining.
            boxes = [ DummyBox() ] 

        type = 1
        for box in boxes:
            #print 'box = ', box
            tensor = numpy.concatenate((box.xVector, box.yVector, box.zVector))
            #print 'tensor = ', tensor
            self.appendLists(posList, box.origin, typeList, type, tensorList=tensorList, 
                    tensor=tensor) 

        return self.vtk_writer.createDoc(posList, colors=typeList, tensors=tensorList )


    def processCylinders( self, cylinders=[] ):
        posList, radiusList, typeList, scaleList, orientationList, tensorList = \
                [], [], [], [], [], []

        if len(cylinders) == 0:
            # Add dummy cylinder to stop tensorGlyph from complaining.
            cylinders = [ DummyCylinder() ] 

        for cylinder in cylinders:
            radius = cylinder.radius
            type = 1
            orientation = cylinder.orientation
            size = cylinder.size

            # Construct scale vector for scaling by vector components. Not 
            # used anymore.
            scale = numpy.array([radius, radius, size])

            # Construct tensor. Use tensor glyph plugin from:
            # http://www.paraview.org/pipermail/paraview/2009-March/011256.html
            # Unset Extract eigenvalues.

            # Select basis vector in which orientation is smallest.
            _, basisVector = min( zip(orientation, [[1,0,0], [0,1,0], 
            [0,0,1]]) )
            # Find 2 vectors perpendicular to orientation.
            perpendicular1 = numpy.cross( orientation, basisVector )
            perpendicular2 = numpy.cross( orientation, perpendicular1 )
            # A 'tensor' is represented as an array of 9 values.
            # Stupid Paraview wants  a normal vector to the cylinder to orient  
            # it. So orientation and perpendicular1 swapped.
            tensor = numpy.concatenate((perpendicular1*radius, 
                orientation*size, perpendicular2*radius))

            self.appendLists( posList, cylinder.origin, radiusList,
                    size, typeList, type, scaleList, scale, 
                    orientationList, orientation, tensorList, tensor ) 

        return self.vtk_writer.createDoc(posList, radii=radiusList, 
                colors=typeList, orientations=orientationList, 
                scales=scaleList, tensors=tensorList )

    

    # Helper.
    def appendLists(self, posList, pos, radiusList=[], radius=None, 
            typeList=[], type=None, scaleList=[], scale=None, 
            orientationList=[], orientation=None, tensorList=[], tensor=None):
        # Convert all lengths to nanometers.
        #factor = 1e9
        factor = 1

        # Multiply radii and lenghts by by 2 because Paraview sets radius to 
        # 0.5 by default, and it wants full lengths for cylinders and we are 
        # storing half lengths.
        posList.append( pos * factor )
        radiusList.append( radius * 2 * factor )
        typeList.append( type )

        if orientation != None:
            orientationList.append( orientation )

        if scale != None:
            scaleList.append( scale * 2 * factor )

        if tensor != None:
            #print 'tensor = ', tensor
            tensor = tensor * 2 * factor
            #print 'tensor = ', tensor
            
            tensorList.append( tensor )




        """
        # If DNA is along the z-axis in userscript, then in Paraview it should 
        # be drawn along the y-asis. No rotation needed, just set origin to 
        # [500,500,500], height to 10000 and radius to 25.
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
        #self.appendLists( posList, pos, radiusList, 0 ,
        #        typeList, 0, lengthList, scale, orientationList, [0,0,0] )
        """

        """
        # Get data from scheduler.
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
                scale = numpy.array([object.radius_b, object.radius_b, 
                    length])
                self.appendLists( posList, object.pos, radiusList, 
                        object.radius_b ,
                        typeList, type, lengthList, scale, orientationList, 
                        object.orientation )
            except AttributeError:
                # So this is not a cylinder.
                pass
        """

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
                        self.appendLists( posList, object.origin, radiusList, 
                                object.radius, typeList, type )
                except:
                    # This is one of those 1dim spheres that I have to think
                    # about still.
                    pass
                try:
                    # Multi.
                    for single in object.shellList:
                        self.appendLists( posList, single.origin, radiusList, 
                                single.size, typeList, type )
                except:
                    pass
        """

        """
        self.vtk_writer.writePVD(self.dir + self.particlesPrefix + '.pvd', 
                [self.fileNamesParticles], self.times)
        self.vtk_writer.writePVD(self.dir + self.spheresPrefix + '.pvd', 
                [self.fileNamesSpheres], self.times)
        self.vtk_writer.writePVD(self.dir + self.cylindersPrefix + '.pvd', 
                [self.fileNamesCylinders], self.times)
        """
