
import os
#import re
#import logging
import numpy
from vtk_xml_serial_unstructured import *
from shape import *

INF = numpy.inf
BROWNIAN = False

class DummyCylinder( Cylinder ):
    def __init__( self ):
        Cylinder.__init__( self, [0,0,0], 1e-20, [0,0,1], 1e-20 ) 


class DummyBox( Box ):
    def __init__( self ):
        Box.__init__( self, [0,0,0], [1,0,0], [0,1,0], [0,0,1], 1e-20, 1e-20, 1e-20 ) 


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
                'cuboidalSurfaces':[], 'cylindricalSurfaces':[], 
                'planarSurfaces':[]}  # For .pvd file.

        self.i = 0          # Counter.
        self.deltaT = 1e-11 # Needed for hack. Note: don't make too small, 
                            # should be relative to max time.
        self.lastTime = 0   # Needed for hack.


    def log( self ):
        time = self.sim.t
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
        cuboidalSurfacesDoc = self.getCuboidalSurfaceData()
        cylindricalSurfacesDoc = self.getCylindricalSurfaceData()
        planarSurfacesDoc = self.getPlanarSurfaceData()
        if BROWNIAN != True:
            spheresDoc   = self.getSphericalShellDataFromScheduler( )
            cylindersDoc = self.getCylindricalShellData( )

            if self.i == 0:
                self.previousShellsDoc = spheresDoc
                self.previousCylindersDoc = cylindersDoc

            # First a snapshot with only the particles updated. Don't use new 
            # docs, but previousShellsDoc and previousCylindersDoc.
            self.makeSnapshot( 'particles', self.i, time, particlesDoc )
            self.makeSnapshot( 'spheres', self.i, time, self.previousShellsDoc )
            self.makeSnapshot( 'cylinders', self.i, time, self.previousCylindersDoc )
            self.makeSnapshot( 'cuboidalSurfaces', self.i, time, cuboidalSurfacesDoc )
            self.makeSnapshot( 'cylindricalSurfaces', self.i, time, cylindricalSurfacesDoc )
            self.makeSnapshot( 'planarSurfaces', self.i, time, planarSurfacesDoc )
            self.i += 1
            time += self.deltaT  # Hack.
            self.previousShellsDoc = spheresDoc
            self.previousCylindersDoc = cylindersDoc


        # Then a normal snapshot.
        self.makeSnapshot( 'particles', self.i, time, particlesDoc )
        self.makeSnapshot( 'cuboidalSurfaces', self.i, time, cuboidalSurfacesDoc )
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
                self.appendLists( posList, particlePos, typeList, 2 + speciesIndex, radiusList, species.radius)

        return self.vtk_writer.createDoc(posList, radii=radiusList, 
                colors=typeList )


    def getSphericalShellDataFromScheduler( self ):
        posList, radiusList, typeList = [], [], []

        topTime = self.sim.scheduler.getTopTime()
        for eventIndex in range( self.sim.scheduler.getSize() ):
            # Get event
            event = self.sim.scheduler.getEventByIndex( eventIndex )
            object = event.getArg()

            # Give different color to singles/paires/multis.
            type = object.multiplicity
            # Highlight topEvent.
            if event.getTime() == topTime:
                type = 0

            try:
                # Don't show sphere for cylinder.
                object.shellList[0].size  # Only cylinders have a 'size'.
            except:
                # Single or Pair or Multi.
                for single in object.shellList:
                    self.appendLists( posList, single.origin, typeList, type, radiusList, single.radius )

        return self.vtk_writer.createDoc(posList, radii=radiusList, 
                colors=typeList )


    def getCylindricalShellData( self ):
        # Get data from object matrix.
        keys = self.sim.cylinderMatrix.getAll( )
        cylinders = [ key[0].shellList[0] for key in keys ]
        return self.processCylinders( cylinders )


    def getCuboidalSurfaceData( self ):
        posList, typeList, tensorList = [], [], []

        try:
            boxes = self.sim.cuboidalSurfaces
        except:
            # Add dummy box to stop tensorGlyph from complaining.
            boxes = [ DummyBox() ] 

        type = 1
        for box in boxes:
            tensor = numpy.concatenate((box.vectorX, box.vectorY, box.vectorZ))
            self.appendLists(posList, box.origin, typeList, type, 
                    tensorList=tensorList, tensor=tensor) 

        return self.vtk_writer.createDoc(posList, colors=typeList, tensors=tensorList )


    def getCylindricalSurfaceData( self ):
        cylinders = [ surface for surface in self.sim.surfaceList if isinstance(surface, Cylinder) ]
        return self.processCylinders( cylinders )


    def getPlanarSurfaceData( self ):
        posList, typeList, tensorList = [], [], []

        boxes = [ surface for surface in self.sim.surfaceList if isinstance(surface, Box) ]
        if len(boxes) == 0:
            # Add dummy box to stop tensorGlyph from complaining.
            boxes = [ DummyBox() ] 

        type = 1
        for box in boxes:
            tensor = numpy.concatenate((box.vectorX, box.vectorY, box.vectorZ))
            self.appendLists(posList, box.origin, typeList, type, tensorList=tensorList, 
                    tensor=tensor) 

        return self.vtk_writer.createDoc(posList, colors=typeList, tensors=tensorList )


    def processCylinders( self, cylinders=[] ):
        posList, typeList, tensorList = [], [], []

        if len(cylinders) == 0:
            # Add dummy cylinder to stop tensorGlyph from complaining.
            cylinders = [ DummyCylinder() ] 

        for cylinder in cylinders:
            radius = cylinder.radius
            type = 1
            orientation = cylinder.unitZ
            size = cylinder.size

            # Construct tensor. Use tensor glyph plugin from:
            # http://www.paraview.org/pipermail/paraview/2009-March/011256.html
            # Unset Extract eigenvalues.

            # Select basis vector in which orientation is smallest.
            _, basisVector = min( zip(abs(orientation), [[1,0,0], [0,1,0], [0,0,1]]) )
            # Find 2 vectors perpendicular to orientation.
            perpendicular1 = numpy.cross( orientation, basisVector )
            perpendicular2 = numpy.cross( orientation, perpendicular1 )
            # A 'tensor' is represented as an array of 9 values.
            # Stupid Paraview wants  a normal vector to the cylinder to orient  
            # it. So orientation and perpendicular1 swapped.
            tensor = numpy.concatenate((perpendicular1*radius, 
                orientation*size, perpendicular2*radius))

            self.appendLists( posList, cylinder.origin, typeList, type, tensorList=tensorList, tensor=tensor ) 

        return self.vtk_writer.createDoc(posList, colors=typeList, tensors=tensorList )


    # Helper.
    def appendLists(self, posList, pos, typeList, type, radiusList=[], radius=None, 
            tensorList=[], tensor=None):
        factor = 1
        # Convert all lengths to nanometers.
        #factor = 1e8
        posList.append( pos * factor )
        typeList.append( type )

        # Multiply radii and tensors by 2 because Paraview sets radius to 0.5 
        # by default, and it wants full lengths for cylinders and we are 
        # storing half lengths.
        if radius != None:
            radiusList.append( radius * 2 * factor )

        if tensor != None:
            tensor = tensor * 2 * factor
            tensorList.append( tensor )

