
import os
#import re
#import logging
import numpy
from vtk_xml_serial_unstructured import *
from shape import *
from bd import BDSimulator

INF = numpy.inf


'''
Logger that can be used to visualize data with Kitware ParaView.


* Setup. Specify bufferSize to only write last 'bufferSize' simulation steps 
  to file:
    vtklogger = VTKLogger( sim=s, name='run' )
    or
    vtklogger = VTKLogger( sim=s, name='run', bufferSize=100 )

* Each step, to write .vtk files:
    vtklogger.log()

* Finalize, to write .pvd file:
    vtklogger.stop()


Inner workings:
    log():
        ( get.......Data = createDoc() + addPiece() ) + writeDoc()
    stop():
        writePvd( name.pvd )


=== Cylinders
To visualize the cylinders a workaround using tensors is used, as explained 
here: http://www.paraview.org/pipermail/paraview/2009-March/011256.html. The 
mentioned tensorGlyph.xml should be supplied with this package.

As explainded in the above link, to give cylinders the right color, there are 
2 options.

1. Build Paraview from source, after adding in file
Paraview3/VTK/Graphics/vtkTensorGlyph.cxx after
    newScalars = vtkFloatArray::New();
a new line:
    newScalars->SetName( "colors" );
2. Do that Python script thing.

Another hack was needed to get the coloring to work. This makes VTKLogger 
write a vector instead of a scalar value for each color to the .vtu files. But 
you shouldn't have to worry about that, just select 'colors' to color the 
Glyph.
'''
class VTKLogger:
    def __init__( self, sim, name, bufferSize=None, brownian=False ):
        self.sim = sim
        self.brownian = brownian

        self.vtk_writer = VTK_XML_Serial_Unstructured()
        self.bufferSize = bufferSize
        self.buffer = []

        # Filename stuff.
        self.name = name
        if not os.path.exists( 'data/' + self.name + '/files' ):
            os.makedirs( 'data/' + self.name + '/files' )

        self.fileList = []
        self.staticList = []

        self.i = 0          # Step counter.
        self.deltaT = 1e-11 # Needed for hack. Note: don't make too small, 
                            # should be relative to max time.
        self.lastTime = 0   # Needed for hack.


    def log( self ):
        time = self.sim.t
        if ( abs( time - self.lastTime ) < 1e-9 ) and not self.brownian:
            # Hack to make Paraview understand this is a different simulator 
            # step but with the same time.
            # 1. During multi global time is not updated.
            # 2. During initialization time is 0.
            # 3. Sometimes shells have mobilityRadius = 0 --> dt = 0.
            # And I want every step recorded, so I can find out what happened 
            # in which step (divide by 2 actually).
            #
            # Now Paraview should perceive a state change.
            # When doing Brownian Dynamics, this is not needed.
            time = self.lastTime + self.deltaT

        # Get data.
        particles = self.getParticleData()
        if not self.brownian:
            spheres, cylinders = self.getShellDataFromScheduler( )
            #cylinders = self.getCylindricalShellData( )

            if self.i == 0:
                self.previousShells = spheres
                self.previousCylinders = cylinders
        else:
            spheres, cylinders = [], []

        # Write to buffer or file.
        if self.bufferSize:
            # Store in buffer, instead of writing to file directly.
            self.buffer.append( ( time, self.i, particles, spheres, 
                                  cylinders ) )

            if self.i >= self.bufferSize:
                # FIFO.
                del self.buffer[0]
        else:
            # Write normal log.
            self.writelog( time, self.i, ( particles, spheres, cylinders ) )

        self.i += 1
        self.lastTime = time


    def writelog( self, time, index, ( particles, spheres, cylinders ) ):
        self.makeSnapshot( 'particles', particles, index, time )
        if not self.brownian:
            self.makeSnapshot( 'spheres', spheres, index, time )
            self.makeSnapshot( 'cylinders', cylinders, index, time )


    # Write data to file.
    def makeSnapshot( self, type, data, index='', time=None ):
        doc = self.vtk_writer.createDoc( data )
        fileName = 'files/' + type + str( index ) + '.vtu'
        self.vtk_writer.writeDoc( doc, 'data/' + self.name + '/' + fileName )

        # Store filename and time in fileList, used by vtk_wrtie.writePVD().
        if time:
            self.fileList.append( ( type, fileName, index, time ) )
        else:
            self.staticList.append( ( type, fileName, None, None ) )


    def stop( self ):
        self.log()
        for newIndex, entry in enumerate( self.buffer ):
            self.writelog( entry[0], newIndex, entry[2:] )
        self.vtk_writer.writePVD( 'data/' + self.name + '/' + 'files.pvd', 
                                  self.fileList )

        # Surfaces don't move.
        self.makeSnapshot( 'cylindricalSurfaces', 
                           self.getCylindricalSurfaceData() )
        self.makeSnapshot( 'planarSurfaces', self.getPlanarSurfaceData() )
        self.makeSnapshot( 'cuboidalSurfaces', self.getCuboidalSurfaceData() )

        self.vtk_writer.writePVD( 'data/' + self.name + '/' + 'static.pvd', 
                                  self.staticList )


    def getParticleData( self ):
        posList, radiusList, typeList = [], [], []

        for speciesIndex, species in enumerate( self.sim.speciesList.values() ):
            for particlePos in species.pool.positions:
                self.appendLists( posList, particlePos, typeList, 
                                  speciesIndex, radiusList, species.radius )

        return ( posList, radiusList, typeList, [] )


    def getShellDataFromScheduler( self ):
        posList, radiusList, typeList = [], [], []
        cylinders, cylinderTypeList = [], []

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
                cylinders.append( object.shellList[0] )
                cylinderTypeList.append( type )
            except:
                # Single or Pair or Multi.
                for shell in object.shellList:
                    self.appendLists( posList, shell.origin, typeList, type, 
                                      radiusList, shell.radius )

        return ( posList, radiusList, typeList, [] ), \
               self.processCylinders( cylinders, cylinderTypeList )


    def getCylindricalShellData( self ):
        # Get data from object matrix.
        keys = self.sim.cylinderMatrix.getAll( )
        cylinders = [ key[0].shellList[0] for key in keys ]
        return self.processCylinders( cylinders )


    def getCuboidalSurfaceData( self ):
        posList, typeList, tensorList = [], [], []

        try:
            boxes = [ self.sim.defaultSurface ]
        except:
            # Add dummy box to stop tensorGlyph from complaining.
            boxes = [ DummyBox() ] 

        for box in boxes:
            tensor = numpy.concatenate( ( box.vectorX, box.vectorY, 
                                          box.vectorZ ) )
            self.appendLists( posList, box.origin, tensorList=tensorList, 
                              tensor=tensor ) 

        return ( posList, [], [], tensorList )


    def getCylindricalSurfaceData( self ):
        # Todo. Make DNA blink when reaction takes place.
        cylinders = [ surface for surface in self.sim.surfaceList if 
                      isinstance( surface, Cylinder ) ]
        return self.processCylinders( cylinders )


    def getPlanarSurfaceData( self ):
        posList, typeList, tensorList = [], [], []

        boxes = [ surface for surface in self.sim.surfaceList if 
                  isinstance( surface, Box ) ]
        if len( boxes ) == 0:
            # Add dummy box to stop tensorGlyph from complaining.
            boxes = [ DummyBox() ] 

        for box in boxes:
            tensor = numpy.concatenate( ( box.vectorX, box.vectorY,
                                          box.vectorZ ) )
            self.appendLists( posList, box.origin, tensorList=tensorList, 
                              tensor=tensor ) 

        return ( posList, [], [], tensorList )


    def processCylinders( self, cylinders=[], typeList=[] ):
        posList, tensorList = [], []

        if len( cylinders ) == 0:
            # Add dummy cylinder to stop tensorGlyph from complaining.
            cylinders = [ DummyCylinder() ] 

        for cylinder in cylinders:
            radius = cylinder.radius
            orientation = cylinder.unitZ
            size = cylinder.size

            '''
            Construct tensor. Use tensor glyph plugin from:
            http://www.paraview.org/pipermail/paraview/2009-March/011256.html
            Unset Extract eigenvalues.
            '''

            # Select basis vector in which orientation is smallest.
            _, basisVector = min( zip( abs( orientation ), 
                                       [ [ 1, 0, 0 ], 
                                         [ 0, 1, 0 ], 
                                         [ 0, 0, 1 ] ] ) )
            # Find 2 vectors perpendicular to orientation.
            perpendicular1 = numpy.cross( orientation, basisVector )
            perpendicular2 = numpy.cross( orientation, perpendicular1 )
            # A 'tensor' is represented as an array of 9 values.
            # Stupid Paraview wants  a normal vector to the cylinder to orient  
            # it. So orientation and perpendicular1 swapped.
            tensor = numpy.concatenate( ( perpendicular1 * radius, 
                orientation * size, perpendicular2 * radius ) )

            self.appendLists( posList, cylinder.origin, tensorList=tensorList, 
                              tensor=tensor ) 

        return ( posList, [], typeList, tensorList )


    # Helper.
    def appendLists( self, posList, pos, typeList=[], type=None, 
                     radiusList=[], radius=None, tensorList=[], tensor=None ):
        factor = 1
        # Convert all lengths to nanometers.
        #factor = 1e8
        posList.append( pos * factor )
        if type != None:
            typeList.append( type )

        # Multiply radii and tensors by 2 because Paraview sets radius to 0.5 
        # by default, and it wants full lengths for cylinders and we are 
        # storing half lengths.
        if radius:
            radiusList.append( radius * 2 * factor )

        if tensor != None:
            tensor = tensor * 2 * factor
            tensorList.append( tensor )


class DummyCylinder( Cylinder ):
    def __init__( self ):
        Cylinder.__init__( self, [ 0, 0, 0 ], 1e-20, [ 0, 0, 1 ], 1e-20 ) 


class DummyBox( Box ):
    def __init__( self ):
        Box.__init__( self, [ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0 ], [ 0, 0, 1 ], 
                      1e-20, 1e-20, 1e-20 ) 

