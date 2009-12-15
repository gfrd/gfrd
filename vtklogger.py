
import os
import shutil
#import re
#import logging
import numpy
from vtk_xml_serial_unstructured import *
from shape import *
from bd import BDSimulator

INF = numpy.inf

class VTKLogger:
    """Logger that can be used to visualize data with Kitware ParaView.


    * Setup. Specify bufferSize to only write last 'bufferSize' simulation 
      steps to file:
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
    To visualize the cylinders a workaround using tensors is used, as 
    explained here: 
        http://www.paraview.org/pipermail/paraview/2009-March/011256.html.
    The mentioned tensorGlyph.xml should be supplied with this package.

    As explainded in the above link, to give cylinders the right color, there 
    are 2 options.

    1. Build Paraview from source, after adding in file
    Paraview3/VTK/Graphics/vtkTensorGlyph.cxx after
        newScalars = vtkFloatArray::New();
    a new line:
        newScalars->SetName( "colors" );
    2. Do that Python script thing.

    Another hack was needed to get the coloring to work. This makes VTKLogger 
    write a vector instead of a scalar value for each color to the .vtu files. 
    But you shouldn't have to worry about that, just select 'colors' to color 
    the Glyph.

    When doing Brownian Dynamics, don't show shells.

    extraParticleStep=True means that for each timestep an extra step is  
    recorded where only the active particle has been updated (it's shell
    stays unchanged).

    """

    def __init__( self, sim, name, bufferSize=None, brownian=False, 
                  extraParticleStep=True ):
        self.sim = sim
        self.brownian = brownian
        self.extraParticleStep = extraParticleStep

        self.vtk_writer = VTK_XML_Serial_Unstructured()
        self.bufferSize = bufferSize
        self.buffer = []

        # Filename stuff.
        assert len( name ) > 0
        self.name = name

        self.fileNameNumberOfSteps = self.name + '/numberOfSteps.dat'
        if os.path.exists( self.fileNameNumberOfSteps ):
            # Simulation with the same name has been run before. Retrieve 
            # number of steps.
            outFile = open( self.fileNameNumberOfSteps, 'r' )
            self.numberOfStepsPreviousRun = int( outFile.read() )
        else:
            self.numberOfStepsPreviousRun = 0

        filesDirectory = self.name + '/files'
        if not os.path.exists( filesDirectory ):
            os.makedirs( filesDirectory ) 

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
        else:
            spheres, cylinders = self.getDummyData(), self.getDummyCylinders()

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


    def writelog( self, time, index, ( particles, spheres, cylinders )  ):
        if index == 0:
            self.previousSpheres = spheres
            self.previousCylinders = cylinders

        if not self.brownian and self.extraParticleStep:
            # Show step where only particles have been updated.
            index *= 2;
            self.makeSnapshot( 'particles', particles, index, time )
            self.makeSnapshot( 'spheres', self.previousSpheres, index, time )
            self.makeSnapshot( 'cylinders', self.previousCylinders, index, time )

            time += self.deltaT
            index += 1

        self.makeSnapshot( 'particles', particles, index, time )
        self.makeSnapshot( 'spheres', spheres, index, time )
        self.makeSnapshot( 'cylinders', cylinders, index, time )

        self.previousSpheres = spheres
        self.previousCylinders = cylinders


    def makeSnapshot( self, type, data, index='', time=None ):
        """Write data to file.

        """
        doc = self.vtk_writer.createDoc( data )
        fileName = 'files/' + type + str( index ) + '.vtu'
        self.vtk_writer.writeDoc( doc, self.name + '/' + fileName )

        # Store filename and time in fileList, used by vtk_wrtie.writePVD().
        if time:
            self.fileList.append( ( type, fileName, index, time ) )
        else:
            self.staticList.append( ( type, fileName, None, None ) )


    def stop( self ):
        self.log()

        # Write contents of buffer.
        for index, entry in enumerate( self.buffer ):
            if index % 10 == 0:
                print 'vtklogger writing step %s from buffer' % index
            self.writelog( entry[0], index, entry[2:] )

        # Surfaces don't move.
        self.makeSnapshot( 'cylindricalSurfaces', 
                           self.getCylindricalSurfaceData() )
        self.makeSnapshot( 'planarSurfaces', self.getPlanarSurfaceData() )
        self.makeSnapshot( 'cuboidalSurfaces', self.getCuboidalSurfaceData() )

        # Fill with dummy files up to number of steps of previous run. This 
        # should have Paraview complain less.
        dummyData = ( self.getDummyData(), self.getDummyData(),
                      self.getDummyCylinders() )

        # Overwrite, so transition to dummy data is more smooth.
        self.previousSpheres = self.getDummyData()
        self.previousCylinders = self.getDummyCylinders()

        if self.numberOfStepsPreviousRun > self.i:
            print ( 'vtklogger writing dummy files for step %d to %d' % 
                    ( self.i, self.numberOfStepsPreviousRun ) )
            for index in range( self.i, self.numberOfStepsPreviousRun ):
                if index % 10 == 0:
                    print 'vtklogger writing step %s from buffer' % index
                self.writelog( index, index, dummyData )

        # Write number of steps to file. Needed for next run, maybe.
        numberOfSteps = self.i - 1
        if self.bufferSize:
            numberOfSteps = min( self.bufferSize, numberOfSteps )
        outFile = open( self.fileNameNumberOfSteps, 'w' )
        outFile.write( str( numberOfSteps ) )
        outFile.close()

        # Finally, write PVD files.
        self.vtk_writer.writePVD( self.name + '/' + 'files.pvd', 
                                  self.fileList )

        self.vtk_writer.writePVD( self.name + '/' + 'static.pvd', 
                                  self.staticList )


    def getDummyData( self ):
        return ( [], [], [], [] )


    def getDummyCylinders( self ):
        return self.processCylinders( [], [] )


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

        numberOfShells = self.sim.scheduler.getSize()
        if numberOfShells > 0:
            topTime = self.sim.scheduler.getTopTime()

        for eventIndex in range( numberOfShells ):
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


    def getCuboidalSurfaceData( self ):
        posList, typeList, tensorList = [], [], []

        try:
            boxes = [ self.sim.defaultSurface ]
        except:
            # Add dummy box to stop tensorGlyph from complaining.
            boxes = [ Box( [ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0 ], [ 0, 0, 1 ], 
                           1e-20, 1e-20, 1e-20 ) ]

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
            boxes = [ Box( [ 0, 0, 0], [ 1, 0, 0], [ 0, 1, 0 ], [ 0, 0, 1 ], 
                           1e-20, 1e-20, 1e-20 ) ]

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
            cylinders = [ Cylinder( [ 0, 0, 0 ], 1e-20, [ 0, 0, 1 ], 1e-20 ) ]

        for cylinder in cylinders:
            radius = cylinder.radius
            orientation = cylinder.unitZ
            size = cylinder.size

            # Construct tensor. Use tensor glyph plugin from:
            # http://www.paraview.org/pipermail/paraview/2009-March/011256.html
            # Unset Extract eigenvalues.

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


    def appendLists( self, posList, pos, typeList=[], type=None, 
                     radiusList=[], radius=None, tensorList=[], tensor=None ):
        """Helper.

        """
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

