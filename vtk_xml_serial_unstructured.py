#!/usr/bin/env python
import xml.dom.minidom
#import xml.dom.ext # python 2.5 and later


class VTK_XML_Serial_Unstructured:
    def __init__( self ):
        pass


    def createDoc( self, ( posList, radii, colors, tensors ) ):
        # Document and root element
        doc = xml.dom.minidom.Document()
        root_element = doc.createElementNS( "VTK", "VTKFile" )
        root_element.setAttribute( "type", "UnstructuredGrid" )
        root_element.setAttribute( "version", "0.1" )
        root_element.setAttribute( "byte_order", "LittleEndian" )
        doc.appendChild( root_element )

        # Unstructured grid element
        unstructuredGrid = doc.createElementNS( "VTK", "UnstructuredGrid" )
        root_element.appendChild( unstructuredGrid )

        # Piece 0 (only one)
        # The "Piece" elements are meant for multiple pieces of *geometry*.  
        # They are meant for streaming computation to reduce memory usage.  
        # All the pieces have to have the same set of data arrays.
        # So we can not use that to group particles, spheres and cylinders 
        # into 1 file. Use .pvd file using parts for that.
        piece = doc.createElementNS( "VTK", "Piece" )
        piece.setAttribute( "NumberOfPoints", str( len( posList ) ) )
        piece.setAttribute( "NumberOfCells", "0" )
        unstructuredGrid.appendChild( piece )


        # Points.
        points = doc.createElementNS( "VTK", "Points" )
        piece.appendChild( points )

        # Point location data.
        point_coords = doc.createElementNS( "VTK", "DataArray" )
        point_coords.setAttribute( "type", "Float32" )
        point_coords.setAttribute( "format", "ascii" )
        point_coords.setAttribute( "NumberOfComponents", "3" )
        points.appendChild( point_coords )

        string = str()
        for pos in posList:
            string = string + repr( pos[0] ) + ' ' + repr( pos[1] ) + \
                     ' ' + repr( pos[2] ) + ' '
        point_coords_data = doc.createTextNode( string )
        point_coords.appendChild( point_coords_data )


        # Cells.
        # Don't remove.
        cells = doc.createElementNS( "VTK", "Cells" )
        piece.appendChild( cells )

        # Cell locations.
        cell_connectivity = doc.createElementNS( "VTK", "DataArray" )
        cell_connectivity.setAttribute( "type", "Int32" )
        cell_connectivity.setAttribute( "Name", "connectivity" )
        cell_connectivity.setAttribute( "format", "ascii" )        
        cells.appendChild( cell_connectivity )

        # Cell location data.
        connectivity = doc.createTextNode( "0" )
        cell_connectivity.appendChild( connectivity )

        cell_offsets = doc.createElementNS( "VTK", "DataArray" )
        cell_offsets.setAttribute( "type", "Int32" )
        cell_offsets.setAttribute( "Name", "offsets" )
        cell_offsets.setAttribute( "format", "ascii" )                
        cells.appendChild( cell_offsets )
        offsets = doc.createTextNode( "0" )
        cell_offsets.appendChild( offsets )

        cell_types = doc.createElementNS( "VTK", "DataArray" )
        cell_types.setAttribute( "type", "UInt8" )
        cell_types.setAttribute( "Name", "types" )
        cell_types.setAttribute( "format", "ascii" )                
        cells.appendChild( cell_types )
        types = doc.createTextNode( "1" )
        cell_types.appendChild( types )


        # Point data.
        point_data = doc.createElementNS( "VTK", "PointData" )
        piece.appendChild( point_data )

        # Radii.
        if len( radii ) > 0:
            radiiNode = doc.createElementNS( "VTK", "DataArray" )
            radiiNode.setAttribute( "Name", "radii" )
            radiiNode.setAttribute( "type", "Float32" )
            radiiNode.setAttribute( "format", "ascii" )
            point_data.appendChild( radiiNode )

            string = str()
            for radius in radii:
                string = string + repr( radius ) + ' '
            radiiData = doc.createTextNode( string )
            radiiNode.appendChild( radiiData )

        # Colors.
        if len( colors ) > 0:
            colorNode = doc.createElementNS( "VTK", "DataArray" )
            colorNode.setAttribute( "Name", "colors" )
            colorNode.setAttribute( "NumberOfComponents", "1" )
            colorNode.setAttribute( "type", "Float32" )
            colorNode.setAttribute( "format", "ascii" )
            point_data.appendChild( colorNode )

            string = str()
            for color in colors:
                string = string + repr( color ) + ' '
            colorData = doc.createTextNode( string )
            colorNode.appendChild( colorData )

        # Tensors.
        if len( tensors ) > 0:

            # Hack to make VTK understand I want to color the TensorGlyphs.

            # I think there is a bug actually in VTK somewhere, when using 
            # tensorGlyph.xml to make a vtkTensorGlyph object with both a 
            # Tensor array and a Scalar array. When you select a Tensor or a 
            # Scalar array from the dropdown menu and click 'apply', 
            # SetInputArrayToProcess is called with idx = 0 both times. This 
            # is wrong.
            # 1. First element Tensor array gets overwritten.
            # 2. Scalar value is never written ( which is accessed using 
            # GetInputArrayToProcess with an idx of 1.
            
            # The workaround here uses an additional Vector array ( 
            # vtkTensorGlyph doesn't have a vector array ), which when it gets 
            # updated actually results in SetInputArrayToProcess to be called 
            # with idx = 1. So by also supplying Paraview with a vector for 
            # each color ( just 3 times the color int ), the Tensor array 
            # doesn't get overwritten and GetInputArrayToProcess with idx is 
            # also happy and updates Scalars ( ! ).

            # Tried to find the root cause using Gdb but failed.

            #Important files:
            #    ParaView3/VTK/Graphics/vtkTensorGlyph.cxx
            #    ParaView3/VTK/Filtering/vtkAlgorithm.cxx

            if len( colors ) > 0:
                colorNode = doc.createElementNS( "VTK", "DataArray" )
                colorNode.setAttribute( "Name", "colors_( really_a_scalar )" ) 
                colorNode.setAttribute( "NumberOfComponents", "3" )
                colorNode.setAttribute( "type", "Float32" )
                colorNode.setAttribute( "format", "ascii" )
                point_data.appendChild( colorNode )

                string = str()
                for color in colors:
                    string += repr( color ) + ' ' + repr( color ) + ' ' + \
                              repr( color ) + ' '
                colorData = doc.createTextNode( string )
                colorNode.appendChild( colorData )

            tensorNode = doc.createElementNS( "VTK", "DataArray" )
            tensorNode.setAttribute( "Name", "tensors" )
            tensorNode.setAttribute( "NumberOfComponents", "9" )
            tensorNode.setAttribute( "type", "Float32" )
            tensorNode.setAttribute( "format", "ascii" )
            point_data.appendChild( tensorNode )

            string = str()
            for tensor in tensors:
                for value in tensor:
                    # A 'tensor' is represented as a list of 9 values.
                    string = string + repr( value ) + ' '
            tensorData = doc.createTextNode( string )
            tensorNode.appendChild( tensorData )


        # Cell data ( Dummy ).
        cell_data = doc.createElementNS( "VTK", "CellData" )
        piece.appendChild( cell_data )

        return doc


    def writeDoc( self, doc, fileName ):
        # Write to file and exit.
        outFile = open( fileName, 'w' )
        # xml.dom.ext.PrettyPrint( doc, file )
        doc.writexml( outFile, newl='\n' )
        outFile.close()


    def writePVD( self, file, fileList ):
        outFile = open( file, 'w' )

        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS( "VTK", "VTKFile" )
        pvd_root.setAttribute( "type", "Collection" )
        pvd_root.setAttribute( "version", "0.1" )
        pvd_root.setAttribute( "byte_order", "LittleEndian" )
        pvd.appendChild( pvd_root )

        collection = pvd.createElementNS( "VTK", "Collection" )
        pvd_root.appendChild( collection )

        for type, fileName, index, time in fileList:
            dataSet = pvd.createElementNS( "VTK", "DataSet" )
            if time:
                # Problem with adding real time is that TimestepValues is not 
                # updated in Proxy group="misc" in .pvsm file after a reload.
                #time = str( time )
                time = str( index )
                dataSet.setAttribute( "timestep", time )
            dataSet.setAttribute( "group", "" )
            dataSet.setAttribute( "part", type ) # Use ExtractBlock.
            dataSet.setAttribute( "file", fileName )
            collection.appendChild( dataSet )

        outFile = open( file, 'w' )
        pvd.writexml( outFile, newl='\n' )
        outFile.close()

