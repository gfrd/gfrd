#include <iostream>
#include "position.hpp"

#include "peer/utils.hpp"

#include "peer/ObjectContainer.hpp"


typedef position<double> position_type;

struct position_to_ndarray_converter
{
    typedef position_type native_type;
    
    static PyObject* convert( const native_type& p )
    {
        static const npy_intp dims[1] = { native_type::size() };
        void* data( PyDataMem_NEW( native_type::size() * sizeof( double ) ) );
        memcpy( data, static_cast<const void*>( p.data() ),
                native_type::size() * sizeof( double ) );
        PyObject* array( PyArray_New( &PyArray_Type, 1, 
                                      const_cast<npy_intp*>( dims ),
                                      peer::util::get_numpy_typecode<double>
                                      ::value, NULL,
                                      data, 0, NPY_CARRAY, NULL ) );
        reinterpret_cast<PyArrayObject*>( array )->flags |= NPY_OWNDATA;
        return array;
    }
};

struct ndarray_to_position_converter
{
    typedef position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PyArray_Check(ptr))
        {
            return NULL;
        }
        
        PyObject* retval(PyArray_CastToType(
                             reinterpret_cast<PyArrayObject*>(ptr),
                             PyArray_DescrFromType(
                                 peer::util::get_numpy_typecode<double>::value ),
                                 0) );
        if (!retval)
        {
            return NULL;
        }
        
        if (PyObject_Size(retval) != 3)
        {
            boost::python::decref(retval);
            return NULL;
        }
        
        return retval;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        PyArrayObject* array_obj = static_cast<PyArrayObject*>(
            data->stage1.convertible);
        data->stage1.convertible = new(data->storage.bytes) native_type(
            reinterpret_cast<double*>(PyArray_DATA(array_obj)));
        boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
    }
};

struct seq_to_position_converter
{
    typedef position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PySequence_Check(ptr))
        {
            return NULL;
        }
        
        if (PySequence_Size(ptr) != 3)
        {
            return NULL;
        }
        
        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        data->stage1.convertible = new(data->storage.bytes) native_type(
            PyFloat_AsDouble(PySequence_GetItem(ptr, 0)),
            PyFloat_AsDouble(PySequence_GetItem(ptr, 1)),
            PyFloat_AsDouble(PySequence_GetItem(ptr, 2)));
    }
};



BOOST_PYTHON_MODULE(object_matrix)
{
    using namespace boost::python;

    to_python_converter<position_type,
        position_to_ndarray_converter>();
    peer::util::to_native_converter<position_type,
        ndarray_to_position_converter>();
    peer::util::to_native_converter<position_type,
        seq_to_position_converter>();

    //register_ptr_to_python<boost::shared_ptr<position_type> >();


    import_array();
#if OBJECTMATRIX_USE_ITERATOR
    peer::util::register_stop_iteration_exc_translator();
#endif
    typedef sphere<double> default_mapped_type;
    peer::ObjectContainer<default_mapped_type>::__register_class();
    peer::SphereContainer::__register_class();
    peer::CylinderContainer::__register_class();
    peer::BoxContainer::__register_class();
}

