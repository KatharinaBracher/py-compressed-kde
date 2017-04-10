#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

template <class T>
std::vector<T> numpy_array_to_vector( py::array_t<T> & array ) {
    auto buf = array.request();
    std::vector<T> vec( buf.size );
    std::copy( (T*) buf.ptr, ((T*) buf.ptr) + buf.size, vec.begin() );
    return vec;
}
