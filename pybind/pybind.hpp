#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

template <class T>
std::vector<T> numpy_array_to_vector( py::array_t<T, py::array::c_style | py::array::forcecast> & array ) {
    auto buf = array.request();
    std::vector<T> vec( buf.size );
    std::copy( (T*) buf.ptr, ((T*) buf.ptr) + buf.size, vec.begin() );
    return vec;
}

// redefine HighFive::File enum for flags
enum Flags : unsigned {
        /// Open flag: Read only access
        ReadOnly = 0x00u,
        /// Open flag: Read Write access
        ReadWrite = 0x01u,
        /// Open flag: Truncate a file if already existing
        Truncate = 0x02u,
        /// Open flag: Open will fail if file already exist
        Excl = 0x04u,
        /// Open flag: Open in debug mode
        Debug = 0x08u,
        /// Open flag: Create non existing file
        Create = 0x10u,
        /// Derived open flag: common write mode (=ReadWrite|Create|Truncate)
        Overwrite = Truncate,
        /// Derived open flag: Opens RW or exclusivelly creates
        OpenOrCreate = ReadWrite | Create
    };