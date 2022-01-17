#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "datatype_generated.h"

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


template <typename T>
py::tuple pickle_get_state(const T &obj) {
    flatbuffers::FlatBufferBuilder builder(1024);
    auto f = obj.to_flatbuffers(builder);
    builder.Finish(f);

    uint8_t *buf = builder.GetBufferPointer();
    int size = builder.GetSize();

    return py::make_tuple(py::bytes((char*)buf,size));
}

template <typename T, typename S>
std::unique_ptr<T> pickle_set_state(py::tuple t) {
    if (t.size() != 1)
        throw std::runtime_error("Invalid state!");

    char *buffer_pointer;
    ssize_t length;
    PYBIND11_BYTES_AS_STRING_AND_SIZE(
        t[0].cast<py::bytes>().ptr(), &buffer_pointer, &length
    );

    // Get a pointer to the root object inside the buffer.
    auto ptr = flatbuffers::GetRoot<S>((uint8_t*)buffer_pointer);

    /* Create a new C++ instance */
    return T::from_flatbuffers(ptr);
}
