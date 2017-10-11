#include "pybind.hpp"

#include "decoder.hpp"

void pybind_decoder(py::module &m) {
    
    py::class_<Decoder>(m, "Decoder",
    R"pbdoc(
        Decoder class.
        
        Construction of a Decoder requires one or more sets of PoissonLikelihood
        objects and optional priors. There are two ways to construct a Decoder,
        depending on whether you would like to decode over multiple stimulus spaces
        or not.
        
        | Decoder( likelihoods, prior )
        | Decoder( likelihoods, priors )
        
        In the first case, one passes a list of likelihood objects that all share the
        same stimulus space and (optionally) an array with prior probabilities for
        the grid points in the stimulus space.
        
        In the second case, one passes a nested list of likelihoods, where the outer
        list represents the data sources (e.g. individual tetrodes) and the inner lists
        represent multiple stimulus spaces one would like to decode over. The (optional)
        prior probabilities for each stimulus space are provided as a list of arrays.
        
        Parameters
        ----------
        likelihoods : list or nested list of PoissonLikelihood objects
        prior(s) : array or list of arrays
        
    )pbdoc")
    
    .def( py::init<std::vector<std::shared_ptr<PoissonLikelihood>> &, std::vector<value> &>(), py::arg("likelihoods"), py::arg("prior") )
    .def( py::init<std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> &, std::vector<std::vector<value>> &>(), py::arg("likelihoods"), py::arg("priors") )
    
    .def_property_readonly("nsources", &Decoder::nsources,
    R"pbdoc(Number of sources (likelihoods).)pbdoc")
    
    .def_property_readonly("is_union", &Decoder::is_union,
    R"pbdoc(Whether decoding is performed over union of stimulus spaces.)pbdoc")
    
    .def_property_readonly("n_union", &Decoder::n_union,
    R"pbdoc(Number of stimulus spaces over which decoding is performed.)pbdoc")
    
    .def_property_readonly("grid_shapes", &Decoder::grid_shapes,
    R"pbdoc(Grid shape for each stimulus space.)pbdoc")
    
    .def("grid_shape", &Decoder::grid_shape, py::arg("index")=0,
    R"pbdoc(
        Grid size.
        
        Parameters
        ----------
        index : int
            Index of stimulus space in union (zero-based).
        
        Returns
        -------
        grid shape
            
    )pbdoc")
    
    .def_property_readonly("grid_sizes", &Decoder::grid_sizes,
    R"pbdoc(Grid size for each stimulus space.)pbdoc")
    
    .def("grid_size", &Decoder::grid_size, py::arg("index")=0,
    R"pbdoc(
        Grid size.
        
        Parameters
        ----------
        index : int
            Index of stimulus space in union (zero-based).
        
        Returns
        -------
        int : grid size
            
    )pbdoc")
    
    .def_property_readonly("nenabled_sources", &Decoder::nenabled_sources,
    R"pbdoc(Number of sources that are used for decoding.)pbdoc")
    
    .def_property_readonly("enabled_sources", &Decoder::enabled_sources,
    R"pbdoc(Enabled state for all sources.)pbdoc")
    
    .def("grid", &Decoder::grid, py::return_value_policy::reference_internal, py::arg("index")=0,
    R"pbdoc(
        Get stimulus space grid.
        
        Parameters
        ----------
        index : int
            Index of stimulus space in union (zero-based).
        
        Returns
        -------
        Grid
            
    )pbdoc")
    
    .def("stimulus", &Decoder::stimulus, py::arg("index")=0,
    R"pbdoc(
        Get stimulus space.
        
        Parameters
        ----------
        index : int
            Index of stimulus space in union (zero-based).
        
        Returns
        -------
        Stimulus
            
    )pbdoc")
    
    .def("likelihood", &Decoder::likelihood, py::arg("source"), py::arg("index")=0,
    R"pbdoc(
        Get likelihood.
        
        Parameters
        ----------
        source : int
            Index of source likelihood (zero-based).
        index : int
            Index of stimulus space in union (zero-based).
        
        Returns
        -------
        PoissonLikelihood
            
    )pbdoc")
    
    .def("enable_source", &Decoder::enable_source, py::arg("source"),
    R"pbdoc(
        Enable source.
        
        Parameters
        ----------
        source : int
            Index of source (zero-based).
            
    )pbdoc")
    
    .def("disable_source", &Decoder::disable_source, py::arg("source"),
    R"pbdoc(
        Disable source.
        
        Parameters
        ----------
        source : int
            Index of source (zero-based).
            
    )pbdoc")
    
    .def("enable_all_sources", &Decoder::enable_all_sources,
    R"pbdoc(Enable all sources.)pbdoc")
    
    .def("enable_one_source", &Decoder::enable_one_source, py::arg("source"),
    R"pbdoc(
        Enable single source and disable all others.
        
        Parameters
        ----------
        source : int
            Index of source (zero-based).
            
    )pbdoc")
    
    .def("enable_sources", &Decoder::enable_sources, py::arg("state"),
    R"pbdoc(
        Set enable state of sources.
        
        Parameters
        ----------
        state : list or 1d array
            For each source the enabled state (True/False).
            
    )pbdoc")
    
    .def("save_to_hdf5", &Decoder::save_to_hdf5,
    py::arg("filename"), py::arg("flags")=19, py::arg("path")="")
    
    .def_static("load_from_hdf5", &Decoder::load_from_hdf5,
    py::arg("filename"), py::arg("path")="",
    R"pbdoc(
        Load poisson likelihood from hdf5 file.
        
        Parameters
        ----------
        path : string
            path to hdf5 file
        
        Returns
        -------
        PoissonLikelihood
        
    )pbdoc" )
    
    .def("decode", [](Decoder & obj, std::vector<py::array_t<value, py::array::c_style | py::array::forcecast>> events, value delta_t, bool normalize)->std::vector<py::array_t<value>> {
        
        std::vector<py::array_t<value>> out;
        std::vector<value*> out_ptr;
        
        // for each union
        for (unsigned int k = 0 ; k<obj.n_union(); ++k) {
            
            std::vector<long unsigned int> strides(obj.grid_shape(k).size(), sizeof(value));
            for (int s=obj.grid_shape(k).size()-2; s>=0; --s) {
                strides[s] = strides[s+1] * obj.grid_shape(k)[s+1];
            }
            
            // construct array buffer
            auto result = py::array( py::buffer_info(
                nullptr,
                sizeof(value),
                py::format_descriptor<value>::value,
                1,
                obj.grid_shape(k),
                strides
            ));
            
            auto result_buf = result.request();
            
            std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + obj.grid_size(k), 0. );
            
            out.push_back( result );
            out_ptr.push_back( (value *) result_buf.ptr );
        }
        
        // construct vector of data pointers
        std::vector<value*> events_data;
        std::vector<unsigned int> events_n;
        
        for (auto & item : events) {
            auto buf = item.request();
            events_data.push_back( (value*) buf.ptr );
            events_n.push_back( buf.size );
        }
        
        obj.decode( events_data, events_n, delta_t, out_ptr, normalize );
        
        return out;
        
    }, py::arg("events"), py::arg("delta"), py::arg("normalize")=true,
    R"pbdoc(
        Compute posterior probability distribution.
        
        Parameters
        ----------
        events : list of (n,ndim) arrays
            A list with for each source the observed event data.
        delta : float
            Time duration over which events were observed.
        normalize : bool
            Normalize posterior distribution such that is sums to one.
        
        Returns
        -------
        list with posterior distribution for each of the union-ed stimulus spaces.
        
    )pbdoc")
    .def("decode_single", [](Decoder & obj, std::vector<py::array_t<value, py::array::c_style | py::array::forcecast>> events, value delta_t, unsigned int index, bool normalize)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(obj.grid_shape(index).size(), sizeof(value));
        for (int s=obj.grid_shape(index).size()-2; s>=0; --s) {
            strides[s] = strides[s+1] * obj.grid_shape(index)[s+1];
        }
        
        // construct array buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            1,
            obj.grid_shape(index),
            strides
        ));
        
        auto result_buf = result.request();
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + obj.grid_size(index), 0. );
        
        // construct vector of data pointers
        std::vector<value*> events_data;
        std::vector<unsigned int> events_n;
        
        for (auto & item : events) {
            auto buf = item.request();
            events_data.push_back( (value*) buf.ptr );
            events_n.push_back( buf.size );
        }
        
        obj.decode( events_data, events_n, delta_t, (value*) result_buf.ptr, index, normalize );
        
        return result;
        
    }, py::arg("events"), py::arg("delta"), py::arg("index")=0, py::arg("normalize")=true,
    R"pbdoc(
        Compute posterior probability distribution for single stimulus space.
        
        Parameters
        ----------
        events : list of (n,ndim) arrays
            A list with for each source the observed event data.
        delta : float
            Time duration over which events were observed.
        index : int
            Index of stimulus space in union that is target of decoding.
        normalize : bool
            Normalize posterior distribution such that is sums to one.
        
        Returns
        -------
        posterior distribution for selected stimulus space.
        
    )pbdoc");

}
