#include "pybind.hpp"

#include "stimulus.hpp"

void pybind_stimulus(py::module &m) {
    
    py::class_<StimulusOccupancy, std::shared_ptr<StimulusOccupancy>>(m, "Stimulus",
    R"pbdoc(
        Stimulus Occupancy class.
        
        Parameters
        ----------
        space : Space
            Description of stimulus space.
        grid : Grid
            Evaluation grid for stimulus space.
        stimulus_duration : float
            Duration (in seconds) of single stimulus.
        compression : float
            Compression threshold.
            
    )pbdoc")
    .def( py::init<Space &, Grid &, double, value>(), py::arg("space"), py::arg("grid"), py::arg("stimulus_duration")=1., py::arg("compression")=1. )
    
    .def_property_readonly("compression", &StimulusOccupancy::compression,
    R"pbdoc(Threshold for compression when merging new stimuli into distribution.)pbdoc")
    
    .def_property_readonly("stimulus_duration", &StimulusOccupancy::stimulus_duration,
    R"pbdoc(Duration (in seconds) of single stimulus.)pbdoc")
    
    .def_property_readonly("stimulus_time", &StimulusOccupancy::stimulus_time,
    R"pbdoc(Total stimulus presentation time.)pbdoc")
    
    .def_property_readonly("ndim", &StimulusOccupancy::ndim,
    R"pbdoc(Dimensionality of stimulus space.)pbdoc")
    
    .def_property("random_insertion", &StimulusOccupancy::random_insertion, &StimulusOccupancy::set_random_insertion,
    R"pbdoc(Whether new stimuli will be merged into the distribution in randomized order.)pbdoc")
    
    .def_property_readonly("space", py::cpp_function(&StimulusOccupancy::space, py::return_value_policy::reference_internal),
    R"pbdoc(Stimulus space.)pbdoc")
    
    .def("to_yaml", [](StimulusOccupancy &m)->std::string {
        YAML::Emitter out;
        YAML::Node node = m.to_yaml();
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;
    },
    R"pbdoc(
        Represent stimulus occupancy as YAML.
        
        Returns
        -------
        string
        
    )pbdoc")
    
    .def("save_to_yaml", [](const StimulusOccupancy& obj, std::string path) { obj.save_to_yaml( path ); }, py::arg("path"),
    R"pbdoc(
        Save stimulus occupancy to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        
    )pbdoc")
    
    .def("save_to_hdf5", &StimulusOccupancy::save_to_hdf5,
    py::arg("filename"), py::arg("flags")=19, py::arg("path")="")
    
    .def_static("load_from_hdf5", [](std::string filename, std::string path) { return std::shared_ptr<StimulusOccupancy>(StimulusOccupancy::load_from_hdf5(filename,path)); },
    py::arg("filename"), py::arg("path")="",
    R"pbdoc(
        Load stimulus occupancy from hdf5 file.
        
        Parameters
        ----------
        path : string
            path to hdf5 file
        
        Returns
        -------
        StimulusOccupancy
        
    )pbdoc" )
    
    .def("add_stimuli", [](StimulusOccupancy & obj, py::array_t<value, py::array::c_style | py::array::forcecast> stimuli, unsigned int repetitions) {
        
        unsigned int ndim = obj.ndim();
        unsigned int nsamples;
        
        auto buf = stimuli.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        obj.add_stimulus( (value *) buf.ptr, nsamples, repetitions );
        
        }, 
    py::arg("stimuli"), py::arg("repetitions")=1,
    R"pbdoc(
        Merge new stimuli into distribution.
        
        Parameters
        ----------
        stimuli : (n,ndim) array
            Array of stimulus values.
        repetitions : int
            The number of repetitions for the stimuli.
        
    )pbdoc")
    
    .def("occupancy", [](StimulusOccupancy & obj)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(obj.grid().ndim(), sizeof(value));
        for (int k=obj.grid().ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * obj.grid().shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            obj.grid().ndim(),
            obj.grid().shape(),
            strides
        ));
        
        auto result_buf = result.request();
        
        //std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples, 0. );
        
        obj.occupancy( (value *) result_buf.ptr );
        
        return result;
        
        },
    R"pbdoc(
        Evaluate stimulus occupancy on grid.
        
        Returns
        -------
        nd array
        
    )pbdoc")
    
    .def("logp", [](StimulusOccupancy & obj)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(obj.grid().ndim(), sizeof(value));
        for (int k=obj.grid().ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * obj.grid().shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            obj.grid().ndim(),
            obj.grid().shape(),
            strides
        ));
        
        auto result_buf = result.request();
        
        //std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples, 0. );
        
        obj.logp( (value *) result_buf.ptr );
        
        return result;
        
        },
    R"pbdoc(
        Evaluate log probability of stimulus distribution on grid.
        
        Returns
        -------
        nd array
        
    )pbdoc")
    
    .def("prob", [](StimulusOccupancy & obj)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(obj.grid().ndim(), sizeof(value));
        for (int k=obj.grid().ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * obj.grid().shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            obj.grid().ndim(),
            obj.grid().shape(),
            strides
        ));
        
        auto result_buf = result.request();
        
        //std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples, 0. );
        
        obj.prob( (value *) result_buf.ptr );
        
        return result;
        
        },
    R"pbdoc(
        Evaluate probability of stimulus distribution on grid.
        
        Returns
        -------
        nd array
        
    )pbdoc");

}
