#include "pybind.hpp"

#include "likelihood.hpp"

void pybind_likelihood(py::module &m) {
    
    py::class_<PoissonLikelihood, std::shared_ptr<PoissonLikelihood>>(m, "PoissonLikelihood",
    R"pbdoc(
        Poisson likelihood class.
        
        There are three ways to construct this object, depending on
        whether you have already constructed a Stimulus object and whether events have 
        attributes (e.g. spike amplitude) or not (e.g. when using sorted spikes).
        
        | PoissonLikelihood( event_space, stimulus_space, grid, stimulus_duration, compression )
        | PoissonLikelihood( stimulus_space, grid, stimulus_duration, compression )
        | PoissonLikelihood( event_space, stimulus )
        | PoissonLikelihood( stimulus )
        
        Parameters
        ----------
        event_space : Space object
            Description of event space.
        stimulus_space : Space object
            Description of stimulus space.
        stimulus : Stimulus object
            Stimulus distribution. The stimulus space, grid and compression threshold will be
            determined by the Stimulus object.
        grid : Grid object
            Evaluation grid for stimulus space.
        stimulus_duration : float
            Duration (in seconds) of single stimulus.
        compression : float
            Compression threshold.
        
    )pbdoc")
    .def( py::init<Space &, Grid &, double, value>(), py::arg("stimulus_space"), py::arg("grid"), py::arg("stimulus_duration")=1., py::arg("compression")=1. )
    .def( py::init<Space &, Space &, Grid &, double, value>(), py::arg("event_space"), py::arg("stimulus_space"), py::arg("grid"), py::arg("stimulus_duration")=1., py::arg("compression")=1. )
    .def( py::init<Space &, std::shared_ptr<StimulusOccupancy>>(), py::arg("event_space"), py::arg("stimulus") )
    .def( py::init<std::shared_ptr<StimulusOccupancy>>(), py::arg("stimulus") )
    
    .def_property_readonly("changed", &PoissonLikelihood::changed,
    R"pbdoc(True of underlying distributions have changed and updated pre-computation is needed.)pbdoc")
    
    .def_property_readonly("mu", &PoissonLikelihood::mu,
    R"pbdoc(Mean event rate.)pbdoc")
    
    .def_property_readonly("ndim", &PoissonLikelihood::ndim,
    R"pbdoc(Combined dimensionality of event and stimulus space.)pbdoc")
    
    .def_property_readonly("ndim_stimulus", &PoissonLikelihood::ndim_stimulus,
    R"pbdoc(Dimensionality of stimulus space.)pbdoc")
    
    .def_property_readonly("ndim_events", &PoissonLikelihood::ndim_events,
    R"pbdoc(Dimensionality of event space.)pbdoc")
    
    .def_property_readonly("grid", py::cpp_function(&PoissonLikelihood::grid, py::return_value_policy::reference_internal),
    R"pbdoc(Evaluation grid in stimulus space.)pbdoc")
    
    .def_property_readonly("event_distribution", py::cpp_function(&PoissonLikelihood::event_distribution, py::return_value_policy::reference_internal),
    R"pbdoc(Underlying (compressed) density of merged events.)pbdoc")
    
    .def_property("random_insertion", &PoissonLikelihood::random_insertion, &PoissonLikelihood::set_random_insertion,
    R"pbdoc(Randomize new samples before merging into distribution.)pbdoc")
    
    //.def_property("rate_offset", &PoissonLikelihood::rate_offset, &PoissonLikelihood::set_rate_offset )
    .def_property("rate_scale", &PoissonLikelihood::rate_scale, &PoissonLikelihood::set_rate_scale,
    R"pbdoc(Event rate scaling factor that is applied during likelihood evaluation.)pbdoc")
    
    .def("to_yaml", [](PoissonLikelihood &m, bool b)->std::string {
        YAML::Emitter out;
        YAML::Node node = m.to_yaml(b);
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;
    }, py::arg("save_stimulus")=true,
    R"pbdoc(
        Represent stimulus occupancy as YAML.
        
        Parameters
        ----------
        save_stimulus : bool
            Convert stimulus occupancy distribution to YAML in addition
            to event distribution.
            
        Returns
        -------
        string
        
    )pbdoc")
    
    .def("save_to_yaml", [](const PoissonLikelihood& obj, std::string path, bool b) { obj.save_to_yaml( path, b ); },
    py::arg("path"), py::arg("save_stimulus")=true,
    R"pbdoc(
        Save Poisson likelihood to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        save_stimulus : bool
            Save stimulus occupancy distribution to YAML file in addition
            to event distribution.
        
    )pbdoc")
    
    .def("save_to_hdf5", &PoissonLikelihood::save_to_hdf5,
    py::arg("filename"), py::arg("save_stimulus")=true, py::arg("flags")=19, py::arg("path")="")
    
    .def_static("load_from_hdf5", [](std::string filename, std::string path, std::shared_ptr<StimulusOccupancy> stimulus) { return std::shared_ptr<PoissonLikelihood>( PoissonLikelihood::load_from_hdf5(filename, path, stimulus) ); },
    py::arg("filename"), py::arg("path")="", py::arg("stimulus")=nullptr,
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
    
    .def("add_events", [](PoissonLikelihood & obj, py::array_t<value, py::array::c_style | py::array::forcecast> events, unsigned int repetitions) {
        
        unsigned int ndim = obj.ndim();
        unsigned int nsamples;
        
        auto buf = events.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        obj.add_events( (value *) buf.ptr, nsamples, repetitions );
        
        },
    py::arg("events"), py::arg("repetitions")=1,
    R"pbdoc(
        Merge new events into event distribution.
        
        Parameters
        ----------
        events : (n,ndim) array
            Array of event data
        repetitions : int
            Number of repetitions for events to be merged.
        
    )pbdoc" )
        
    .def("precompute", &PoissonLikelihood::precompute, 
    R"pbdoc(Execute and cache intermediate computations.)pbdoc")
    
    .def_property_readonly("stimulus_logp", [](const PoissonLikelihood & obj)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(obj.grid().ndim(), sizeof(value));
        for (int k=obj.grid().ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * obj.grid().shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            const_cast<value*>(obj.stimulus_logp().data()),
            sizeof(value),
            py::format_descriptor<value>::value,
            obj.grid().ndim(),
            obj.grid().shape(),
            strides
        ));
        
        return result;
        
        },
    R"pbdoc(Log probability of stimulus distribution evaluated on grid.)pbdoc")
    
    .def_property_readonly("event_rate", [](const PoissonLikelihood & obj)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(obj.grid().ndim(), sizeof(value));
        for (int k=obj.grid().ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * obj.grid().shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            const_cast<value*>(obj.event_rate().data()),
            sizeof(value),
            py::format_descriptor<value>::value,
            obj.grid().ndim(),
            obj.grid().shape(),
            strides
        ));
        
        return result;
        
        },
    R"pbdoc(Marginal event rate evaluated on grid.)pbdoc")
    
    //.def_property_readonly("offset", [](const PoissonLikelihood & obj)->py::array_t<value> {
        
        //std::vector<long unsigned int> strides(obj.grid().ndim(), sizeof(value));
        //for (int k=obj.grid().ndim()-2; k>=0; --k) {
            //strides[k] = strides[k+1] * obj.grid().shape()[k+1];
        //}
        
        //// create output buffer
        //auto result = py::array( py::buffer_info(
            //const_cast<value*>(obj.offset().data()),
            //sizeof(value),
            //py::format_descriptor<value>::value,
            //obj.grid().ndim(),
            //obj.grid().shape(),
            //strides
        //));
        
        //return result;
        
        //} )
    
    .def("logL", [](PoissonLikelihood & obj, py::array_t<value, py::array::c_style | py::array::forcecast> events, value delta_t)->py::array_t<value> {
        
        unsigned int ndim = obj.ndim_events();
        unsigned int nsamples;
        
        auto buf = events.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
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
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + obj.grid().size(), 0. );
        
        obj.logL( (value *) buf.ptr, nsamples, delta_t, (value *) result_buf.ptr );
        
        return result;
        
        },
    py::arg("events"), py::arg("delta"),
    R"pbdoc(
        Evaluate log likelihood on grid given observed events.
        
        Parameters
        ----------
        events : (n,ndim) array
            Array with event data.
        delta : float
            Time duration over which events were observed.
        
        Returns
        -------
        nd array
        
    )pbdoc")
    
    .def("likelihood", [](PoissonLikelihood & obj, py::array_t<value, py::array::c_style | py::array::forcecast> events, value delta_t)->py::array_t<value> {
        
        unsigned int ndim = obj.ndim_events();
        unsigned int nsamples;
        
        auto buf = events.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
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
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + obj.grid().size(), 0. );
        
        obj.likelihood( (value *) buf.ptr, nsamples, delta_t, (value *) result_buf.ptr );
        
        return result;
        
        },
    py::arg("events"), py::arg("delta"),
    R"pbdoc(
        Evaluate likelihood on grid given observed events.
        
        Parameters
        ----------
        events : (n,ndim) array
            Array with event data.
        delta : float
            Time duration over which events where observed.
        
        Returns
        -------
        nd array
        
    )pbdoc")
    
    .def("event_prob", [](PoissonLikelihood & obj, py::array_t<value, py::array::c_style | py::array::forcecast> events)->py::array_t<value> {
        
        unsigned int ndim = obj.ndim_events();
        unsigned int nsamples;
        
        auto buf = events.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
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
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + obj.grid().size(), 0. );
        
        obj.event_prob( (value *) buf.ptr, nsamples, (value *) result_buf.ptr );
        
        return result;
        
        },
    py::arg("events"),
    R"pbdoc(
        Probability of observing events evaluated on grid.
        
        Parameters
        ----------
        events : (n,ndim) array
            Array with event data
        
        Returns
        -------
        nd array
        
    )pbdoc")
    
    .def("event_logp", [](PoissonLikelihood & obj, py::array_t<value, py::array::c_style | py::array::forcecast> events)->py::array_t<value> {
        
        unsigned int ndim = obj.ndim_events();
        unsigned int nsamples;
        
        auto buf = events.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
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
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + obj.grid().size(), 0. );
        
        obj.event_logp( (value *) buf.ptr, nsamples, (value *) result_buf.ptr );
        
        return result;
        
        },
    py::arg("events"),
    R"pbdoc(
        Log probability of observing events evaluated on grid.
        
        Parameters
        ----------
        events : (n,ndim) array
            Array with event data
        
        Returns
        -------
        nd array
        
    )pbdoc");

}
