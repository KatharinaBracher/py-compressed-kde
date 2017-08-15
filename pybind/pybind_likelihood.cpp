#include "pybind.hpp"

#include "likelihood.hpp"

void pybind_likelihood(py::module &m) {
    
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
    
    .def("to_yaml", [](StimulusOccupancy &m)->std::string {
        YAML::Emitter out;
        YAML::Node node = m.toYAML();
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
    
    .def("save", [](const StimulusOccupancy& obj, std::string path) { obj.save( path ); }, py::arg("path"),
    R"pbdoc(
        Save stimulus occupancy to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        
    )pbdoc")
    
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
        YAML::Node node = m.toYAML(b);
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
    
    .def("save", [](const PoissonLikelihood& obj, std::string path, bool b) { obj.save( path, b ); },
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
    
    .def_property_readonly("grid_sizes", &Decoder::grid_sizes,
    R"pbdoc(Grid size for each stimulus space.)pbdoc")
    
    .def("grid_size", &Decoder::grid_size, py::arg("source"),
    R"pbdoc(
        Grid size.
        
        Parameters
        ----------
        source : int
            Index of source likelihood (zero-based).
        
        Returns
        -------
        int : grid size
            
    )pbdoc")
    
    .def_property_readonly("nenabled_sources", &Decoder::nenabled_sources,
    R"pbdoc(Number of sources that are used for decoding.)pbdoc")
    
    .def_property_readonly("enabled_sources", &Decoder::enabled_sources,
    R"pbdoc(Enabled state for all sources.)pbdoc")
    
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
    
    .def("decode", [](Decoder & obj, std::vector<py::array_t<value, py::array::c_style | py::array::forcecast>> events, value delta_t, bool normalize)->std::vector<py::array_t<value>> {
        
        std::vector<py::array_t<value>> out;
        std::vector<value*> out_ptr;
        
        // for each union
        for (unsigned int k = 0 ; k<obj.n_union(); ++k) {
            
            // construct array buffer
            auto result = py::array( py::buffer_info(
                nullptr,
                sizeof(value),
                py::format_descriptor<value>::value,
                1,
                {obj.grid_size(k)},
                {sizeof(value)}
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
        
        // construct array buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            1,
            {obj.grid_size(index)},
            {sizeof(value)}
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
