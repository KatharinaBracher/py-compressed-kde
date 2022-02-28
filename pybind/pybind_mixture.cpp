
#include "pybind.hpp"

#include "mixture.hpp"

void pybind_mixture(py::module &m) {
    
    // MIXTURE CLASS
    py::class_<Mixture>(m, "Mixture",
    R"pbdoc(
        Mixture class for (compressed) kernel density estimation.
        
        Parameters
        ----------
        space : Space
            Description of mixture space.
        threshold : scalar
            Compression threshold
            
    )pbdoc")
    .def(py::init<const Space&, value>(), py::arg("space"), py::arg("threshold")=THRESHOLD)
    
    .def_property_readonly("space", &Mixture::space, py::return_value_policy::reference_internal,
    R"pbdoc(Mixture space.)pbdoc")
    .def_property_readonly("sum_of_weights", &Mixture::sum_of_weights,
    R"pbdoc(Sum of weights of all samples that were added to the density.)pbdoc")
    .def_property_readonly("sum_of_nsamples", &Mixture::sum_of_nsamples,
    R"pbdoc(Number of samples that were added to the density.)pbdoc")
    .def_property("threshold", &Mixture::threshold, &Mixture::set_threshold,
    R"pbdoc(Compression threshold.)pbdoc")
    .def_property_readonly("ncomponents", &Mixture::ncomponents,
    R"pbdoc(Number of components in (compressed) density.)pbdoc")
    .def_property_readonly("weights", &Mixture::weights,
    R"pbdoc(Weights of all components.)pbdoc")
    .def_property_readonly("kernel_locations", [](const Mixture& obj) -> py::array_t<value> {
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            2,
            { obj.ncomponents(), obj.space().ndim() },
            { obj.space().ndim() * sizeof(value), sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        value* ptr = (value*) result_buf.ptr;
        
        for (auto const & c : obj.components()) {
            std::copy( c->location.begin(), c->location.end(), ptr );
            ptr += obj.space().ndim();
        }
        
        return result;
        
    },
    R"pbdoc(Locations of all components.)pbdoc")
    .def_property_readonly("kernel_bandwidths", [](const Mixture& obj) -> py::array_t<value> {
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            2,
            { obj.ncomponents(), obj.space().nbw() },
            { obj.space().nbw() * sizeof(value), sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        value* ptr = (value*) result_buf.ptr;
        
        for (auto const & c : obj.components()) {
            std::copy( c->bandwidth.begin(), c->bandwidth.end(), ptr );
            ptr += obj.space().nbw();
        }
        
        return result;
        
    },
    R"pbdoc(Bandwidths of all components.)pbdoc")
    .def_property_readonly("kernel_scale_factors", [](const Mixture& obj) -> py::array_t<value> {
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            1,
            { obj.ncomponents() },
            { sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        value* ptr = (value*) result_buf.ptr;
        
        for (auto const & c : obj.components()) {
            *ptr++ = c->scale_factor;
        }
        
        return result;
        
    },
    R"pbdoc(Scale factors of all components.)pbdoc")
    .def("clear", &Mixture::clear, R"pbdoc(Remove all components and reset mixture.)pbdoc")
    
    .def("save_to_hdf5", &Mixture::save_to_hdf5,
    py::arg("filename"), py::arg("flags")=Flags::OpenOrCreate|Flags::Truncate, py::arg("path")="",
    R"pbdoc(
        save_to_hdf5(filename, flags, path) -> None

        Save mixture to hdf5 file.

        Parameters
        ----------
        filename : str
            path to hdf5 file
        flags : int
            flags for file creation
        path : str
            path inside hdf5 file

    )pbdoc" )
    
    .def("save_to_yaml", [](const Mixture& obj, std::string path) { obj.save_to_yaml( path ); }, py::arg("path"),
    R"pbdoc(
        save_to_yaml(path) -> None
        
        Save mixture to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        
    )pbdoc")
    
    .def("to_yaml", [](Mixture &m)->std::string {
        YAML::Emitter out;
        YAML::Node node = m.to_yaml();
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;
    },
    R"pbdoc(
        to_yaml() -> str

        Represent mixture as YAML.
        
        Returns
        -------
        string
        
    )pbdoc")
    
    .def_static("from_yaml", [](std::string s)->std::unique_ptr<Mixture> {
        YAML::Node node = YAML::Load( s );
        return Mixture::from_yaml( node );
    }, py::arg("string"),
    R"pbdoc(
        from_yaml(str)-> Mixture

        Create mixture from YAML.
        
        Parameters
        ----------
        string : string
            YAML string space representation
        
        Returns
        -------
        Mixture
        
    )pbdoc")
    
    .def_static("load_from_yaml", [](std::string path) { return std::unique_ptr<Mixture>( Mixture::load_from_yaml(path) ); },
    py::arg("path"),
    R"pbdoc(
        load_from_yaml(path) -> Mixture

        Load mixture from file.
        
        Parameters
        ----------
        path : string
            path to YAML file
        
        Returns
        -------
        Mixture
        
    )pbdoc" )
    
    .def_static("load_from_hdf5", [](std::string path) { return std::unique_ptr<Mixture>( Mixture::load_from_hdf5(path) ); },
    py::arg("path"),
    R"pbdoc(
        load_from_hdf5(filename, path) -> Mixture

        Load mixture from hdf5 file.
        
        Parameters
        ----------
        filename : string
            path to hdf5 file
        path : string
            path inside hdf5 file
        
        Returns
        -------
        Mixture
        
    )pbdoc" )
    
    .def(py::pickle(
        &(pickle_get_state<Mixture>),
        &(pickle_set_state<Mixture, fb_serialize::Mixture>)
    ))

    .def("add", [](Mixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples) {
        
        unsigned int ndim = m.space().ndim();
        unsigned int nsamples;
        
        // check sample array
        auto buf = samples.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        m.add_samples( (value *) buf.ptr, nsamples );
        
    }, py::arg("samples"),
    R"pbdoc(
        add(samples) -> None

        Add samples to the mixture.
        
        New mixture components are added at the sample location with the
        default kernel bandwidth. No merging with existing components is
        performed.
        
        Parameters
        ----------
        samples : (n,ndim) array
            Array of samples
        
    )pbdoc")
    
    .def("merge", [](Mixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples, bool random=true) {
        
        unsigned int ndim = m.space().ndim();
        unsigned int nsamples;
        
        // check sample array
        auto buf = samples.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        m.merge_samples( (value *) buf.ptr, nsamples, random );
        
    }, py::arg("samples"), py::arg("random")=true,
    R"pbdoc(
        merge(samples, random) -> None

        Merge samples into the mixture.
        
        New mixture components are added at the sample location with the
        default kernel bandwidth and merged with existing components if 
        the mahalanobis distance is below the threshold.
        
        Parameters
        ----------
        samples : (n,ndim) array
            Array of samples
        random : bool
            Whether the samples will be randomized before merging.
        
    )pbdoc")
    
    .def("evaluate", [](Mixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples)->py::array_t<value> {
        
        unsigned int ndim = m.space().ndim();
        unsigned int nsamples;
        
        // check sample array
        auto buf = samples.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            1,
            { nsamples },
            { sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        //std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples, 0. );
        
        m.evaluate( (value *) buf.ptr, nsamples, (value *) result_buf.ptr );
        
        return result;
        
    }, py::arg("samples"))
    
    .def("evaluate", [](Mixture &m, Grid & grid)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(grid.ndim(), sizeof(value));
        for (int k=grid.ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * grid.shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            grid.ndim(),
            grid.shape(),
            strides
        ));
        
        auto result_buf = result.request();
        
        //std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples, 0. );
        
        m.evaluate( grid, (value *) result_buf.ptr );
        
        return result;
        
    }, py::arg("grid"),
    R"pbdoc(
        evaulate(*args,**kwargs) -> array

        Evaluate mixture at samples.

        .. py:function:: evaluate(samples)
                         evaluate(grid)

        Parameters
        ----------
        samples : (n,ndim) array
            Array of samples
        grid : Grid
            grid specification
        
        Returns
        -------
        ndarray
            Evaluated probabilities at sample locations
        
    )pbdoc" )
    
    .def("partial", [](Mixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples, py::array_t<bool, py::array::c_style | py::array::forcecast> selection)->py::array_t<value> {
        
        auto vec_selection = numpy_array_to_vector( selection );
        
        if (vec_selection.size() != m.space().ndim()) {
            throw std::runtime_error("Invalid selection.");
        }
        
        unsigned int ndim = std::count( vec_selection.begin(), vec_selection.end(), true );
        unsigned int nsamples;
        
        // check sample array
        auto buf = samples.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            2,
            { m.ncomponents(), nsamples },
            { sizeof(value)*nsamples, sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + m.ncomponents() * nsamples, 0. );
        
        m.partial( (value *) buf.ptr, nsamples, vec_selection, (value *) result_buf.ptr );
        
        return result;
        
    }, py::arg("samples"), py::arg("selection"),
    R"pbdoc(
        Partially evaluate mixture at samples for select dimensions.
        
        Parameters
        ----------
        samples : (n,nselect) array
            Array of samples
        selection : (ndim,) boolean array
            Selected dimensions for partial evaluation.
        
        Returns
        -------
        (ncomp,nsamples) array
            Partial log probabilities
        
    )pbdoc" )
    
    .def("partialize", [](Mixture &m, Grid & grid)->PartialMixture* {
        
        return m.partial( grid );
        
    }, py::arg("grid"))
    
    .def("partialize", [](Mixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples, py::array_t<bool, py::array::c_style | py::array::forcecast> selection)->PartialMixture* {
        
        auto vec_selection = numpy_array_to_vector( selection );
        
        unsigned int ndim = std::count( vec_selection.begin(), vec_selection.end(), true );
        unsigned int nsamples;
        
        // check sample array
        auto buf = samples.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        return m.partial( (value *) buf.ptr, nsamples, vec_selection );
        
    }, py::arg("samples"), py::arg("selection"),
    R"pbdoc(
        partialize(*args, **kwargs) -> PartialMixture

        Partially evaluate mixture at grid points for select dimensions.

        .. py:function:: partialize(grid)
                         partialize(samples, selection)

        Parameters
        ----------
        grid : Grid
            Grid on subspace of mixture
        samples : (n,nselect) array
            Array of samples
        selection : (ndim,) boolean array
            Selected dimensions for partial evaluation.
        
        Returns
        -------
        PartialMixture
        
    )pbdoc")
    
    .def("marginal", [](Mixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples, py::array_t<bool, py::array::c_style | py::array::forcecast> selection)->py::array_t<value> {
        
        auto vec_selection = numpy_array_to_vector( selection );
        
        if (vec_selection.size() != m.space().ndim()) {
            throw std::runtime_error("Invalid selection.");
        }
        
        unsigned int ndim = std::count( vec_selection.begin(), vec_selection.end(), true );
        unsigned int nsamples;
        
        // check sample array
        auto buf = samples.request();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            1,
            { nsamples },
            { sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples, 0. );
        
        m.marginal( (value *) buf.ptr, nsamples, vec_selection, (value *) result_buf.ptr );
        
        return result;
        
    }, py::arg("samples"), py::arg("selection"))
    
    .def("marginal", [](Mixture &m, Grid & grid)->py::array_t<value> {
        
        std::vector<long unsigned int> strides(grid.ndim(), sizeof(value));
        for (int k=grid.ndim()-2; k>=0; --k) {
            strides[k] = strides[k+1] * grid.shape()[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            grid.ndim(),
            grid.shape(),
            strides
        ));
        
        auto result_buf = result.request();
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + grid.size(), 0. );
        
        m.marginal( grid, (value *) result_buf.ptr );
        
        return result;
        
    }, py::arg("grid"),
    R"pbdoc(
        marginal(*args, **kwargs) -> array

        Evaluate marginal at samples for select dimensions.

        .. py:function:: marginal(samples, selection)
                         marginal(grid)

        Parameters
        ----------
        samples : (n,nselect) array
            Array of samples
        selection : (ndim,) boolean array
            Selected dimensions for marginal evaluation.
        grid : Grid
            Grid on subspace of mixture
        
        Returns
        -------
        ndarray
            Marginal probabilities. 
        
    )pbdoc" );
    
    
    
    
    
    py::class_<PartialMixture>(m, "PartialMixture",
    R"pbdoc(Partially evaluated mixture class.)pbdoc")
    .def_property_readonly("ncomponents", &PartialMixture::ncomponents,
    R"pbdoc(Number of components in partially evaluated density.)pbdoc")
    .def_property_readonly("nsamples", &PartialMixture::nsamples,
    R"pbdoc(Number of partially evaluated samples.)pbdoc")
    .def_property_readonly("partial_shape", &PartialMixture::partial_shape,
    R"pbdoc(Array shape of partially evaluated samples.)pbdoc")
    
    .def("mixture", &PartialMixture::mixture, py::return_value_policy::reference_internal,
    R"pbdoc(Parent mixture.)pbdoc")
    
    .def("partial_logp", [](PartialMixture &m)->py::array_t<value> {
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            2,
            { m.ncomponents(), m.nsamples() },
            { sizeof(value)*m.nsamples(), sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        // copy values
        std::copy( m.partial_logp().cbegin(), m.partial_logp().cend(), (value*) result_buf.ptr );
        
        return result;
        
    },
    R"pbdoc(
        partial_logp() -> array

        Retrieve precomputed partial log probabilities.
        
        Returns
        -------
        (ncomponents,nsamples) array
        
    )pbdoc")
    
    .def("complete", [](PartialMixture &m, py::array_t<value, py::array::c_style | py::array::forcecast> samples )->py::array_t<value> {
        
        // check sample array
        auto buf = samples.request();
        unsigned int nsamples;
        unsigned int ndim = m.mixture().space().ndim() - m.partial_shape().size();
        
        if ( (buf.ndim==1 && ndim==1) || (buf.ndim==2 && buf.shape[1]==ndim) ) {
            nsamples = buf.shape[0];
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        std::vector<long unsigned int> shape = m.partial_shape();
        shape.insert( shape.begin(), nsamples );
        
        std::vector<long unsigned int> strides(shape.size(), sizeof(value));
        for (int k=shape.size()-2; k>=0; --k) {
            strides[k] = strides[k+1] * shape[k+1];
        }
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            shape.size(), //2
            shape, //{ nsamples, m.nsamples() },
            strides //{ sizeof(value) * m.nsamples(), sizeof(value) }
        ));
        
        auto result_buf = result.request();
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + nsamples*m.nsamples(), 0. );
        
        m.complete( (value*) buf.ptr, nsamples, (value*) result_buf.ptr );
        
        return result;
        
    }, py::arg("samples"),
    R"pbdoc(
        complete(samples) -> array

        Complete partial probabilities with remaining part of samples.
        
        Parameters
        ----------
        samples : (n,ndim) array
            Partial samples to complete evaluation of mixture
        
        Returns
        -------
        (n, ...) ndarray
            
    )pbdoc")
    
    .def("marginal", [](PartialMixture &m)->py::array_t<value> {
        
        std::vector<long unsigned int> shape = m.partial_shape();
        
        std::vector<long unsigned int> strides(shape.size(), sizeof(value));
        for (int k=shape.size()-2; k>=0; --k) {
            strides[k] = strides[k+1] * shape[k+1];
        }
        
        
        // create output buffer
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            shape.size(),
            shape,
            strides
        ));
        
        auto result_buf = result.request();
        
        std::fill( (value*) result_buf.ptr, ((value*) result_buf.ptr) + m.nsamples(), 0. );
        
        m.marginal( (value*) result_buf.ptr );
        
        return result;
        
    },
    R"pbdoc(
        marginal() -> array
        
        Compute marginal probabilities.
        
        Returns
        -------
        ndarray
        
    )pbdoc");
    
}
