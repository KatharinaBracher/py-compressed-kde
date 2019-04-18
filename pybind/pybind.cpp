#include "pybind.hpp"

#include <highfive/H5File.hpp>

void pybind_component( py::module & );
void pybind_kernel( py::module & );
void pybind_grid( py::module & );
void pybind_space( py::module & );
void pybind_mixture( py::module & );
void pybind_stimulus( py::module & );
void pybind_likelihood( py::module & );
void pybind_decoder( py::module & );

PYBIND11_MODULE(compressed_kde, m) {
    py::options options;
    options.disable_function_signatures();

    m.doc() = \
    R"pbdoc(
        ===================================================
        Compressed KDE (:mod:`fklab.decode.compressed_kde`)
        ===================================================
        
        .. currentmodule:: fklab.decode.compressed_kde
        
        Classes for compressed kernel density estimation and decoding.
        
        .. autosummary::
            :toctree: generated/
    
            GaussianKernel
            EpanechnikovKernel
            BoxKernel
            EuclideanSpace
            CategoricalSpace
            CircularSpace
            EncodedSpace
            MultiSpace
            Mixture
            PartialMixture
            Stimulus
            PoissonLikelihood
            Decoder
            
    )pbdoc";
    
    py::enum_<Flags>(m,"Flags",py::arithmetic())
    .value("ReadOnly", Flags::ReadOnly)
    .value("ReadWrite", Flags::ReadWrite)
    .value("Truncate", Flags::Truncate)
    .value("Excl", Flags::Excl)
    .value("Debug", Flags::Debug)
    .value("Create", Flags::Create)
    .value("Overwrite", Flags::Overwrite)
    .value("OpenOrCreate", Flags::OpenOrCreate)
    .export_values();

    pybind_component(m);
    pybind_kernel(m);
    pybind_grid(m);
    pybind_space(m);
    pybind_mixture(m);
    pybind_stimulus(m);
    pybind_likelihood(m);
    pybind_decoder(m);
    
}
