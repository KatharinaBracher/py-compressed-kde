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
        ======================================
        Compressed KDE (:mod:`compressed_kde`)
        ======================================
        
        .. currentmodule:: compressed_kde
        
        Classes for compressed kernel density estimation.
        
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

    py::module subm = m.def_submodule("decode",
                                      R"pbdoc(
                                          =====================================
                                          Decode (:mod:`compressed_kde.decode`)
                                          =====================================

                                          .. currentmodule:: compressed_kde.decode

                                          Classes for decoding.

                                          .. autosummary::
                                              :toctree: generated/

                                              Stimulus
                                              PoissonLikelihood
                                              Decoder

                                      )pbdoc");

    pybind_stimulus(subm);
    pybind_likelihood(subm);
    pybind_decoder(subm);
    
}
