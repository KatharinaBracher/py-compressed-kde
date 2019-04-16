#include "pybind.hpp"

void pybind_component( py::module & );
void pybind_kernel( py::module & );
void pybind_grid( py::module & );
void pybind_space( py::module & );
void pybind_mixture( py::module & );
void pybind_stimulus( py::module & );
void pybind_likelihood( py::module & );
void pybind_decoder( py::module & );

PYBIND11_MODULE(compressed_kde, m) {

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
    
    pybind_component(m);
    pybind_kernel(m);
    pybind_grid(m);
    pybind_space(m);
    pybind_mixture(m);
    pybind_stimulus(m);
    pybind_likelihood(m);
    pybind_decoder(m);
    
}
