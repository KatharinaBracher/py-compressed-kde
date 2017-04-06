#include "kernel_vonmises.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include <cmath>

value vonmises_scale_factor( value kappa, bool log ) {
    
    value s;
    
    if (kappa>KAPPA_GAUSS_APPROX) {
        if (log) {
            s = -0.5 * fastlog( 2*M_PI/kappa );
        } else {
            s = std::pow( 2*M_PI/kappa, -0.5 );
        }
    } else {
        if (log) {
            s = -0.5 * fastlog( 2. * M_PI * boost::math::cyl_bessel_i( 0., kappa )); 
        } else {
            s = 1. / (2. * M_PI * boost::math::cyl_bessel_i( 0., kappa )); 
        }
    }
    
    return s;
    
}

