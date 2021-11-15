// ---------------------------------------------------------------------
// This file is part of the compressed decoder library.
//
// Copyright (C) 2020 - now Neuro-Electronics Research Flanders
//
// The compressed decoder library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The compressed decoder library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with falcon-core. If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------
#include "kernel_vonmises.hpp"
//#include <boost/math/special_functions/bessel.hpp>
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
            s = -0.5 * fastlog( 2. * M_PI *std::cyl_bessel_i( 0., kappa )); 
        } else {
            s = 1. / (2. * M_PI * std::cyl_bessel_i( 0., kappa )); 
        }
    }
    
    return s;
    
}

