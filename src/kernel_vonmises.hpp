#pragma once

#include "common.hpp"

static const value KAPPA_GAUSS_APPROX = 40.;
value vonmises_scale_factor( value kappa, bool log );
