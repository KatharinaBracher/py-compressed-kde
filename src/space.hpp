#include "space_base.hpp"

#include "space_categorical.hpp"
#include "space_euclidean.hpp"
#include "space_multi.hpp"
#include "space_circular.hpp"
#include "space_encoded.hpp"


Space * space_from_YAML( const YAML::Node & node );

Space * load_space( std::string path );
