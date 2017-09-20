#include "space_base.hpp"

#include "space_categorical.hpp"
#include "space_euclidean.hpp"
#include "space_multi.hpp"
#include "space_circular.hpp"
#include "space_encoded.hpp"

// yaml
std::unique_ptr<Space> space_from_yaml( const YAML::Node & node );
std::unique_ptr<Space> load_space_from_yaml( std::string path );

// hdf5
std::unique_ptr<Space> space_from_hdf5(const HighFive::Group & group);
