#pragma once

#include "common.hpp"

#include <string>
#include <stdexcept>
#include <vector>

#include "yaml-cpp/yaml.h"

enum class KernelType { Gaussian, Epanechnikov, Box };

std::string kerneltype_tostring( KernelType k );
KernelType kerneltype_fromstring( std::string k );

class Kernel {
public:
    Kernel( KernelType k );
    
    virtual Kernel * clone() const { throw std::runtime_error("not implemented"); }
    
    KernelType type() const;
    
    virtual std::string to_string() const { return kerneltype_tostring(type()); }
    
    YAML::Node toYAML() const;
    virtual YAML::Node asYAML() const;
    
    virtual value scale_factor( unsigned int n, value * bw, bool log ) const;
    virtual value scale_factor( unsigned int n, value * bw, bool log, std::vector<bool>::const_iterator selection ) const;
    virtual value probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value probability( value dsquared ) const;
    virtual value log_probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value log_probability( value dsquared ) const;
    virtual value partial_logp( unsigned int n, const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection) const;
    
protected:
    KernelType type_;
};

template <typename T>
class KernelBase : public Kernel {
public:
    
    using Kernel::Kernel;

    virtual Kernel* clone() const override {
        return new T(static_cast<T const&>(*this));
    }
};    

